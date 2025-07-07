import json
import os
import sys

import numpy as np

from brick import SolidElement, FluidElement
from geodetic import GeoPoints
from prem import PREM


def build_solid_element(node_indices, face_index, name):
    node_xyz, lambda_n, mu_n = [], [], []

    for t_, p_, d_ in node_indices:
        r = R_earth - depths[d_]
        th = thetas[t_]
        ph = phis[p_]

        x = r * np.sin(th) * np.cos(ph)
        y = r * np.sin(th) * np.sin(ph)
        z = r * np.cos(th)
        node_xyz.append([x, y, z])

        lam, mu, rho = prem.query_solid(r)
        lambda_n.append(lam)
        mu_n.append(mu)

    element = SolidElement(node_xyz, lambda_n, mu_n, gamma_face_index=face_index)
    element.save_to(f"{save_dir}/{medium}_elements/{name}.pt")


def build_fluid_element(node_indices, face_index, name):
    node_xyz, rho_n = [], []

    for t_, p_, d_ in node_indices:
        r = R_earth - depths[d_]
        th = thetas[t_]
        ph = phis[p_]

        x = r * np.sin(th) * np.cos(ph)
        y = r * np.sin(th) * np.sin(ph)
        z = r * np.cos(th)
        node_xyz.append([x, y, z])

        kappa, rho = prem.query_fluid(r)
        rho_n.append(rho)

    element = FluidElement(node_xyz, rho_n, gamma_face_index=face_index)
    element.save_to(f"{save_dir}/{medium}_elements/{name}.pt")


if __name__ == "__main__":
    medium = sys.argv[1]
    build_element = build_solid_element if medium == "solid" else build_fluid_element

    # 路径设置
    grid_path = '../wave_extrapolation/outputs'
    save_dir = "./precomputed"
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(save_dir + f"/{medium}_elements", exist_ok=True)

    # 1D模型
    R_earth = 6371000.0  # 地球半径，单位 m
    prem = PREM('./prem_1s_mantle.csv') if medium == "solid" else PREM('./prem_1s_core.csv')

    # ULVZ 位置
    with open(f"{grid_path}/args.json", "r") as fs:
        meta = json.load(fs)
    ulvz_lat, ulvz_lon = meta["ulvz_lat"], meta["ulvz_lon"]

    # 读取网格
    phis = np.loadtxt(f'{grid_path}/grid_azim_ml_ready.txt')  # (Np,)
    thetas = np.loadtxt(f'{grid_path}/grid_dist_ml_ready.txt')  # (Nt,)
    depths = np.loadtxt(f'{grid_path}/grid_depth_ml_ready_{medium}.txt') * 1000  # m
    Np = phis.shape[0]
    Nt = thetas.shape[0]
    Nd = depths.shape[0]

    # Connectivity
    connectivity = []
    print(f'Grid: depth={Nd}, theta={Nt}, phi={Np}')

    ### 1. 最外层 depth 层单元 ###
    if medium == "solid":
        deeper = 1
        shallower = 0
        face_depth = 6
    elif medium == "fluid":
        deeper = Nd - 1
        shallower = Nd - 2
        face_depth = 5
    else:
        assert 0, "Invalid medium"
    for ip in range(Np):
        for it in range(Nt - 1):
            ip_next = (ip + 1) % Np

            indices = [
                (it, ip, deeper),
                (it + 1, ip, deeper),
                (it + 1, ip_next, deeper),
                (it, ip_next, deeper),
                (it, ip, shallower),
                (it + 1, ip, shallower),
                (it + 1, ip_next, shallower),
                (it, ip_next, shallower)
            ]
            build_element(indices, face_depth, f"{medium}_{len(connectivity)}")
            connectivity.append(indices)

    print(f'✅ 最外层depth层单元共 {Np} × {Nt - 1} 个完成')

    ### 2. 最大 theta 层单元 ###
    for ip in range(Np):
        for id_ in range(0, Nd - 1):
            ip_next = (ip + 1) % Np

            indices = [
                (Nt - 2, ip, id_ + 1),
                (Nt - 1, ip, id_ + 1),
                (Nt - 1, ip_next, id_ + 1),
                (Nt - 2, ip_next, id_ + 1),
                (Nt - 2, ip, id_),
                (Nt - 1, ip, id_),
                (Nt - 1, ip_next, id_),
                (Nt - 2, ip_next, id_)
            ]
            build_element(indices, 2, f"{medium}_{len(connectivity)}")
            connectivity.append(indices)

    print(f'✅ 最大theta层单元共 {Np} × {Nd - 1} 个完成')

    # 存储 Connectivity
    connectivity = np.array(connectivity)
    np.savetxt(f"{save_dir}/connectivity_{medium}.txt",
               connectivity.reshape(len(connectivity), -1), fmt="%d")

    # Unique points
    idx_unique = []
    idx_seen = set()
    for connect in connectivity:
        for point in connect:
            t, p, d = point
            key = (t, p, d)
            if key not in idx_seen:
                idx_unique.append(key)
                idx_seen.add(key)
    idx_unique = np.array(idx_unique)

    # 构造 GeoPoints
    r_vals = R_earth - depths[idx_unique[:, 2]]
    t_vals = thetas[idx_unique[:, 0]]
    p_vals = phis[idx_unique[:, 1]]
    rtp = np.stack([r_vals / 1000., t_vals, p_vals], axis=1)  # r in km
    gp = GeoPoints.create_rtp_src_centered(rtp, src_lat=ulvz_lat, src_lon=ulvz_lon)

    # Stations 字符串拼接
    lines = [
        f"{t}_{p}_{d} {medium.upper()} {lat:.8f} {lon:.8f} dummy {dep * 1e3:.6f}"
        for (t, p, d), (lat, lon, dep) in zip(
            idx_unique, gp.lld
        )
    ]

    with open(f"{save_dir}/reciprocal_stations_{medium}.txt", "w") as fs:
        fs.write("\n".join(lines) + "\n")
    np.savetxt(f"{save_dir}/reciprocal_nodes_{medium}.txt", idx_unique, fmt="%d")
