import os

import numpy as np

from brick import FluidElement
from prem import PREM


def build_fluid_element(node_indices, face_index, name, save_path):
    """
    构建并保存 FluidElement
    """
    nodes, rho_n = [], []

    for t, p, d in node_indices:
        r = R_earth - depths[d]
        th = thetas[t]
        ph = phis[p]

        x = r * np.sin(th) * np.cos(ph)
        y = r * np.sin(th) * np.sin(ph)
        z = r * np.cos(th)
        nodes.append([x, y, z])

        _, rho = prem.query_fluid(r)
        rho_n.append(rho)

    elem = FluidElement(nodes, rho_n, gamma_face_index=face_index)
    elem.save_to(f"{save_path}/{name}.pt")
    np.savetxt(f"{save_path}/{name}.connect.txt", np.array(node_indices), fmt="%d")


if __name__ == "__main__":
    # 路径配置
    grid_path = '../wave_extrapolation/outputs'
    prem_path = './prem_1s.csv'
    save_dir = "./element_cache"
    os.makedirs(save_dir, exist_ok=True)

    R_earth = 6371000.0  # m

    # 读取网格
    phis = np.loadtxt(f'{grid_path}/grid_azim_ml_ready.txt')  # (Np,)
    thetas = np.loadtxt(f'{grid_path}/grid_dist_ml_ready.txt')  # (Nt,)
    depths = np.loadtxt(f'{grid_path}/grid_depth_ml_ready_fluid.txt') * 1000  # m

    Np = phis.shape[0]
    Nt = thetas.shape[0]
    Nd = depths.shape[0]

    print(f'Grid: depth={Nd}, theta={Nt}, phi={Np}')
    prem = PREM(prem_path)

    ### 1. 最深 depth 层单元 ###
    for ip in range(Np):
        for it in range(Nt - 1):
            ip_next = (ip + 1) % Np

            indices = [
                (it, ip, -1),
                (it + 1, ip, -1),
                (it, ip_next, -1),
                (it + 1, ip_next, -1),
                (it, ip, -2),
                (it + 1, ip, -2),
                (it, ip_next, -2),
                (it + 1, ip_next, -2)
            ]
            build_fluid_element(indices, face_index=5, name=f"fluid_bottom_p{ip}_t{it}", save_path=save_dir)

    print(f'✅ Fluid最深depth层单元共 {Np} × {Nt - 1} 个完成')

    ### 2. 最大 theta 层单元 ###
    for ip in range(Np):
        for id_ in range(0, Nd - 2):
            ip_next = (ip + 1) % Np

            indices = [
                (-2, ip, id_),
                (-1, ip, id_),
                (-2, ip_next, id_),
                (-1, ip_next, id_),
                (-2, ip, id_ + 1),
                (-1, ip, id_ + 1),
                (-2, ip_next, id_ + 1),
                (-1, ip_next, id_ + 1)
            ]
            build_fluid_element(indices, face_index=2, name=f"fluid_side_p{ip}_t{id_}", save_path=save_dir)

    print(f'✅ Fluid最大theta层单元共 {Np} × {Nd - 2} 个完成')
