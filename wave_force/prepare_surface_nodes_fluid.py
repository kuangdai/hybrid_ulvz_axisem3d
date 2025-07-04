import json
import numpy as np
import os

from geodetic import GeoPoints


def main():
    """
    预处理 Fluid 严格表面节点，避免重复，保存唯一节点列表
    输出：
        surface_nodes_cache/index_fluid.txt
        surface_nodes_cache/reciprocal_stations_fluid.txt
    """

    grid_path = "../wave_extrapolation/outputs"
    save_dir = "./surface_nodes_cache/"
    os.makedirs(save_dir, exist_ok=True)

    # 读取网格与参数
    phi = np.loadtxt(f"{grid_path}/grid_azim_ml_ready.txt")  # φ 方向
    theta = np.loadtxt(f"{grid_path}/grid_dist_ml_ready.txt")  # θ 方向
    depth = np.loadtxt(f"{grid_path}/grid_depth_ml_ready_fluid.txt")  # Fluid 深度
    with open(f"{grid_path}/args.json", "r") as fs:
        meta = json.load(fs)

    Np, Nt, Nd = phi.shape[0], theta.shape[0], depth.shape[0]

    node_list = []

    # 1. 最深 depth 层
    node_list.extend((Nd - 1, it, ip) for ip in range(Np) for it in range(Nt))

    # 2. 最大 theta 层，排除重叠
    node_list.extend((id, Nt - 1, ip) for ip in range(Np) for id in range(Nd - 1))

    node_array = np.array(node_list, dtype=int)
    print(f"✅ 识别 Fluid 严格表面节点共 {node_array.shape[0]} 个")

    np.savetxt(f"{save_dir}/index_fluid.txt", node_array, fmt="%d")

    # 计算坐标
    r_vals = 6371.0 - depth[node_array[:, 0]]
    t_vals = theta[node_array[:, 1]]
    p_vals = phi[node_array[:, 2]]
    rtp = np.stack([r_vals, t_vals, p_vals], axis=1)

    gp = GeoPoints.create_rtp_src_centered(rtp, src_lat=meta["ulvz_lat"], src_lon=meta["ulvz_lon"])

    # 输出文件
    lines = [
        f"{d}_{t}_{p} FLUID {lat:.8f} {lon:.8f} dummy {dep * 1e3:.6f}"
        for idx, (d, t, p), (lat, lon, dep) in zip(
            range(len(node_array)), node_array, gp.lld
        )
    ]

    with open(f"{save_dir}/reciprocal_stations_fluid.txt", "w") as fs:
        fs.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
