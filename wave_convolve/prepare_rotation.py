import os

import numpy as np


def compute_rtz_to_xyz_matrix(theta, phi):
    """
    RTZ 分量 → XYZ 坐标
    - 输入：theta, phi (rad)，theta 极角，phi 方位角
    - RTZ = [theta方向分量, phi方向分量, r方向分量]
    - 输出：3x3 矩阵 R，满足 vec_xyz = R @ vec_rtz
    """
    sin_th = np.sin(theta)
    cos_th = np.cos(theta)
    sin_ph = np.sin(phi)
    cos_ph = np.cos(phi)

    R = np.array([
        [cos_th * cos_ph, -sin_ph, sin_th * cos_ph],
        [cos_th * sin_ph,  cos_ph, sin_th * sin_ph],
        [-sin_th,           0.0,    cos_th]
    ])
    return R


def main():
    # 路径
    grid_path = '../wave_extrapolation/outputs'
    save_dir = './precomputed'
    os.makedirs(save_dir, exist_ok=True)

    # 读取网格
    phi = np.loadtxt(f'{grid_path}/grid_azim_ml_ready.txt')  # (Np,)
    theta = np.loadtxt(f'{grid_path}/grid_dist_ml_ready.txt')  # (Nt,)
    Nt = theta.shape[0]
    Np = phi.shape[0]

    print(f'Grid: theta={Nt}, phi={Np}')

    # 预分配：每个 (theta, phi) 一个 3x3 旋转矩阵
    rotations = np.zeros((Nt, Np, 3, 3), dtype=np.float64)

    for it in range(Nt):
        for ip in range(Np):
            R_mat = compute_rtz_to_xyz_matrix(theta[it], phi[ip])
            rotations[it, ip] = R_mat

    # 存储
    np.save(f'{save_dir}/rotation_matrices.npy', rotations)
    print(f'✅ 旋转矩阵存储完成：{save_dir}/rotation_matrices.npy')


if __name__ == '__main__':
    main()
