import numpy as np
import torch

from brick import SolidElement, FluidElement


def load_rotation_matrices():
    rotations = np.load(f"precomputed/rotation_matrices.npy")  # (Ntheta, Nphi, 3, 3)
    return torch.tensor(rotations, dtype=torch.float32)


def load_elements(medium, device='cpu'):
    connectivity = np.loadtxt(f'precomputed/connectivity_{medium}.txt', dtype=int)
    connectivity = connectivity.reshape(connectivity.shape[0], 8, -1)
    ElementClass = SolidElement if medium == "solid" else FluidElement
    elements = []
    for i in range(len(connectivity)):
        element = ElementClass.load_from(f"precomputed/{medium}_elements/{medium}_{i}.pt",
                                         device=device)
        elements.append(element)
    return elements, torch.tensor(connectivity).to(device)


def load_reciprocal(u_near, station_name, medium):
    data = np.load(f"../wave_reciprocal/reciprocal_simulations/"
                   f"{station_name}/output/rotated_{medium}.npz")["arr_0"]
    nodes = np.loadtxt(f"precomputed/reciprocal_nodes_{medium}.txt", dtype=int)
    full_data = np.zeros((len(data), *u_near.shape[1:5]))
    for i, (t, p, d) in enumerate(nodes):
        full_data[:, :, t, p, d] = data[:, :, i]
    return torch.tensor(full_data, dtype=torch.float32, device=u_near.device)

def gather_nodes(u, connectivity):
    """
    u: [T, C, D0, D1, D2]  张量场，T时刻，C通道，三维网格
    connectivity: [N, 8, 3]  N个单元，每单元8个点，点的3D索引
    返回: [N, T, C, 8]  每个单元提取的点数据
    """
    N, _, _ = connectivity.shape
    T, C, D0, D1, D2 = u.shape

    # 三维索引
    d0 = connectivity[:, :, 0]  # [N, 8]
    d1 = connectivity[:, :, 1]
    d2 = connectivity[:, :, 2]

    # 计算一维flatten索引
    flat_idx = d0 * (D1 * D2) + d1 * D2 + d2  # [N, 8]

    # 展平u
    u_flat = u.reshape(T, C, -1)  # [T, C, D0*D1*D2]

    # 扩展索引
    flat_idx = flat_idx.unsqueeze(1).unsqueeze(2)  # [N, 1, 1, 8]
    u_flat = u_flat.unsqueeze(0).expand(N, T, C, -1)  # [N, T, C, D0*D1*D2]

    # gather操作
    output = torch.gather(u_flat, -1, flat_idx.expand(N, T, C, 8))  # [N, T, C, 8]
    output = output.permute(0, 1, 3, 2)  # [N, T, 8, C]
    return output

