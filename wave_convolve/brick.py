import numpy as np
import torch


def _shape_function(xi, eta, zeta):
    """
    计算8节点线性单元标准形函数值
    """
    N = torch.zeros(8)
    signs = torch.tensor([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                          [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
    coeff = 0.125
    for i in range(8):
        sx, sy, sz = signs[i]
        N[i] = coeff * (1 + sx * xi) * (1 + sy * eta) * (1 + sz * zeta)
    return N


def _shape_function_derivatives(xi, eta, zeta):
    """
    计算形函数对局部坐标偏导
    """
    dN = torch.zeros((8, 3))
    coeff = 0.125
    signs = torch.tensor([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                          [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
    for i in range(8):
        sx, sy, sz = signs[i]
        dN[i, 0] = coeff * sx * (1 + sy * eta) * (1 + sz * zeta)
        dN[i, 1] = coeff * (1 + sx * xi) * sy * (1 + sz * zeta)
        dN[i, 2] = coeff * (1 + sx * xi) * (1 + sy * eta) * sz
    return dN


def _build_B_matrix(dN_dx):
    """
    构造6x24标准B矩阵，严格符合小应变张量Voigt定义：

    Voigt顺序：
    0 - ε_xx
    1 - ε_yy
    2 - ε_zz
    3 - γ_yz = 2 ε_yz
    4 - γ_xz = 2 ε_xz
    5 - γ_xy = 2 ε_xy
    """
    B = torch.zeros((6, 24))
    for i in range(8):
        idx = i * 3
        dNxi, dNyi, dNzi = dN_dx[i]

        # 正确的主对角线分量
        B[0, idx] = dNxi  # ∂u_x / ∂x
        B[1, idx + 1] = dNyi  # ∂u_y / ∂y
        B[2, idx + 2] = dNzi  # ∂u_z / ∂z

        # 对称剪切分量
        B[3, idx + 1] = 2 * dNzi  # ∂u_y / ∂z
        B[3, idx + 2] = 2 * dNyi  # ∂u_z / ∂y

        B[4, idx] = 2 * dNzi  # ∂u_x / ∂z
        B[4, idx + 2] = 2 * dNxi  # ∂u_z / ∂x

        B[5, idx] = 2 * dNyi  # ∂u_x / ∂y
        B[5, idx + 1] = 2 * dNxi  # ∂u_y / ∂x
    return B


def fft_convolve_multidim(a, b, sum_dim=True):
    """
    多通道一维卷积，支持 [time, d] 结构，d为通道数，自动零填充至2的次幂长度，提升FFT性能

    参数：
        a: [L1, d]  时间、通道，代表 force 或 displacement
        b: [L2, d]  时间、通道，代表 reciprocal green

    返回：
        result:
            若 sum_dim=True : [L1 + L2 - 1]，所有通道卷积并求和
            若 sum_dim=False: [L1 + L2 - 1, d]，各通道独立卷积
    """
    L1, d = a.shape
    L2 = b.shape[0]
    N = L1 + L2 - 1  # 理论卷积长度
    n_fft = 2 ** int(np.ceil(np.log2(N)))  # 2的幂次，保障FFT高效

    # 执行FFT，使用norm='forward'确保严格能量一致性
    A = torch.fft.rfft(a, n=n_fft, dim=0, norm='forward')  # [fft_len, d]
    B = torch.fft.rfft(b, n=n_fft, dim=0, norm='forward')
    C = A * B  # 逐通道频域乘积

    result = torch.fft.irfft(C, n=n_fft, dim=0, norm='forward')[:N]  # [N, d]

    if sum_dim:
        result = torch.sum(result, dim=1)  # 所有通道累加，返回 [N]

    return result


def _compute_stress2traction(face_node_pos, gp_index):
    """
    计算给定面上第gp_index个高斯点的单位法向及对应的应力-牵引力投影矩阵P

    P 矩阵作用：Voigt格式下应力映射至牵引力 t_i = σ_ij * n_j
    Voigt编码顺序：
        0 → σ_xx
        1 → σ_yy
        2 → σ_zz
        3 → σ_yz = σ_zy
        4 → σ_xz = σ_zx
        5 → σ_xy = σ_yx

    :param face_node_pos: (4, 3) 面上四个节点物理坐标
    :param gp_index: int，高斯点编号 0~3
    :return: (3, 6) 投影矩阵 P
    """
    # 计算单位法向
    point_self = face_node_pos[gp_index]
    point_next = face_node_pos[(gp_index + 1) % 4]
    point_prev = face_node_pos[gp_index - 1]
    vec1 = point_self - point_prev
    vec2 = point_next - point_self
    # 计算面单元当前高斯点贡献的带面积法向
    # cross结果为平行四边形面积，乘0.5对应四边形的四分之一，可直接用于积分权重
    normal = 0.5 * torch.cross(vec1, vec2)

    # 构造 P矩阵： σ → t 转换，以满足Voigt格式
    n0, n1, n2 = normal[0], normal[1], normal[2]
    P = torch.zeros((3, 6), dtype=torch.float32)

    # t_x = σ_xx * n_x + σ_xy * n_y + σ_xz * n_z
    P[0, 0] = n0  # σ_xx n_x
    P[0, 5] = n1  # σ_xy n_y
    P[0, 4] = n2  # σ_xz n_z

    # t_y = σ_yx * n_x + σ_yy * n_y + σ_yz * n_z
    P[1, 5] = n0  # σ_xy n_x
    P[1, 1] = n1  # σ_yy n_y
    P[1, 3] = n2  # σ_yz n_z

    # t_z = σ_zx * n_x + σ_zy * n_y + σ_zz * n_z
    P[2, 4] = n0  # σ_xz n_x
    P[2, 3] = n1  # σ_yz n_y
    P[2, 2] = n2  # σ_zz n_z

    return P


class SolidElement:
    def __init__(self, node_pos, lambda_g, mu_g, gamma_face_index, device="cpu"):
        """
        初始化SolidElement
        :param node_pos: (8, 3) 节点坐标
        :param lambda_g: (8,) Lamé 参数 λ
        :param mu_g: (8,) 剪切模量 μ
        :param gamma_face_index: 1~6 指定面编号
        """
        self.gamma_face_index = gamma_face_index
        self.device = device
        self.face_disp2stress = None  # [4, 6, 24]
        self.face_stress2traction = None  # [4, 3, 6]
        self.face_node2gauss = None  # [4, 8]
        self._precompute(node_pos, lambda_g, mu_g)

    def _precompute(self, node_pos, lambda_g, mu_g):
        """
        预计算面上所有矩阵
        """
        # 定义三维 8 个 Gauss 点位置
        sqrt3_inv = 1.0 / np.sqrt(3)
        gauss_1d = torch.tensor([-sqrt3_inv, sqrt3_inv])
        gauss_3d = torch.stack(torch.meshgrid(
            gauss_1d, gauss_1d, gauss_1d, indexing='ij'), dim=-1).reshape(-1, 3)

        # 定义面

        # 节点编号示意：
        #   7 -------- 6
        #  /|         /|
        # 4 -------- 5 |
        # | |        | |
        # | 3 -------| 2
        # |/         |/
        # 0 -------- 1
        face_node_dict = {
            1: [0, 4, 7, 3],  # x = -1
            2: [1, 2, 6, 5],  # x = +1
            3: [0, 1, 5, 4],  # y = -1
            4: [2, 3, 7, 6],  # y = +1
            5: [0, 3, 2, 1],  # z = -1
            6: [4, 5, 6, 7]  # z = +1
        }
        face_dim_direction_dict = {
            1: [0, -1],  # [dim, face position]
            2: [0, 1],
            3: [1, -1],
            4: [1, 1],
            5: [2, -1],
            6: [2, 1],
        }

        face_node_idx = face_node_dict[self.gamma_face_index]
        dim, pos = face_dim_direction_dict[self.gamma_face_index]

        # 面上节点坐标
        node_pos = torch.tensor(node_pos, dtype=torch.float32)
        node_pos_face = node_pos[face_node_idx]

        # 面上Gauss点坐标
        gauss_pos_face = gauss_3d[face_node_idx]  # 快速从三维获取
        gauss_pos_face[:, dim] = pos  # 将三维投影到面上

        # 对Gauss点循环
        face_disp2stress, face_stress2traction, face_node2gauss = [], [], []
        for gp in range(4):
            # 形函数，将节点位移map到面上Gauss点位移
            xi, eta, zeta = gauss_pos_face[gp]
            N = _shape_function(xi, eta, zeta)
            face_node2gauss.append(N)

            # 应变矩阵
            dN_dxi = _shape_function_derivatives(xi, eta, zeta)
            J = dN_dxi.T @ node_pos
            invJ = torch.inverse(J)
            dN_dx = dN_dxi @ invJ.T
            B = _build_B_matrix(dN_dx)

            # 刚度矩阵
            lam, mu = lambda_g[gp], mu_g[gp]
            D = torch.zeros((6, 6))
            D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
            D[3, 3] = D[4, 4] = D[5, 5] = mu
            D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = lam

            # 位移->应力矩阵
            face_disp2stress.append(D @ B)

            # 应力->牵引力矩阵
            face_stress2traction.append(_compute_stress2traction(node_pos_face, gp))

        self.face_disp2stress = torch.stack(face_disp2stress, dim=0).to(self.device)
        self.face_stress2traction = torch.stack(face_stress2traction, dim=0).to(self.device)
        self.face_node2gauss = torch.stack(face_node2gauss, dim=0).to(self.device)

    def _compute_traction(self, u):
        """
        计算面上牵引力
        :param u: [time, 8, 3]  多时刻、节点、方向的位移
        :return: traction [4, 3, time]  Gauss点、方向、时间
        """
        u = u.to(self.device)  # u: [time, 8, 3]
        u = u.reshape(u.shape[0], 24)  # u: [time, 24]

        # face_disp2stress: [4, 6, 24]
        # u: [time, 24]
        # stress: [time, 4, 6]
        stress = torch.einsum('gsn, tn -> tgs', self.face_disp2stress, u)

        # face_stress2traction: [4, 3, 6]
        # stress: [time, 4, 6]
        # traction: [time, 4, 3]
        traction = torch.einsum('gcs, tgs -> tgc', self.face_stress2traction, stress)

        return traction

    def _compute_disp_gauss(self, u):
        """
        计算面上高斯点的位移
        :param u: [time, 8, 3]  多时刻、节点、方向的位移
        :return: disp [time, 4, 3]  时间、Gauss点、方向
        """
        u = u.to(self.device)  # u: [time, 8, 3]

        # face_node2gauss: [4, 8]
        # u: [time, 8, 3]
        # disp: [time, 4, 3]
        disp = torch.einsum('gn, tnc -> tgc', self.face_node2gauss, u)

        return disp

    def compute_convolve(self, u_near, u_recip):
        t_near = self._compute_traction(u_near).reshape(-1, 12)  # [tn, 4 * 3]
        ug_near = self._compute_disp_gauss(u_near).reshape(-1, 12)  # [tn, 4 * 3]
        t_recip = self._compute_traction(u_recip).reshape(-1, 12)  # [tr, 4 * 3]
        ug_recip = self._compute_disp_gauss(u_recip).reshape(-1, 12)  # [tr, 4 * 3]
        first = fft_convolve_multidim(t_near, ug_recip)
        second = fft_convolve_multidim(ug_near, t_recip)
        convolved = first - second
        return convolved
