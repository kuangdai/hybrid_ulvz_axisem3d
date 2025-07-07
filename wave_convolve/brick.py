import math

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
        B[3, idx + 1] = dNzi  # ∂u_y / ∂z
        B[3, idx + 2] = dNyi  # ∂u_z / ∂y

        B[4, idx] = dNzi  # ∂u_x / ∂z
        B[4, idx + 2] = dNxi  # ∂u_z / ∂x

        B[5, idx] = dNyi  # ∂u_x / ∂y
        B[5, idx + 1] = dNxi  # ∂u_y / ∂x
    return B


def compute_normal(face_node_pos, gp_index, voigt_solid):
    """
    计算面上第 gp_index 个高斯点的面积法向与应力-牵引力投影矩阵（或法向量）

    适配两种模式：
    1. Solid（voigt_solid=True）:
        - 返回 [3, 6] 投影矩阵 P
        - 满足 Voigt 格式下：t_i = σ_ij * n_j
        - Voigt编码顺序：
            0 → σ_xx
            1 → σ_yy
            2 → σ_zz
            3 → σ_yz = σ_zy
            4 → σ_xz = σ_zx
            5 → σ_xy = σ_yx

    2. Fluid（voigt_solid=False）:
        - 返回 [1, 3] 面法向（已包含面积权重）
        - 用于 t = n ⋅ ∇χ / ρ，其中 t、n、∇χ均为矢量

    注意：
    - 输出的 normal 已包含1/4面面积系数（四边形总面积/4），便于Gauss积分直接累加
    - 正确法向方向由单元节点逆时针排序确定，面外侧为正方向

    参数：
        face_node_pos : (4, 3)
            面上4个节点的三维物理坐标
        gp_index : int
            Gauss点编号，取值 0 ~ 3
        voigt_solid : bool
            Solid模式返回 [3, 6] 投影矩阵，Fluid模式返回 [1, 3] 法向

    返回：
        torch.Tensor
            - Solid模式： [3, 6] 投影矩阵 P
            - Fluid模式： [1, 3] 面法向，含面积权重
    """
    # 计算面积法向（含系数 1/4）
    point_self = face_node_pos[gp_index]
    point_next = face_node_pos[(gp_index + 1) % 4]
    point_prev = face_node_pos[gp_index - 1]
    vec1 = point_self - point_prev
    vec2 = point_next - point_self
    normal = 0.25 * torch.cross(vec1, vec2, dim=-1)  # 四边形1/4面积法向

    if not voigt_solid:
        return normal.unsqueeze(0)  # (1, 3)，Fluid模式直接返回矢量

    # Solid模式：构造Voigt投影矩阵
    n0, n1, n2 = normal[0], normal[1], normal[2]
    P = torch.zeros((3, 6), dtype=torch.float32)

    # t_x = σ_xx * n_x + σ_xy * n_y + σ_xz * n_z
    P[0, 0] = n0
    P[0, 5] = n1
    P[0, 4] = n2

    # t_y = σ_yx * n_x + σ_yy * n_y + σ_yz * n_z
    P[1, 5] = n0
    P[1, 1] = n1
    P[1, 3] = n2

    # t_z = σ_zx * n_x + σ_zy * n_y + σ_zz * n_z
    P[2, 4] = n0
    P[2, 3] = n1
    P[2, 2] = n2

    return P


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
    n_fft = 2 ** int(math.ceil(math.log2(N)))  # 2的幂次，保障FFT高效

    # 执行FFT，使用forward确保严格能量一致性
    A = torch.fft.rfft(a, n=n_fft, dim=0)  # [fft_len, d]
    B = torch.fft.rfft(b, n=n_fft, dim=0)
    C = A * B  # 逐通道频域乘积

    result = torch.fft.irfft(C, n=n_fft, dim=0)[:N]  # [N, d]

    if sum_dim:
        result = torch.sum(result, dim=1)  # 所有通道累加，返回 [N]

    return result


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

# 定义三维 8 个 Gauss 点位置
sqrt3_inv = 1.0 / math.sqrt(3)
gauss_3d = torch.tensor([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                         [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]) * sqrt3_inv


#######################################
#######################################
#######################################


class Element:
    def __init__(self, gamma_face_index, device="cpu"):
        self.gamma_face_index = gamma_face_index
        self.device = device
        self.face_disp2traction = None
        self.face_node2gauss = None

    def save_to(self, filename):
        """
        保存 Element 关键矩阵到文件，便于重用、跨平台部署
        :param filename: 保存路径（.pt 文件）
        """
        state = {
            'gamma_face_index': self.gamma_face_index,
            'face_disp2traction': self.face_disp2traction.cpu(),  # 存盘统一转CPU，兼容性好
            'face_node2gauss': self.face_node2gauss.cpu()
        }
        torch.save(state, filename)

    @classmethod
    def load_from(cls, filename, device="cpu"):
        """
        从文件加载 Element 矩阵，返回新对象，支持切换计算设备
        :param filename: .pt 文件路径
        :param device: 目标计算设备
        :return: Element 实例
        """
        state = torch.load(filename, map_location=device, weights_only=False)

        # 跳过 __init__ 直接创建对象
        obj = cls.__new__(cls)

        # 手动填充属性
        obj.gamma_face_index = state['gamma_face_index']
        obj.face_disp2traction = state['face_disp2traction'].to(device)
        obj.face_node2gauss = state['face_node2gauss'].to(device)
        obj.device = device

        return obj

    def _compute_traction(self, u):
        """
        计算面上牵引力
        :param u: [time, 8, 3]  时间、节点、方向
        :return: traction [time, 4, 3]  时间、节点、方向
        """
        u = u.to(self.device)  # u: [time, 8, 3]
        u = u.reshape(u.shape[0], -1)  # u: [time, 24]

        # face_disp2traction: [4, 3, 24]
        # u: [time, 24]
        # traction: [time, 4, 3]
        traction = torch.einsum('gcn, tn -> tgc', self.face_disp2traction, u)
        return traction

    def _compute_disp_gauss(self, u):
        """
        计算面上高斯点的位移
        :param u: [time, 8, 3]  时间、节点、方向
        :return: disp [time, 4, 3]  时间、Gauss点、方向
        """
        u = u.to(self.device)  # u: [time, 8, 3]

        # face_node2gauss: [4, 8]
        # u: [time, 8, 3]
        # disp: [time, 4, 3]
        disp = torch.einsum('gn, tnc -> tgc', self.face_node2gauss, u)

        return disp

    def compute_convolve(self, u_near, u_recip):
        """
        计算面上双场卷积积分，返回理论表达式中的最终结果：

            s(t) = ∑_g ∑_i ( t_i^near ∗ u_i^recip - u_i^near ∗ t_i^recip )

        其中：
        - ∗ 表示时间卷积
        - g 表示面上高斯点编号（共4个）
        - i 表示方向（3个方向）
        - u_near, t_near：近场位移、牵引力
        - u_recip, t_recip：互反场位移、牵引力

        参数：
            u_near : [time, 8, 3]
                近场单元8节点、3方向的多时刻位移
            u_recip : [time, 8, 3]
                互反场单元8节点、3方向的多时刻位移

        返回：
            convolved : [total_time]
                面上积分、方向累加后的最终卷积结果
        """

        # Step 1：计算面上牵引力，拉平为 [time, 12]，即 4 Gauss点 × 3方向
        t_near = self._compute_traction(u_near).reshape(u_near.shape[0], -1)
        t_recip = self._compute_traction(u_recip).reshape(u_recip.shape[0], -1)

        # Step 2：计算面上Gauss点位移，拉平为 [time, 12]
        ug_near = self._compute_disp_gauss(u_near).reshape(u_near.shape[0], -1)
        ug_recip = self._compute_disp_gauss(u_recip).reshape(u_recip.shape[0], -1)

        # Step 3：时间卷积
        # t_near ∗ u_recip：近场牵引力与互反场位移卷积
        first = fft_convolve_multidim(t_near, ug_recip, sum_dim=True)

        # u_near ∗ t_recip：近场位移与互反场牵引力卷积
        second = fft_convolve_multidim(t_recip, ug_near, sum_dim=True)

        # Step 4：计算最终贡献，两个方向相减
        convolved = first - second  # [total_time]

        return convolved


class SolidElement(Element):
    def __init__(self, node_pos, lambda_n, mu_n, gamma_face_index, device="cpu"):
        """
        初始化 SolidElement，仅输入节点物理参数，自动插值至Gauss点
        :param node_pos: [8, 3] 节点坐标
        :param lambda_n: [8] 节点Lamé参数 λ
        :param mu_n: [8] 节点剪切模量 μ
        :param gamma_face_index: int，面编号 1~6
        """
        super().__init__(gamma_face_index, device)
        self._precompute(node_pos, lambda_n, mu_n)

    def _precompute(self, node_pos, lambda_n, mu_n):
        """
        预计算面上各矩阵：形函数、位移→应力、应力→牵引力
        物理参数通过节点值自动插值至Gauss点
        """
        # 面定义
        face_node_idx = face_node_dict[self.gamma_face_index]
        dim, pos = face_dim_direction_dict[self.gamma_face_index]

        # 节点坐标
        if not isinstance(node_pos, torch.Tensor):
            node_pos = torch.tensor(node_pos, dtype=torch.float32)
        node_pos_face = node_pos[face_node_idx]

        # 面上Gauss点局部坐标
        gauss_pos_face = gauss_3d[face_node_idx]
        gauss_pos_face[:, dim] = pos

        # 预计算矩阵
        face_disp2traction, face_node2gauss = [], []
        for gp in range(4):
            xi, eta, zeta = gauss_pos_face[gp]

            # 形函数：节点值插值到Gauss点
            N = _shape_function(xi, eta, zeta)
            face_node2gauss.append(N)

            lam = torch.sum(N * torch.tensor(lambda_n, dtype=torch.float32))
            mu = torch.sum(N * torch.tensor(mu_n, dtype=torch.float32))

            # 位移->应变
            dN_dxi = _shape_function_derivatives(xi, eta, zeta)
            J = dN_dxi.T @ node_pos
            invJ = torch.inverse(J)
            dN_dx = dN_dxi @ invJ.T
            B = _build_B_matrix(dN_dx)  # [6, 24]

            # 位移->应力
            D = torch.zeros((6, 6))  # [6, 6]
            D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
            D[3, 3] = D[4, 4] = D[5, 5] = mu
            D[0, 1] = D[0, 2] = D[1, 2] = D[1, 0] = D[2, 0] = D[2, 1] = lam

            # 应力->牵引力
            F = compute_normal(node_pos_face, gp, voigt_solid=True)  # [3, 6]

            # 合并
            face_disp2traction.append(F @ D @ B)  # [3, 6, 24]

        # 存储
        self.face_disp2traction = torch.stack(face_disp2traction, dim=0).to(self.device)
        self.face_node2gauss = torch.stack(face_node2gauss, dim=0).to(self.device)


class FluidElement(Element):
    def __init__(self, node_pos, rho_n, gamma_face_index, device="cpu"):
        """
        初始化 FluidElement，节点输入，自动插值密度
        :param node_pos: (8, 3) 节点坐标
        :param rho_n: (8,) 节点密度 ρ
        :param gamma_face_index: 1~6 指定面编号
        """
        super().__init__(gamma_face_index, device)
        self._precompute(node_pos, rho_n)

    def _precompute(self, node_pos, rho_n):
        """
        预计算面上所有矩阵，密度通过形函数自动插值
        """
        # 面节点
        face_node_idx = face_node_dict[self.gamma_face_index]
        dim, pos = face_dim_direction_dict[self.gamma_face_index]

        # 坐标准备
        if not isinstance(node_pos, torch.Tensor):
            node_pos = torch.tensor(node_pos, dtype=torch.float32)
        node_pos_face = node_pos[face_node_idx]

        # 面上 Gauss 点局部坐标
        gauss_pos_face = gauss_3d[face_node_idx]
        gauss_pos_face[:, dim] = pos

        face_disp2traction, face_node2gauss = [], []

        for gp in range(4):
            xi, eta, zeta = gauss_pos_face[gp]

            # 形函数：节点值插值到 Gauss 点
            N = _shape_function(xi, eta, zeta)
            face_node2gauss.append(N)

            rho = torch.sum(N * torch.tensor(rho_n, dtype=torch.float32))  # 密度插值

            # 势函数 -> 位移梯度
            dN_dxi = _shape_function_derivatives(xi, eta, zeta)
            J = dN_dxi.T @ node_pos
            invJ = torch.inverse(J)
            dN_dx = dN_dxi @ invJ.T
            B = dN_dx.T  # [3, 8]

            # 位移梯度 -> 牵引力
            F = compute_normal(node_pos_face, gp, voigt_solid=False)  # [1, 3]

            # 合并
            face_disp2traction.append(F @ (B / rho))  # [1, 8]

        self.face_disp2traction = torch.stack(face_disp2traction, dim=0).to(self.device)
        self.face_node2gauss = torch.stack(face_node2gauss, dim=0).to(self.device)


def main():
    """
    生成单位立方体，顶面施加均匀竖直位移，打印面上Gauss点的位移与牵引力
    """
    # 单元节点坐标，立方体 [-0.5, 0.5]
    node_pos = torch.tensor([
        [-0.5, -0.5, -0.5],
        [0.5, -0.5, -0.5],
        [0.5, 0.5, -0.5],
        [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5],
        [0.5, -0.5, 0.5],
        [0.5, 0.5, 0.5],
        [-0.5, 0.5, 0.5],
    ], dtype=torch.float32)

    # 物性参数，均匀各向同性
    lambda_n = torch.ones(8) * 2.0
    mu_n = torch.ones(8) * 1.0

    # 选择顶面 z = +0.5
    gamma_face = 6

    # 初始化Solid单元
    elem = SolidElement(node_pos, lambda_n, mu_n, gamma_face_index=gamma_face)

    # 生成位移场：顶面节点（4,7,6,5）z方向均匀+0.01，其他为0
    u = torch.zeros(1, 8, 3)
    top_nodes = [4, 5, 6, 7]
    u[0, top_nodes, 0] = 1.0

    # 计算Gauss点位移与牵引力
    ug = elem._compute_disp_gauss(u)[0]  # [4,3]
    tg = elem._compute_traction(u)[0]  # [4,3]

    print("面上4个Gauss点位移 [x, y, z]：")
    print(ug)

    print("\n面上4个Gauss点牵引力 [x, y, z]：")
    print(tg)


if __name__ == "__main__":
    main()
