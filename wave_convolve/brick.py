import numpy as np
import torch

# 节点编号示意：
#   7 -------- 6
#  /|         /|
# 4 -------- 5 |
# | |        | |
# | 3 -------| 2
# |/         |/
# 0 -------- 1

# 定义三维 8 个 Gauss 点位置
sqrt3_inv = 1.0 / np.sqrt(3)
gauss_1d = torch.tensor([-sqrt3_inv, sqrt3_inv])
gauss_3d = torch.stack(torch.meshgrid(
    gauss_1d, gauss_1d, gauss_1d, indexing='ij'), dim=-1).reshape(-1, 3)


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
    3 - ε_yz
    4 - ε_xz
    5 - ε_xy
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
        self.face_stress_mat = None  # [4, 6, 24]
        self.face_project_matrix = None  # [4, 3, 6]
        self.face_shape_func = None  # [4, 8]
        self._precompute(node_pos, lambda_g, mu_g)

    def _precompute(self, node_pos, lambda_g, mu_g):
        """
        预计算面上所有矩阵
        """
        node_pos = torch.tensor(node_pos, dtype=torch.float32)
        face_node_dict = {
            1: [0, 4, 7, 3],  # x = -1
            2: [1, 2, 6, 5],  # x = +1
            3: [0, 1, 5, 4],  # y = -1
            4: [2, 3, 7, 6],  # y = +1
            5: [0, 3, 2, 1],  # z = -1
            6: [4, 5, 6, 7]  # z = +1
        }
        face_dim_direction_dict = {
            1: [0, -1],
            2: [0, 1],
            3: [1, -1],
            4: [1, 1],
            5: [2, -1],
            6: [2, 1],
        }
        face_node_idx = face_node_dict[self.gamma_face_index]
        node_pos_face = node_pos[face_node_idx]
        gauss_pos_face = gauss_3d[face_node_idx]
        dim, pos = face_dim_direction_dict[self.gamma_face_index]
        gauss_pos_face[:, dim] = pos

        face_stress_mat, face_project_matrix, face_shape_func = [], [], []

        for gp in range(4):
            xi, eta, zeta = gauss_pos_face[gp]
            N = _shape_function(xi, eta, zeta)
            face_shape_func.append(N)

            dN_dxi = _shape_function_derivatives(xi, eta, zeta)
            J = dN_dxi.T @ node_pos
            invJ = torch.inverse(J)
            dN_dx = dN_dxi @ invJ.T
            B = _build_B_matrix(dN_dx)

            lam, mu = lambda_g[gp], mu_g[gp]
            D = torch.zeros((6, 6))
            D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
            D[3, 3] = D[4, 4] = D[5, 5] = mu
            D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = lam

            face_stress_mat.append(D @ B)
            face_project_matrix.append(self._compute_face_projection_matrix(node_pos_face, gp))

        self.face_stress_mat = torch.stack(face_stress_mat, dim=0).to(self.device)
        self.face_project_matrix = torch.stack(face_project_matrix, dim=0).to(self.device)
        self.face_shape_func = torch.stack(face_shape_func, dim=0).to(self.device)

    @staticmethod
    def _compute_face_projection_matrix(face_node_pos, gp_index):
        """
        计算面上单位法向与投影矩阵
        """
        point_self = face_node_pos[gp_index]
        point_next = face_node_pos[(gp_index + 1) % 4]
        point_prev = face_node_pos[gp_index - 1]
        vec1 = point_self - point_prev
        vec2 = point_next - point_self
        normal = torch.cross(vec1, vec2)
        normal = normal / torch.norm(normal)

        n0, n1, n2 = normal[0], normal[1], normal[2]
        P = torch.zeros((3, 6), dtype=torch.float32)
        P[0, 0] = n0
        P[0, 5] = n1
        P[0, 4] = n2
        P[1, 5] = n0
        P[1, 1] = n1
        P[1, 3] = n2
        P[2, 4] = n0
        P[2, 3] = n1
        P[2, 2] = n2
        return P

    def _compute_traction(self, u):
        """
        计算面上牵引力
        :param u: [time, 8, 3]
        :return: [4, 3, time]
        """
        u = u.to(self.device)
        n_time = u.shape[0]
        u = u.reshape(n_time, 24).T
        stress = self.face_stress_mat @ u  # [4, 6, time]
        traction = torch.einsum('gij, gjt -> git', self.face_project_matrix, stress)
        return traction

    def _compute_displacement_gauss(self, u):
        """
        计算面上高斯点位移
        :param u: [time, 8, 3]
        :return: [time, 4, 3]  时间、Gauss点、方向
        """
        u = u.to(self.device)
        disp = torch.einsum('gn, tnc -> tgc', self.face_shape_func, u)
        return disp
