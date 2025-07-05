import os
import torch
from math import sqrt


class BaseElement:
    def __init__(self, device='cpu'):
        self.device = device
        self.K_total = None
        self.M_total = None

    @staticmethod
    def _shape_function_derivatives(xi, eta, zeta, device='cpu'):
        """
        计算8节点Brick单元形函数对自然坐标的偏导 (8, 3)
        """
        dN = torch.zeros((8, 3), device=device)
        coeff = 0.125
        signs = torch.tensor([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                              [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], device=device)
        for i in range(8):
            sx, sy, sz = signs[i]
            dN[i, 0] = coeff * sx * (1 + sy * eta) * (1 + sz * zeta)
            dN[i, 1] = coeff * (1 + sx * xi) * sy * (1 + sz * zeta)
            dN[i, 2] = coeff * (1 + sx * xi) * (1 + sy * eta) * sz
        return dN

    @staticmethod
    def _shape_function_value(i, xi, eta, zeta):
        """
        计算第i个节点的标准形函数值 (标量)
        """
        coeff = 0.125
        signs = [-1, 1, 1, -1, -1, 1, 1, -1]
        sx = signs[i % 4]
        sy = signs[(i // 2) % 2]
        sz = signs[i // 4]
        return coeff * (1 + sx * xi) * (1 + sy * eta) * (1 + sz * zeta)

    def save_to(self, name, save_dir):
        """
        保存刚度与质量矩阵
        """
        os.makedirs(save_dir, exist_ok=True)
        if not name.endswith('.pt'):
            name += '.pt'
        torch.save({'K': self.K_total.cpu(), 'M': self.M_total.cpu()}, f'{save_dir}/{name}')


class SolidElement(BaseElement):
    def __init__(self, node_position, lambda_g, mu_g, rho_g, device='cpu'):
        super().__init__(device)
        self.K_total = torch.zeros(24, 24, device=self.device)
        self.M_total = torch.zeros(24, 24, device=self.device)
        self._precompute(node_position, lambda_g, mu_g, rho_g)

    def _precompute(self, node_position, lambda_g, mu_g, rho_g):
        node_pos = torch.tensor(node_position, dtype=torch.float32, device=self.device)
        lambda_g = torch.tensor(lambda_g, dtype=torch.float32, device=self.device)
        mu_g = torch.tensor(mu_g, dtype=torch.float32, device=self.device)
        rho_g = torch.tensor(rho_g, dtype=torch.float32, device=self.device)

        sqrt3_inv = 1.0 / sqrt(3)
        gauss_pos = torch.tensor([-sqrt3_inv, sqrt3_inv], device=self.device)
        gauss_points = torch.stack(torch.meshgrid(gauss_pos, gauss_pos, gauss_pos, indexing='ij'),
                                   dim=-1).reshape(-1, 3)

        for gp in range(8):
            xi, eta, zeta = gauss_points[gp]
            dN_dxi = self._shape_function_derivatives(xi, eta, zeta, device=self.device)
            J = dN_dxi.T @ node_pos
            invJ = torch.inverse(J)
            detJ = torch.det(J)
            dN_dx = dN_dxi @ invJ.T
            B = self._build_B_matrix(dN_dx, device=self.device)

            lam = lambda_g[gp]
            mu = mu_g[gp]
            rho = rho_g[gp]

            # 刚度矩阵
            D = torch.zeros((6, 6), device=self.device)
            D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
            D[3, 3] = D[4, 4] = D[5, 5] = mu
            D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = lam

            K_local = B.T @ D @ B * detJ
            self.K_total += K_local

            # 质量矩阵
            N_matrix = torch.zeros(3, 24, device=self.device)
            for i in range(8):
                idx = i * 3
                N_val = self._shape_function_value(i, xi, eta, zeta)
                N_matrix[0, idx] = N_val
                N_matrix[1, idx + 1] = N_val
                N_matrix[2, idx + 2] = N_val

            M_local = N_matrix.T @ N_matrix * detJ * rho
            self.M_total += M_local

    @staticmethod
    def _build_B_matrix(dN_dx, device='cpu'):
        B = torch.zeros((6, 24), device=device)
        for i in range(8):
            idx = i * 3
            dNxi, dNyi, dNzi = dN_dx[i]
            B[0, idx] = dNxi
            B[1, idx + 1] = dNyi
            B[2, idx + 2] = dNzi
            B[3, idx] = dNyi
            B[3, idx + 1] = dNxi
            B[4, idx + 1] = dNzi
            B[4, idx + 2] = dNyi
            B[5, idx] = dNzi
            B[5, idx + 2] = dNxi
        return B

    def compute_internal_force(self, displacement):
        """
        :param displacement: (t, 8, 3)
        :return: (t, 8, 3) 内力
        """
        t = displacement.shape[0]
        u = displacement.reshape(t, -1)
        f_int = (self.K_total @ u.T.to(self.device)).T
        return f_int.reshape(t, 8, 3)

    def compute_mass_force(self, displacement, dt):
        """
        :param displacement: (t, 8, 3) 多时刻位移
        :param dt: 时间步长，标量
        :return: (t-2, 8, 3) 质量内力，差分法近似加速度
        """
        assert displacement.shape[0] >= 3, "至少提供3个时刻的位移用于差分"

        # 差分二阶时间导数
        u_prev = displacement[:-2]  # (t-2, 8, 3)
        u_curr = displacement[1:-1]
        u_next = displacement[2:]

        acc = (u_next - 2 * u_curr + u_prev) / (dt ** 2)  # (t-2, 8, 3)

        t_eff = acc.shape[0]
        acc_flat = acc.reshape(t_eff, -1)

        f_mass = (self.M_total @ acc_flat.T.to(self.device)).T
        return f_mass.reshape(t_eff, 8, 3)

    def compute_total_force(self, displacement, dt):
        return self.compute_internal_force(displacement)[1:-1] + self.compute_mass_force(displacement, dt)

    @staticmethod
    def load_from(path, device='cpu'):
        data = torch.load(path, map_location=device, weights_only=False)
        obj = SolidElement.__new__(SolidElement)
        obj.device = device
        obj.K_total = data['K'].to(device)
        obj.M_total = data['M'].to(device)
        return obj


class FluidElement(BaseElement):
    def __init__(self, node_position, rho_g, kappa_g, device='cpu'):
        super().__init__(device)
        self.K_total = torch.zeros(8, 8, device=self.device)
        self.M_total = torch.zeros(8, 8, device=self.device)
        self._precompute(node_position, rho_g, kappa_g)

    def _precompute(self, node_position, rho_g, kappa_g):
        node_pos = torch.tensor(node_position, dtype=torch.float32, device=self.device)
        rho_g = torch.tensor(rho_g, dtype=torch.float32, device=self.device)
        kappa_g = torch.tensor(kappa_g, dtype=torch.float32, device=self.device)

        sqrt3_inv = 1.0 / sqrt(3)
        gauss_pos = torch.tensor([-sqrt3_inv, sqrt3_inv], device=self.device)
        gauss_points = torch.stack(torch.meshgrid(gauss_pos, gauss_pos, gauss_pos, indexing='ij'),
                                   dim=-1).reshape(-1, 3)

        for gp in range(8):
            xi, eta, zeta = gauss_points[gp]
            dN_dxi = self._shape_function_derivatives(xi, eta, zeta, device=self.device)
            J = dN_dxi.T @ node_pos
            invJ = torch.inverse(J)
            detJ = torch.det(J)
            dN_dx = dN_dxi @ invJ.T  # (8, 3)

            rho = rho_g[gp]
            kappa = kappa_g[gp]

            # 刚度矩阵：∫ (1/rho) ∇N^T ∇N dV
            gradN = dN_dx  # (8, 3)
            K_local = gradN @ gradN.T * detJ / rho
            self.K_total += K_local

            # 质量矩阵：∫ (1/kappa) N^T N dV
            N_vec = torch.zeros(8, device=self.device)
            for i in range(8):
                N_vec[i] = self._shape_function_value(i, xi, eta, zeta)

            M_local = torch.outer(N_vec, N_vec) * detJ / kappa
            self.M_total += M_local

    def compute_internal_force(self, potential):
        """
        :param potential: (t, 8, 1) 位势场
        :return: (t, 8, 1) 内力
        """
        t = potential.shape[0]
        chi = potential.reshape(t, -1)
        f_int = (self.K_total @ chi.T.to(self.device)).T
        return f_int.reshape(t, 8, 1)

    def compute_mass_force(self, potential, dt):
        """
        :param potential: (t, 8, 1) 位势场
        :param dt: 时间步长
        :return: (t-2, 8, 1) 质量内力
        """
        assert potential.shape[0] >= 3, "至少提供3个时刻位势用于差分"

        u_prev = potential[:-2].reshape(-1, 8)
        u_curr = potential[1:-1].reshape(-1, 8)
        u_next = potential[2:].reshape(-1, 8)

        acc = (u_next - 2 * u_curr + u_prev) / (dt ** 2)
        f_mass = (self.M_total @ acc.T.to(self.device)).T
        return f_mass.reshape(-1, 8, 1)

    def compute_total_force(self, potential, dt):
        return self.compute_internal_force(potential)[1:-1] + self.compute_mass_force(potential, dt)

    @staticmethod
    def load_from(path, device='cpu'):
        data = torch.load(path, map_location=device, weights_only=False)
        obj = FluidElement.__new__(FluidElement)
        obj.device = device
        obj.K_total = data['K'].to(device)
        obj.M_total = data['M'].to(device)
        return obj
