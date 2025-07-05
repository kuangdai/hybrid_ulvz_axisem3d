import numpy as np
import os

from brick import SolidElement
from prem import PREM


def main():
    # 路径设置
    grid_path = '../wave_extrapolation/outputs'
    prem_path = './prem_1s.csv'
    save_dir = "./element_cache"
    os.makedirs(save_dir, exist_ok=True)

    R_earth = 6371000.0  # 地球半径，单位 m

    # 读取网格
    phi = np.loadtxt(f'{grid_path}/grid_azim_ml_ready.txt')  # (Np,)
    theta = np.loadtxt(f'{grid_path}/grid_dist_ml_ready.txt')  # (Nt,)
    depth = np.loadtxt(f'{grid_path}/grid_depth_ml_ready_solid.txt') * 1000  # m

    Np = phi.shape[0]
    Nt = theta.shape[0]
    Nd = depth.shape[0]

    print(f'Grid: depth={Nd}, theta={Nt}, phi={Np}')

    prem = PREM(prem_path)

    ### 1. 最浅 depth 层单元 ###
    for ip in range(Np):
        for it in range(Nt - 1):
            ip_next = (ip + 1) % Np

            indices = [
                (0, it, ip),
                (0, it, ip_next),
                (0, it + 1, ip),
                (0, it + 1, ip_next),
                (1, it, ip),
                (1, it, ip_next),
                (1, it + 1, ip),
                (1, it + 1, ip_next)
            ]

            nodes, lambda_g, mu_g, rho_g = [], [], [], []

            for d, t, p in indices:
                r_m = R_earth - depth[d]
                th = theta[t]
                ph = phi[p]

                x = r_m * np.sin(th) * np.cos(ph)
                y = r_m * np.sin(th) * np.sin(ph)
                z = r_m * np.cos(th)
                nodes.append([x, y, z])

                r_m_other = R_earth - depth[indices[0][0] + indices[4][0] - d]
                delta_r = r_m - r_m_other
                factor = np.sqrt(3) / 2
                r_gauss = r_m_other + factor * delta_r
                lam, mu, rho = prem.query_solid(r_gauss)
                lambda_g.append(lam)
                mu_g.append(mu)
                rho_g.append(rho)

            nodes = np.array(nodes)
            element = SolidElement(nodes, lambda_g, mu_g, rho_g)
            name = f'solid_top_p{ip}_t{it}'
            element.save_to(name, save_dir)

            # 存储connectivity
            connect_path = f"{save_dir}/{name}.connect.txt"
            with open(connect_path, 'w') as f:
                for d, t, p in indices:
                    f.write(f"{d} {t} {p}\n")

    print(f'✅ 最浅depth层单元共 {Np} × {Nt - 1} 个完成')

    ### 2. 最大 theta 层单元 ###
    for ip in range(Np):
        for id in range(1, Nd - 1):  # 从id=1开始，避免与最浅层重叠
            ip_next = (ip + 1) % Np

            indices = [
                (id, -1, ip),
                (id, -1, ip_next),
                (id, -2, ip),
                (id, -2, ip_next),
                (id + 1, -1, ip),
                (id + 1, -1, ip_next),
                (id + 1, -2, ip),
                (id + 1, -2, ip_next)
            ]

            nodes, lambda_g, mu_g, rho_g = [], [], [], []

            for d, t, p in indices:
                r_m = R_earth - depth[d]
                th = theta[t]
                ph = phi[p]

                x = r_m * np.sin(th) * np.cos(ph)
                y = r_m * np.sin(th) * np.sin(ph)
                z = r_m * np.cos(th)
                nodes.append([x, y, z])

                r_m_other = R_earth - depth[indices[0][0] + indices[4][0] - d]
                delta_r = r_m - r_m_other
                factor = np.sqrt(3) / 2
                r_gauss = r_m_other + factor * delta_r
                lam, mu, rho = prem.query_solid(r_gauss)
                lambda_g.append(lam)
                mu_g.append(mu)
                rho_g.append(rho)

            nodes = np.array(nodes)
            element = SolidElement(nodes, lambda_g, mu_g, rho_g)
            name = f'solid_side_p{ip}_d{id}'
            element.save_to(name, save_dir)

            # 存储connectivity
            connect_path = f"{save_dir}/{name}.connect.txt"
            with open(connect_path, 'w') as f:
                for d, t, p in indices:
                    f.write(f"{d} {t} {p}\n")

    print(f'✅ 最大theta层单元共 {Np} × {Nd - 2} 个完成')


if __name__ == '__main__':
    main()
