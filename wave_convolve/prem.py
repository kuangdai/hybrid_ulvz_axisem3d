import numpy as np

CMB_RADIUS = 3480000.0  # m
EPSILON = 0.01  # m


class PREM:
    def __init__(self, prem_file):
        """
        读取 PREM 数据，全部转为 SI 单位：
        - 半径 r: m
        - Vp, Vs: m/s
        - rho: kg/m³
        """
        data = np.loadtxt(prem_file, delimiter=',')

        radius_km = data[:, 0]  # 单位 km
        self.radius = radius_km * 1000  # 单位 m

        self.vp = data[:, 3] * 1000  # m/s
        self.vs = data[:, 4] * 1000  # m/s
        self.rho = data[:, 2] * 1000  # kg/m³

    CMB_RADIUS = 3480000.0  # m
    EPSILON = 0.01  # m

    def query_solid(self, r_m):
        """
        输入：半径 r，单位 m
        输出：lambda, mu, rho，全为 SI 单位
        """
        if abs(r_m - CMB_RADIUS) < EPSILON:
            r_m += EPSILON
        vp = np.interp(r_m, self.radius, self.vp)
        vs = np.interp(r_m, self.radius, self.vs)
        rho = np.interp(r_m, self.radius, self.rho)

        mu = rho * vs ** 2
        lam = rho * vp ** 2 - 2 * mu
        return lam, mu, rho

    def query_fluid(self, r_m):
        """
        输入：半径 r，单位 m
        输出：kappa, rho，单位 SI
        """
        if abs(r_m - CMB_RADIUS) < EPSILON:
            r_m -= EPSILON
        vp = np.interp(r_m, self.radius, self.vp)
        rho = np.interp(r_m, self.radius, self.rho)
        kappa = vp ** 2 / rho
        return kappa, rho
