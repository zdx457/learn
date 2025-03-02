# constellation.py

import numpy as np
from satellite import Satellite
from constant import R_e, mu, time_step


class Constellation:
    def __init__(self, h, P, S, i_deg = 46, F = 0):
        """
        初始化星座对象。

        参数：
        - h: 卫星高度（km）
        - P: 轨道平面数
        - S: 每个平面的卫星数
        - i_deg: 轨道倾角（度）
        - F: 相位因子
        """
        self.h = h
        self.P = P
        self.S = S
        self.i_deg = i_deg
        self.F = F
        self.satellites = []
        self._initialize_constellation()

    def _initialize_constellation(self):
        """
        初始化星座，生成所有卫星对象。
        """
        default_wavelength = 0.15
        default_pt = 44.3
        a = R_e + self.h
        n = np.sqrt(mu / a ** 3)
        n_deg = np.degrees(n)

        # 遍历每个轨道平面和每个卫星
        for j in range(1, self.P + 1):
            for k in range(1, self.S + 1):
                Omega_jk = (360.0 / self.P * (j - 1)) % 360.0
                f_jk = (360.0 * (self.F / (self.P * self.S) * (j - 1) + 1.0 / self.S * (k - 1))) % 360.0

                orbital_elements = {
                    'a': a,
                    'e': 0.0,
                    'i': np.radians(self.i_deg),
                    'Omega': np.radians(Omega_jk),
                    'omega': 0.0,
                    'f': np.radians(f_jk),
                    'n': n
                }

                satellite_id = (j - 1) * self.S + k - 1
                satellite = Satellite(satellite_id, orbital_elements,default_wavelength, default_pt)
                self.satellites.append(satellite)
    def update_velocity(self):
        for satellite in self.satellites:
            satellite.calculate_velocity(time_step)

    def update_positions(self, delta_t):
        """
        更新所有卫星的位置。

        参数：
        - delta_t: 时间增量（秒）
        """
        for satellite in self.satellites:
            satellite.update_position(delta_t)
