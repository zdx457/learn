import numpy as np
from tool import enu_to_ecef,ecef_to_lat_lon_alt
from constant_grounp import user_receive_antenna_gain

class User:
    def __init__(self, bs_x, bs_y, bs_r, reference_lat_lon_alt):
        self.transmitted_power = -7  # 23 dBm-30=dBw
        self.antenna_gain = user_receive_antenna_gain  # dBi 接收天线增益
        self.h = 1.5  # 人的身高 1.5m
        self.coordinate_enu = self.generate_user(bs_x, bs_y, bs_r)  # ENU
        self.coordinate_ecef = enu_to_ecef(self.coordinate_enu[0], self.coordinate_enu[1], self.h, reference_lat_lon_alt)  # ECEF
        self.coordinate_lat_lon_alt = ecef_to_lat_lon_alt(self.coordinate_ecef[0], self.coordinate_ecef[1], self.coordinate_ecef[2])  # lat lon alt
        # self.interferece_range = np.sqrt(3)*2 * bs_r # 受扰范围
        self.interferece_from_bs = [] # 受扰的基站对象
        self.beam_index=None

    # 生成用户坐标
    def generate_user(self, bs_x, bs_y, bs_r = 0):
        r_k_m = np.sqrt(3)/2 * bs_r * np.random.rand() #保证用户在基站里面
        a_k_m = 2 * np.pi * np.random.rand()
        x1 = bs_x + r_k_m * np.cos(a_k_m)
        y2 = bs_y + r_k_m * np.sin(a_k_m)
        return [x1, y2, self.h]  # 单个用户的坐标