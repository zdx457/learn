import numpy as np
from tool import enu_to_ecef,ecef_to_lat_lon_alt,ecef_to_enu,lat_lon_alt_to_ecef

class Hexagon:
    def __init__(self, x_initial, y_initial, r, reference_lat_lon_alt):
        self.h = 0.0  # 六边形在平面上，高度为 0 m
        self.coordinate_enu = self.generate_Hexagon(x_initial, y_initial, r) # ENU
        self.coordinate_ecef = [enu_to_ecef(E, N, U, reference_lat_lon_alt) for E, N, U in self.coordinate_enu]  # ECEF
        self.coordinate_lat_lon_alt = [ecef_to_lat_lon_alt(X ,Y ,Z) for X ,Y ,Z in self.coordinate_ecef] # lat lon alt

        # self.coordinate_ecef_to_enu = [ecef_to_enu(x ,y ,z, reference_lat_lon_alt) for x ,y ,z in self.coordinate_ecef]
        # self.coordinate_lat_lon_alt_to_ecef = [lat_lon_alt_to_ecef(lat ,lon ,alt) for lat ,lon ,alt in self.coordinate_lat_lon_alt]
        # self.coordinate_ecef_to_lat_lon_alt = [ecef_to_lat_lon_alt(x ,y ,z) for x ,y ,z in self.coordinate_ecef]


    def generate_Hexagon(self, x_center, y_center, r):
        n = 6
        a = np.linspace(0, 2 * np.pi, n + 1)  # 基准蜂窝的六角定位（七个点,便于绘图）
        x0 = x_center + r * np.cos(a)
        y0 = y_center + r * np.sin(a)  # 转化为坐标
        hexagon = [[x, y, self.h] for x, y in zip(x0, y0)]
        return hexagon  # list
