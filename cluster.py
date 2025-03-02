from constant_grounp import *
from create_BS import BS
from create_Hexagon import Hexagon
from create_User import User
from tool import *
import numpy as np
import matplotlib.pyplot as plt





class Cluster:
    def __init__(self, region):
        self.region_lat_lon_alt = region  # 纬度，经度，海拔
        self.region_center_lat_lon_alt = None  # 区域中心经纬度海拔[lat, lon, alt]
        self.region_center_ecef = None  # 区域中心 ECEF 坐标[X, Y, Z]
        self.region_center_enu = None  # 区域中心 ENU 坐标
        self.region_ecef = None  # 区域 ECEF 坐标 list [X,Y,Z]
        self.region_enu = None  # list [E,N,U]
        self.Base_Station = []  # 基站对象
        self.hexagon_r = bs_radius
        self.user_num = user_num_in_bs

    def generate_all_BS(self):
        """
        生成此 enu 区域所有基站
        :return: None
        """
        # 估计区域中心的经纬度海拔
        self.region_center_lat_lon_alt = calculate_center(self.region_lat_lon_alt)  # 区域中心经纬度海拔[lat, lon, alt]
        self.region_center_ecef = lat_lon_alt_to_ecef(*self.region_center_lat_lon_alt)  # 区域中心 ECEF 坐标[X, Y, Z]
        self.region_center_enu = ecef_to_enu(*self.region_center_ecef, reference_lat_lon_alt=self.region_center_lat_lon_alt)  # 区域中心 ENU 坐标
        self.region_ecef = [lat_lon_alt_to_ecef(lat, lon, h) for lat, lon, h in self.region_lat_lon_alt]  # 区域 ECEF 坐标 list [X,Y,Z]
        self.region_enu = [ecef_to_enu(X, Y, Z, reference_lat_lon_alt=self.region_center_lat_lon_alt) for X, Y, Z in self.region_ecef]  # list [E,N,U]

        region_enu_array = np.array(self.region_enu)
        max_x = max(region_enu_array[:, 0])
        min_x = min(region_enu_array[:, 0])
        max_y = max(region_enu_array[:, 1])
        min_y = min(region_enu_array[:, 1])

        x0 = (max_x + min_x) / 2
        y0 = (max_y + min_y) / 2

        dx = 3 / 2 * self.hexagon_r  # 横向两个 cell 中心的距离
        dy = np.sqrt(3) * self.hexagon_r  # 纵向两个 cell 中心的距离
        m = int((x0 - min_x) / dx)  # 区群列数    以区域经纬度坐标为中心，原点
        c = 0  # 基站记数
        # 下面计算的是 enu 坐标
        for i in range(-m, m + 1):
            if abs(i) % 2 == 0:
                n = int((max_y - y0) / dy)  # 行数
                for j in range(-n, n + 1):
                    bs_xtest = x0 + i * dx
                    bs_ytest = y0 + j * dy
                    point = [bs_xtest, bs_ytest]
                    if is_point_in_polygon(point, region_enu_array):
                        hexagon = Hexagon(*point, r=self.hexagon_r, reference_lat_lon_alt=self.region_center_lat_lon_alt)
                        bs = BS(c + 1, point, self.hexagon_r, hexagon, reference_lat_lon_alt=self.region_center_lat_lon_alt)
                        self.Base_Station.append(bs)
                        c += 1
            if abs(i) % 2!= 0:
                n = int((max_y - y0 - dy) / dy + 1)  # 行数
                for j in list(range(-n, 0)) + list(range(1, n + 1)):  # 去掉 0
                    t = 0.5 if j > 0 else -0.5
                    bs_xtest = x0 + i * dx
                    bs_ytest = y0 + (j - t) * dy
                    point = [bs_xtest, bs_ytest]
                    if is_point_in_polygon(point, region_enu_array):
                        hexagon = Hexagon(*point, r=self.hexagon_r, reference_lat_lon_alt=self.region_center_lat_lon_alt)
                        bs = BS(c + 1, point, self.hexagon_r, hexagon, reference_lat_lon_alt=self.region_center_lat_lon_alt)
                        self.Base_Station.append(bs)
                        c += 1

    def find_and_store_nearby_bs_for_all_users(self):
        """
        对所有用户找到其附近的基站，并存储在用户对象的 interferece_from_bs 属性中
        """
        for bs in self.Base_Station:
            self.find_nearby_bs_for_user(bs.radius, bs.user, self.Base_Station)


    def draw_ENU_hexagons_with_users(self, ax=None, draw_bs=False,draw_hexagon=False,draw_user=False, draw_local_coordinate_system=False):
        """
        绘制所有基站的六边形，基站中心的位置标上基站的 index，以及对应基站下的用户。
        :param ax: 可选的坐标轴对象，若未传入则使用默认的当前坐标轴
        :return: 绘制所用的坐标轴对象
        """
        if ax is None:
            ax = plt.gca()
        region = self.region_enu.copy()
        region = np.array(region + [region[0]])
        ax.plot(region[:, 0], region[:, 1], 'k--')

        for bs in self.Base_Station:
            if draw_hexagon:
                # 绘制六边形
                hexagon_coords = np.array(bs.hexagon.coordinate_enu)
                ax.plot(hexagon_coords[:, 0], hexagon_coords[:, 1], 'k-')
            if draw_bs:
                # 绘制基站中心
                ax.plot(bs.coordinate_enu[0], bs.coordinate_enu[1], 'r', marker='*')
                ax.text(bs.coordinate_enu[0], bs.coordinate_enu[1], str(bs.index), fontsize=5)

            if draw_user:
                # 正确生成用户坐标数组
                user_coords = np.array([user.coordinate_enu for user in bs.user])
                # 绘制用户
                # if user_coords.size > 0:
                ax.plot(user_coords[:, 0], user_coords[:, 1], 'r.', markersize=1)

        ax.axis('equal')
        ax.set_xlabel('East (m)')
        ax.set_ylabel('North (m)')
        ax.set_title('ENU Coordinates')
        return ax

    def draw_lat_lon_hexagons_with_users(self, ax=None,draw_bs=False,draw_hexagon=False,draw_user=False, draw_local_coordinate_system=False):
        """
        绘制经纬度下所有基站的六边形，基站中心的位置标上基站的 index，以及对应基站下的用户（二维，不体现海拔）
        :param ax: 可选的坐标轴对象，若未传入则使用默认的当前坐标轴
        :return: 绘制所用的坐标轴对象
        """
        if ax is None:
            ax = plt.gca()
        region_lat_lon = [coord for coord in self.region_lat_lon_alt]
        region_lat_lon = np.array(region_lat_lon + [region_lat_lon[0]])
        ax.plot(region_lat_lon[:, 1], region_lat_lon[:, 0], 'k--')  # (经度，纬度)

        for bs in self.Base_Station:
            if draw_hexagon:
                # 获取基站六边形顶点的经纬度坐标（忽略海拔）
                hexagon_lat_lon = [coord for coord in bs.hexagon.coordinate_lat_lon_alt]
                hexagon_lat_lon = np.array(hexagon_lat_lon)
                ax.plot(hexagon_lat_lon[:, 1], hexagon_lat_lon[:, 0], 'k-')  # (经度，纬度)
            if draw_bs:
                # 绘制基站中心（经纬度坐标）并标注索引
                ax.plot(bs.coordinate_lat_lon_alt[1], bs.coordinate_lat_lon_alt[0], 'r', marker='*')
                ax.text(bs.coordinate_lat_lon_alt[1], bs.coordinate_lat_lon_alt[0], str(bs.index), fontsize=5)
            if draw_user:
                # 获取用户的经纬度坐标（忽略海拔）
                user_lat_lon = np.array([user.coordinate_lat_lon_alt for user in bs.user])
                # if user_lat_lon.size > 0:
                ax.plot(user_lat_lon[:, 1], user_lat_lon[:, 0], 'r.', markersize=1)

        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Latitude (deg)')
        ax.set_title('Latitude-Longitude Coordinates')
        return ax

    def find_nearby_bs_for_user(self, bs_r, user_object_list, BS_list):
        d1 = np.sqrt(3)*5# 指定范围
        d2 = np.inf  # 所有基站 （会遍历所有）
        d3 = np.sqrt(3) * 3
        d4 = np.sqrt(3) * 4

        for user_object in user_object_list:
            nearby_bs_list = []
            user_coordinate_enu = user_object.coordinate_enu  # 获取该用户对象的 enu 坐标
            for bs in BS_list:
                bs_coordinate_enu = bs.coordinate_enu[:2]  # 取 x,y
                user_coordinate_enu_x_y = user_coordinate_enu[:2]  # 取 x,y
                d = np.linalg.norm(np.array(bs_coordinate_enu) - np.array(user_coordinate_enu_x_y))
                if d <= d2 * bs_r:
                    nearby_bs_list.append(bs)
            user_object.interferece_from_bs = nearby_bs_list


if __name__ == "__main__":
    region = {
        'lat_min': 27,
        'lat_max': 28,
        'lon_min': 78,
        'lon_max': 79
    }
    alt = 0  # 海拔为  0 m
    region = [[region['lat_min'], region['lon_min'], alt],
              [region['lat_max'], region['lon_min'], alt],
              [region['lat_max'], region['lon_max'], alt],
              [region['lat_min'], region['lon_max'], alt],
              ]
    cluster1 = Cluster(region)
    cluster1.generate_all_BS()
    cluster1.find_and_store_nearby_bs_for_all_users()


    # 先绘制 ENU 图，但不显示
    plt.figure()
    cluster1.draw_ENU_hexagons_with_users(draw_bs=False,draw_hexagon=True,draw_user=True,draw_local_coordinate_system=True)
    fig1 = plt.gcf()  # 获取当前图形对象

    # 再绘制经纬度图，但不显示
    plt.figure()
    cluster1.draw_lat_lon_hexagons_with_users(draw_bs=True,draw_hexagon=True,draw_user=True)
    fig2 = plt.gcf()  # 获取当前图形对象

    # 同时显示两个图形
    plt.show()

