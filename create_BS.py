from tool import enu_to_ecef, ecef_to_lat_lon_alt, lat_lon_alt_to_ecef
import numpy as np
from constant_grounp import bs_down_angle, bs_vertical_beam_width, bs_h, bs_power, bs_radius, user_num_in_bs
from create_User import User
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class BS:
    def __init__(self, index, coordinate_enu, radius, hexagon, reference_lat_lon_alt):
        self.h = bs_h  # m 基站高度
        self.index = index
        self.coordinate_enu = coordinate_enu + [self.h]  # ENU
        self.reference_lat_lon_alt = reference_lat_lon_alt
        self.coordinate_ecef = enu_to_ecef(*self.coordinate_enu, reference_lat_lon_alt=reference_lat_lon_alt)  # ECEF
        self.coordinate_lat_lon_alt = ecef_to_lat_lon_alt(*self.coordinate_ecef)  # lat lon alt
        self.power = bs_power
        self.vertical_beam_width = bs_vertical_beam_width
        self.da = bs_down_angle
        self.radius = bs_radius
        self.hexagon = hexagon
        self.user = self.create_user_list() if user_num_in_bs > 0 else None
        self.activation = np.random.binomial(1, 0.5)  # 基站50%的激活率
        self.acceptable_interference_range = None
        self.local_coordinate_system = self.create_ecef_coordinate_system(self.coordinate_lat_lon_alt,
                                                                          self.coordinate_ecef)  # x,y,z 对应 lat,lon,alt
    # enu平面上生成用户
    def create_user_list(self):
        user_list = []
        bs_x, bs_y = self.coordinate_enu[:2]
        for _ in range(user_num_in_bs):
            user = User(bs_x, bs_y, self.radius, reference_lat_lon_alt=self.reference_lat_lon_alt)
            user_list.append(user)
        return user_list

    # 计算局部坐标系
    def create_ecef_coordinate_system(self, coordinate_lat_lon_alt, coordinate_ecef):
        """
        生成局部坐标系
        :param coordinate_lat_lon_alt: 基站lla坐标
        :param coordinate_ecef: 基站ecef坐标
        :return: 三扇区，三个坐标系
        """
        lat = coordinate_lat_lon_alt[0]
        lon = coordinate_lat_lon_alt[1]
        alt = coordinate_lat_lon_alt[2]
        # x轴方向
        lat_plus1 = lat
        lon_plus1 = lon + 0.1
        alt_plus1 = alt
        # y轴方向
        lat_plus2 = lat + 0.1
        lon_plus2 = lon
        alt_plus2 = alt

        x0, y0, z0 = np.array(coordinate_ecef)  # 基站ecef坐标
        x1, y1, z1 = np.array(lat_lon_alt_to_ecef(lat_plus1, lon_plus1, alt_plus1)) #
        x2, y2, z2 = np.array(lat_lon_alt_to_ecef(lat_plus2, lon_plus2, alt_plus2)) #

        v1 = np.array([x1 - x0, y1 - y0, z1 - z0])  # x  纬度
        v2 = np.array([x2 - x0, y2 - y0, z2 - z0])  # y  经度
        v3 = np.cross(v1, v2)  # z 海拔
        # 施密特正交化
        u1 = v1
        u2 = v2 - (np.dot(v2, u1) / np.dot(u1, u1)) * u1
        u3 = v3 - (np.dot(v3, u2) / np.dot(u2, u2)) * u2 - (np.dot(u1, v3) / np.dot(u1, u1)) * u1

        # 单位化
        u1 = u1 / np.linalg.norm(u1)  # x
        u2 = u2 / np.linalg.norm(u2)  # y
        u3 = u3 / np.linalg.norm(u3)  # z

        x = u1
        y = u2
        z = u3

        # fig = plt.figure()
        # ax = fig.add_subplot(221, projection='3d')
        #
        # # 绘制向量u1
        # ax.quiver(0, 0, 0, u1[0], u1[1], u1[2], color='r', label='u1')
        #
        # # 绘制向量u2
        # ax.quiver(0, 0, 0, u2[0], u2[1], u2[2], color='g', label='u2')
        #
        # # 绘制向量u3
        # ax.quiver(0, 0, 0, u3[0], u3[1], u3[2], color='b', label='u3')
        #
        # ax.set_xlabel('X')
        # ax.set_ylabel('Y')
        # ax.set_zlabel('Z')
        #
        # ax.legend()



        x1,y1,z1 = self.rotate_the_coordinate_system(x,y,z,0, self.da)
        # print(np.dot(x1,y1))
        # print(np.dot(x1,z1))
        # print(np.dot(z1,y1))
        # ax1 = fig.add_subplot(222, projection='3d')
        # ax1.quiver(0, 0, 0, x1[0], x1[1], x1[2], color='r', label='x1')
        # ax1.quiver(0, 0, 0, y1[0], y1[1], y1[2], color='g', label='y1')
        # ax1.quiver(0, 0, 0, z1[0], z1[1], z1[2], color='b', label='z1')
        # ax1.set_xlabel('X1')
        # ax1.set_ylabel('Y1')
        # ax1.set_zlabel('Z1')
        # ax1.legend()

        x2,y2,z2 = self.rotate_the_coordinate_system(x,y,z,120, self.da)
        # print(np.dot(x2,y2))
        # print(np.dot(x2,z2))
        # print(np.dot(z2,y2))
        # ax2 = fig.add_subplot(223, projection='3d')
        # ax2.quiver(0, 0, 0, x2[0], x2[1], x2[2], color='r', label='x2')
        # ax2.quiver(0, 0, 0, y2[0], y2[1], y2[2], color='g', label='y2')
        # ax2.quiver(0, 0, 0, z2[0], z2[1], z2[2], color='b', label='z2')
        # ax2.set_xlabel('X2')
        # ax2.set_ylabel('Y2')
        # ax2.set_zlabel('Z2')
        # ax2.legend()
        x3,y3,z3 = self.rotate_the_coordinate_system(x,y,z,240, self.da)
        # print(np.dot(x3,y3))
        # print(np.dot(x3,z3))
        # print(np.dot(z3,y3))


        # ax3 = fig.add_subplot(224, projection='3d')
        # ax3.quiver(0, 0, 0, x3[0], x3[1], x3[2], color='r', label='x3')
        # ax3.quiver(0, 0, 0, y3[0], y3[1], y3[2], color='g', label='y3')
        # ax3.quiver(0, 0, 0, z3[0], z3[1], z3[2], color='b', label='z3')
        # ax3.set_xlabel('X3')
        # ax3.set_ylabel('Y3')
        # ax3.set_zlabel('Z3')
        # ax3.legend()
        # plt.show()

        return [[x1, y1, z1],[x2,y2,z2],[x3,y3,z3]]

    def rotate_the_coordinate_system(self,x0,y0,z0,theta1,da):
        """
        逆时针旋转 theta度
        :param x: x方向
        :param y: y方向
        :param z: z方向
        :param theta1: 旋转角度
        :param da: 下倾角
        :return: 旋转后的坐标系x_da,y_da,z_da
        """
        theta1 = np.radians(theta1)
        x = np.cos(theta1) * x0 + np.sin(theta1) * y0
        y = -np.sin(theta1) * x0 + np.cos(theta1) * y0
        z = z0  # 保持不变

        x_da = np.cos(np.radians(da)) * x - np.sin(np.radians(da)) * z
        y_da = y
        z_da = np.cross(x_da, y_da)

        return [x_da, y_da, z_da]




