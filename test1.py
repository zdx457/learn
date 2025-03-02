"""
下行，卫星对基站用户的干扰
"""

import numpy as np
import matplotlib.pyplot as plt

from constant import satellite_f
from constant_grounp import beam_radius, B
from tool import *
from constant import R_e, mu
from cluster import Cluster
from satellite import Satellite
from channel_h import Channel_h

def calculate_up_angle(satellite_ecef, groundpoint_ecef):
    """
        函数功能：计算卫星和波位点之间的仰角
        返回：仰角alpha，单位：度
        - X, Y, Z
    """
    satellite_ecef = np.array(satellite_ecef)
    groundpoint_ecef = np.array(groundpoint_ecef) / 1000  # 统一单位为km
    d = np.linalg.norm(satellite_ecef - groundpoint_ecef)
    r = np.linalg.norm(satellite_ecef)
    # 计算 α₁
    numerator = R_e**2 + d**2 - r**2
    denominator = 2 * R_e * d
    cos_alpha = numerator / denominator
    # 防止数值误差导致 cos_alpha 超出 [-1, 1]
    cos_alpha = np.clip(cos_alpha, -1.0, 1.0)
    alpha = np.degrees(np.arccos(cos_alpha)) - 90
    return alpha

# 卫星单位都是 km
h = 600  # km
a = R_e + h  # km
n = np.sqrt(mu / a ** 3)
default_wavelength = 0.15
default_pt = 14.3  # W

orbital_elements = {
    'a': a,
    'e': 0.0,
    'i': np.radians(90),
    'Omega': 0,
    'omega': 0.0,
    'f': 0,
    'n': n
}
satellite_id = 0
satellite_init = Satellite(satellite_id, orbital_elements, default_wavelength, default_pt)

# 定义要检查的角度列表
target_angles = list(range(30, 91))
target_angles.sort(reverse=True)

# 获得符合要求的卫星对象
satellite_list = []
# 初始化卫星位置
satellite_init.update_position(0)

while True:
    satellite_ecef = satellite_init.position  # 卫星ecef坐标

    theta_check = calculate_up_angle(satellite_ecef, [R_e * 1000, 0, 0])  # 对应倾角
    # 检查 theta 是否靠近目标角度
    for target_angle in target_angles:
        i = theta_check - target_angle  # 差值
        if abs(i) <= 0.1:  # 这里假设靠近的范围是差值小于 0.1 度
            new_satellite = Satellite(satellite_id, orbital_elements, default_wavelength, default_pt)
            satellite_init.calculate_velocity(0.1)  # 计算速度
            new_satellite.velocity = satellite_init.velocity
            new_satellite.position = satellite_init.position.copy()
            satellite_list.append(new_satellite)
            target_angles.remove(target_angle)
        else:
            break

    if len(satellite_list) == 61:
        break

    # 更新卫星位置
    satellite_init.update_position(0.1)

# 定义一个范围来遍历不同的 abc 值
abc_values = np.linspace(-10, 0, 100)  # 可以根据需要调整范围和步数

# 用于保存满足条件的倾角和距离
saved_angles = []
saved_distances = []

for abc in abc_values:
    # 地面
    region = {
        'lat_min': 0.1 + abc,
        'lat_max': 0.2 + abc,
        'lon_min': -0.05,
        'lon_max': 0.05
    }
    alt = 0  # 海拔为  0 m
    region = [[region['lat_min'], region['lon_min'], alt],
              [region['lat_max'], region['lon_min'], alt],
              [region['lat_max'], region['lon_max'], alt],
              [region['lat_min'], region['lon_max'], alt],
              ]

    cluster = Cluster(region)
    cluster.generate_all_BS()

    beam_position_lla = [0, 0, 0]
    # beam_position_ecef = lat_lon_alt_to_ecef(*beam_position_lla)
    beam_position_ecef = [R_e * 1000, 0, 0]
    beam_position_enu = ecef_to_enu(*beam_position_ecef, reference_lat_lon_alt=cluster.region_center_lat_lon_alt)

    # region = cluster.region_enu.copy()
    # region = np.array(region + [region[0]])
    #
    # angle = np.linspace(0, 2 * np.pi, 100)
    # x = beam_position_enu[0] + beam_radius * np.cos(angle)
    # y = beam_position_enu[1] + beam_radius * np.sin(angle)
    # z = np.zeros_like(x)

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot(beam_position_enu[0], beam_position_enu[1], 0, 'g.', markersize=5, label='beam_center')
    # ax.plot(region[:, 0], region[:, 1], 'k--', label='Region')  # 绘制区域
    #
    # ax.plot(x, y, z, c='b', label='Beam_range')  # 波束范围

    I = 0
    N0 = -167  # dBm / Hz
    N0B_mW = 10 ** (N0 / 10) * B
    N_dBm = 10 * np.log10(N0B_mW)
    for bs in cluster.Base_Station:
        # # 绘制六边形
        # hexagon_coords = np.array(bs.hexagon.coordinate_enu)
        # ax.plot(hexagon_coords[:, 0], hexagon_coords[:, 1], hexagon_coords[:, 2], 'k-')
        # # 绘制用户
        # user_coords = np.array([user.coordinate_enu for user in bs.user])
        # ax.plot(user_coords[:, 0], user_coords[:, 1], 0, 'r.', markersize=1)
        for user in bs.user:  # 目前是一个基站一个用户
            for satellite in satellite_list:  # 不同倾角的卫星
                channel_h = Channel_h(satellite_f, 300, -5, beam_position_ecef, satellite.position * 1000,
                                      satellite.velocity, user.coordinate_ecef)
                channel_gain = channel_h.calculate_channel_gain_from_table()
                # 卫星发射功率 + 信道增益(包括发送天线增益，接收天线增益，自由传播损耗，大气，闪烁，去极化)
                I = 10 * np.log10(satellite.pt) + channel_gain  # dB
                I_dBm = I + 30
                INR = I_dBm - N_dBm
                distance = np.linalg.norm(np.array(user.coordinate_ecef) - np.array(beam_position_ecef))

                # 检查 I/N 是否接近 -6
                if abs(INR - (-6)) <= 0.1:  # 可以根据需要调整接近的范围
                    theta_check = calculate_up_angle(satellite.position * 1000, beam_position_ecef)
                    saved_angles.append(theta_check)
                    saved_distances.append(distance)

            break

    # ax.axis('equal')
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # ax.set_zlabel("Z")
    # ax.legend()
    # plt.close(fig)  # 关闭当前的图形，避免过多图形窗口

# 保存满足条件的倾角和距离到文件
with open('saved_data.txt', 'w') as file:
    for angle, distance in zip(saved_angles, saved_distances):
        file.write(f"倾角: {angle} 度, 距离: {distance} m\n")

print("满足条件的倾角和距离已保存到 saved_data.txt")