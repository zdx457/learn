
import numpy as np
from main_P452 import com_452
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from cluster import Cluster
from create_User import User
from tool import *
from constant_grounp import bs_f, user_receive_antenna_gain,bs_h,B
from antenna.ue_2101 import IMT_2101
from channel_h import Channel_h
import time

time1 = time.time()


# a = 0.05
# b = 0.05 # 1个基站，用时 1.8 s 距离:  8755.038841846148 m 用户纬度经度海拔坐标:  [0.0787353515625, 0, 0]

a = 0.05
b = 0.26 # 3个基站，用时 64 s 距离:  9162.245207098533 m 用户纬度经度海拔坐标:  [0.0823974609375, 0, 0]

# a = 0.26
# b = 0.26 # 45个基站，用时: 852.38秒，距离: 10125.66645424675 m 用户纬度经度海拔坐标 [0.32470703125, 0, 0]

# a = 0.4 #
# b = 0.4 # 115个基站

region = {
    'lat_min':  -a,
    'lat_max':  a,
    'lon_min':  -b,
    'lon_max':  b
}
alt = 0  # 海拔为  0 m
region = [[region['lat_min'], region['lon_min'], alt],
          [region['lat_max'], region['lon_min'], alt],
          [region['lat_max'], region['lon_max'], alt],
          [region['lat_min'], region['lon_max'], alt],
          ]

cluster = Cluster(region)
cluster.generate_all_BS()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
for bs in cluster.Base_Station:
    # 绘制六边形
    hexagon_coords = np.array(bs.hexagon.coordinate_enu)
    ax.plot(hexagon_coords[: , 0], hexagon_coords[: , 1], 'k-')

# ax.axis('equal')
# ax.set_xlabel('East (m)')
# ax.set_ylabel('North (m)')
# ax.set_title('ENU Coordinates')
# ax.legend()
# plt.show()



# 定义目标干扰功率
target_I = -174

# 定义纬度搜索范围
lat_min = 0
lat_max = 5


while lat_max - lat_min > 0.001:
    # 计算中间纬度
    mid_lat = (lat_min + lat_max) / 2

    # 设置用户的经度为 0，纬度为中间纬度
    user_lat_lon_alt = [mid_lat, 0, 0]

    # 计算用户的 ECEF 坐标和 ENU 坐标
    user_ecef = lat_lon_alt_to_ecef(*user_lat_lon_alt)

    # 计算干扰功率
    I_W = 0

    distances = []
    for bs in cluster.Base_Station:
        # 计算干扰
        vertice = np.array(user_ecef) - np.array(bs.coordinate_ecef)  # 基站指向用户的 ecef 方向

        # 对三个扇区计算总增益和路损
        tx_antenna_gain_W = 0  # 基站发射增益 W

        for i in bs.local_coordinate_system:
            x1, y1, z1 = i  # 三个轴的方向
            theta1, phi1 = calculate_angle(x1,y1,z1,bs.coordinate_ecef,user_ecef)# 这里 x 是经度方向（阵列方向），y 是纬度方向，z 是海拔，原点是基站坐标，目标是用户坐标
            tx_dB = IMT_2101(theta1, phi1)
            tx_antenna_gain_W += 10 ** (tx_dB / 10)  # 基站发射增益 W

        tx_antenna_gain_dB = 10 * np.log10(tx_antenna_gain_W)  # 基站发射增益 dB
        # 频率（GHz），概率（默认 50%），极化方式（1 水平极化），
        # 发送经度，发送维度，接收经度，接收纬度，发送天线增益，接收天线增益，发射天线高度，接收天线高度
        Lb = com_452(bs_f, 50, 1,
                     bs.coordinate_lat_lon_alt[1], bs.coordinate_lat_lon_alt[0],
                     user_lat_lon_alt[1], user_lat_lon_alt[0],
                     tx_antenna_gain_dB, user_receive_antenna_gain,
                     bs_h, 1.5)

        I_W += 10 ** ((bs.power + tx_antenna_gain_dB + user_receive_antenna_gain - Lb - 3) / 10)  # 对目标的干扰 W
        distance = calculate_distance(user_ecef, bs.coordinate_ecef)
        distances.append(distance)


    I_dBm_Hz = 10 * np.log10(I_W) + 30 - 10 * np.log10(B)  # dBm/Hz

    if abs(I_dBm_Hz - target_I) < 0.1:
        min_distance = min(distances)
        user_enu = ecef_to_enu(*user_ecef, reference_lat_lon_alt=cluster.region_center_lat_lon_alt)  # 用户ENU坐标
        ax.scatter(user_enu[0], user_enu[1], 0, s=5, c='red', label='User in beam center')
        print("I =",I_dBm_Hz,"dBm/Hz","距离: ",min_distance,"m","用户纬度经度海拔坐标: ", user_lat_lon_alt)
        break

    # 根据干扰功率调整搜索范围
    if I_dBm_Hz > target_I: # 干扰太大则远离
        lat_min = mid_lat
    else:
        lat_max = mid_lat


time2 = time.time()
print("用时:" + str(time2 - time1) + "秒", "基站数量：",len(cluster.Base_Station))

ax.axis('equal')
ax.set_xlabel('East (m)')
ax.set_ylabel('North (m)')
ax.set_title('ENU Coordinates')
ax.legend()
plt.show()


