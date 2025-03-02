"""
上行，蜂窝用户对卫星的集总干扰
"""

from link.rain_618 import com_618_rain
from link.cloud_840 import cloud_loss
from link.atm_676 import atmosphere_loss
import numpy as np
from main_P452 import com_452
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from cluster import Cluster

from tool import *
from constant_grounp import bs_f, B
from antenna.ue_2101 import IMT_2101
from channel_h import Channel_h
import time

time1 = time.time()


a = 0.05
b = 0.05 # 1个基站

# a = 0.05
# b = 0.26 # 3个基站

# a = 0.26
# b = 0.26 # 45个基站

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
ax = fig.add_subplot(111,projection='3d')
region.append(region[0])
region_draw = np.array(region)
ax.plot(region_draw[:,0],region_draw[:,1],region_draw[:,2],color='black')



# 卫星干扰门限
target_I = -174  # dBm/Hz

# 定义纬度搜索范围
lat_min = 0
lat_max = 1


while lat_max - lat_min > 0.001: # 最多迭代12次
    # 计算中间纬度
    mid_lat = (lat_min + lat_max) / 2

    # 设置波位和卫星的经度为 0，纬度为中间纬度
    beam_lat_lon_alt = [mid_lat, 0, 0]
    satellite_lat_lon_alt = [mid_lat,0,600*1e3]  # 仰角为90°
    # 计算波位和卫星的 ECEF 坐标
    beam_ecef = lat_lon_alt_to_ecef(*beam_lat_lon_alt)
    satellite_ecef = lat_lon_alt_to_ecef(*satellite_lat_lon_alt)

    # 计算干扰功率
    I_W = 0

    distances = []
    for bs in cluster.Base_Station: # 每个基站中每个用户对波位点的干扰
        for user in bs.user:
            # 计算用户对卫星方向的增益，也就是卫星发射增益
            channel_h = Channel_h(5e9, 300, -5, beam_ecef, satellite_ecef,[-1.02321947e-05 , 2.31392886e-17 , 3.77893260e-01],
                                  user.coordinate_ecef)
            channel_gain,angle_off_axis = channel_h.calculate_channel_gain() # 包括用户发射增益，卫星天线增益，FSPL，去极化损耗
            # 云雾雨衰减
            daqi = atmosphere_loss(5e9 / 1e9, user.coordinate_lat_lon_alt[0], 5, 1)  #
            yunwu = cloud_loss(satellite_lat_lon_alt[1], satellite_lat_lon_alt[0],
                               user.coordinate_lat_lon_alt[1], user.coordinate_lat_lon_alt[0],
                               600, 5e9 / 1e9,
                               0.01)  # p是概率
            jiangshui = com_618_rain(satellite_lat_lon_alt[1], satellite_lat_lon_alt[0],
                                     user.coordinate_lat_lon_alt[1], user.coordinate_lat_lon_alt[0],
                                     600, 5e9 / 1e9)


            I_dB = user.transmitted_power + channel_gain - daqi-yunwu-jiangshui # dB
            I_W += 10 ** (I_dB / 10) # 每个用户对波位干扰集总功率 W
            distances.append(calculate_distance(user.coordinate_ecef,beam_ecef))



    I_dBm_Hz = 10 * np.log10(I_W) + 30 - 10 * np.log10(B)  # dBm/Hz
    # print(I_dBm_Hz)

    if abs(I_dBm_Hz - target_I) < 0.1:
        min_distance = min(distances)
        beam_enu = ecef_to_enu(*beam_ecef, reference_lat_lon_alt=[0,0,0])

        print("I =",I_dBm_Hz,"dBm/Hz","用户和波位的最近距离: ",min_distance,"m","波位坐标: ", beam_lat_lon_alt)
        ax.scatter(beam_enu[0],beam_enu[1],beam_enu[2],color='red')
        break

    # 根据干扰功率调整搜索范围
    if I_dBm_Hz > target_I: # 干扰太大则远离
        lat_min = mid_lat
    else:
        lat_max = mid_lat
else:
    print("未找到合适的距离")


time2 = time.time()
print("用时:" + str(time2 - time1) + "秒", "基站数量：",len(cluster.Base_Station))




