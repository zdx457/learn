"""
下行，卫星对基站用户的干扰
"""

import numpy as np
import matplotlib.pyplot as plt

from link.rain_618 import com_618_rain
from link.cloud_840 import cloud_loss
from link.atm_676 import atmosphere_loss
from constant import satellite_f, R_e, mu # 卫星参数
from constant_grounp import beam_radius, B # 地面参数
from tool import *
from cluster import Cluster
from satellite import Satellite
from channel_h import Channel_h
from create_User import User
import pandas as pd

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
a = R_e  + h # km
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


beam_position_lla = [0,0,0] # 波位lla坐标
# beam_position_ecef = lat_lon_alt_to_ecef(*beam_position_lla)
beam_position_ecef = [R_e * 1000,0,0] # 波位ecef坐标
beam_position_enu = ecef_to_enu(*beam_position_ecef, reference_lat_lon_alt=[0,0,0])

# 定义要检查的角度列表
target_angles = list(range(0,91))
target_angles.sort(reverse=True)

# angles_2_to_6 = list(range(2, 7))
# angles_7_to_15 = [i * 0.05 + 7 for i in range(int((15 - 7) / 0.05 + 1))]
#
# angles_16_to_91 = list(range(16, 91))
# target_angles = angles_2_to_6 + angles_7_to_15 + angles_16_to_91
# target_angles.sort(reverse=True)

# 获得符合要求的卫星对象
satellite_list = []

# 初始化卫星位置
satellite_init.update_position(0)


while True:

    satellite_ecef = satellite_init.position  # 卫星ecef坐标

    theta_check = calculate_up_angle(satellite_ecef, beam_position_ecef)  # 对应倾角
    # 检查 theta 是否靠近目标角度
    for target_angle in target_angles:
        i = theta_check - target_angle  # 差值
        if abs(i) <= 0.02:  # 这里假设靠近的范围是差值小于 0.1 度
            new_satellite = Satellite(satellite_id, orbital_elements, default_wavelength, default_pt)
            satellite_init.calculate_velocity(0.05)  # 计算速度
            new_satellite.velocity = satellite_init.velocity
            new_satellite.position = satellite_init.position.copy()
            satellite_list.append([new_satellite, theta_check])

            target_angles.remove(target_angle)
        else:
            break

    if len(target_angles) == 0:
        break

    # 更新卫星位置
    satellite_init.update_position(0.05)





data = []

threshold = -173 # dBm/hz
for satellite, angle in satellite_list:
    left, right = -5, 0
    closest_distance = None
    closest_angle = None
    closest_I_dBm_Hz = None
    closest_user_lla = None
    closest_angle_off_axis = None
    min_diff = float('inf')

    while right - left > 0.001:  # 最大找12次
        mid = (left + right) / 2
        user_lla = [mid, 0, 0]  # 用户lla坐标
        user_ecef = lat_lon_alt_to_ecef(*user_lla) # 用户ecef坐标
        # user_enu = ecef_to_enu(*user_ecef, reference_lat_lon_alt=[0, 0, 0])
        # user = User(user_enu[0], user_enu[1], 0, beam_position_lla)
        channel_h = Channel_h(satellite_f, 300, -5, beam_position_ecef, satellite.position * 1000,
                              satellite.velocity, user_ecef)
        # channel_gain = channel_h.calculate_channel_gain_from_table() # 查表
        channel_gain,angle_off_axis = channel_h.calculate_channel_gain()
        satellite_lla = ecef_to_lat_lon_alt(*satellite.position)

        daqi = atmosphere_loss(satellite_f/1e9,user_lla[0],5, 1) #
        yunwu = cloud_loss(satellite_lla[1],satellite_lla[0],user_lla[1],user_lla[0],h,satellite_f/1e9,0.01) # p是概率
        jiangshui = com_618_rain(satellite_lla[1],satellite_lla[0],user_lla[1],user_lla[0],h,satellite_f/1e9)

        # 卫星发射功率 + 信道增益(包括发送天线增益，接收天线增益，自由传播损耗，去极化) -大气损耗-云损耗- 阴影损耗
        I_dBm_Hz = 10 * np.log10(satellite_init.pt) + channel_gain -daqi-yunwu-jiangshui + 30 - 10 * np.log10(B) # dBm/hz
        # I_dBm_Hz = 10 * np.log10(satellite_init.pt) + channel_gain  + 30 - 10 * np.log10(B)# dBm/hz


        tolerance = 0.1

        diff = abs(I_dBm_Hz - threshold) # dBm/hz
        # 保存最近的结果
        if diff < min_diff:
            min_diff = diff
            closest_distance = np.linalg.norm(np.array(user_ecef) - np.array(beam_position_ecef))
            closest_angle = angle
            closest_I_dBm_Hz = I_dBm_Hz
            closest_user_lla = user_lla
            closest_angle_off_axis = angle_off_axis

        if diff <= tolerance:  # 干扰和门限的差值
            distance = np.linalg.norm(np.array(user_ecef) - np.array(beam_position_ecef))  # 波束和用户的距离
            angle = round(angle, 3)
            I_dBm_Hz = round(I_dBm_Hz, 3)
            distance = round(distance, 3)
            print("仰角:", angle, "度，\t对应I=", I_dBm_Hz, "\t保护距离:", distance, "m", "\t用户lla:", user_lla)
            print("离轴角:", angle_off_axis, "度")
            data.append([angle, I_dBm_Hz, distance])
            break
        elif I_dBm_Hz < threshold:  # 干扰太大，远离波束
            left = mid
        else:
            right = mid
    else:
        # 如果二分法没有找到满足条件的值，使用最接近的值
        closest_angle = round(closest_angle, 3)
        closest_I_dBm_Hz = round(closest_I_dBm_Hz, 3)
        closest_distance = round(closest_distance, 3)
        print("仰角:", closest_angle, "度，\t对应I=", closest_I_dBm_Hz, "\t近似保护距离:", closest_distance, "m","\t用户lla:", closest_user_lla)
        print("离轴角:", closest_angle_off_axis, "度")
        data.append([closest_angle, closest_I_dBm_Hz, closest_distance])

df = pd.DataFrame(data, columns=['angle', 'INR', 'distance'])
df.to_excel("倾角对应距离.xlsx", index=False, sheet_name='Sheet1')
