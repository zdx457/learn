from tool import calculate_angle
from tool import calculate_distance, ecef_to_enu
import numpy as np
from antenna.ue_2101 import IMT_2101
from channel_h import Channel_h
from sg_downlink_single import SGDownlink
import sys
# 将文件所在目录添加到 Python 解释器的搜索路径中
sys.path.append(r"D:\python projects\pythonProject\Py452-main\src")
from main_P452 import com_452



def calcultion_grounp_interference(f,rx_antenna_gain, grounp_point_ecef, user_object):
    I = 0.0
    for bs in user_object.interferece_from_bs:
        distance = calculate_distance(bs.coordinate_ecef,grounp_point_ecef)  # 距离
        x,y,z = bs.local_coordinate_system # 三个轴的方向
        theta, phi = calculate_angle(x,y,z,bs.coordinate_ecef,user_object.coordinate_ecef) # 这里 x 是经度方向（阵列方向），y 是纬度方向，z 是海拔，原点是基站坐标，目标是用户坐标
        tx_antenna_gain = IMT_2101(theta, phi) # 增益
        FSPL = 32.45 + 20 * np.log10(distance / 1e3) + 20 * np.log10(f / 1e6)
        remove_polar = 3 # dB 去极化损耗
        I += 10**((bs.power +tx_antenna_gain + rx_antenna_gain- FSPL-remove_polar)/10)  * bs.activation # 对目标的干扰
    return I

