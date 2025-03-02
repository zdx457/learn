import numpy as np

bs_h = 30
bs_power= 10.77 # 40.77dBm-30=10.77dBw
bs_vertical_beam_width = 65 # 65度，波束垂直宽度

# bs_da = 34.218358 # 下倾角34.2°，半径为1000m
# bs_da = 33.07293869768 # 半径3000m
bs_down_angle = 32.84377055  # 下倾角33.84°，半径为5000m
bs_radius = bs_h / np.tan(np.radians(bs_down_angle)-np.radians(bs_vertical_beam_width)/2) # 覆盖半径h/tan(theta-alpha/2), alpha=60是垂直波束宽度，theta是下倾角
user_num_in_bs = 1 # 基站用户数

bs_f = 2 # GHz
k = 1.380649e-23 # J/K
T = 300 # K

B = 5e6 # Hz
user_receive_antenna_gain = -5 # dBi

beam_radius = 24000 # m
satellite_height = 600000 # m