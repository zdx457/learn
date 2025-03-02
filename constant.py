# constants.py
import numpy as np
import math

R_e = 6378137/1000  # km
h = 600 # km
beam_diameter = 24 # km
# 定义星座常量
P = 50
S = 50
a = R_e + h
region = {
    'lat_min': 28-1.5,
    'lat_max': 28+1.5,
    'lon_min': 79-1.5,
    'lon_max': 79+1.5
}
enlarged_region = {
    'lat_min': region['lat_min'] - 10.0,
    'lat_max': region['lat_max'] + 10.0,
    'lon_min': region['lon_min'] - 10.0,
    'lon_max': region['lon_max'] + 10.0
}
# 地球引力常数（单位：km^3/s^2）
mu = 398600.4418
DATABASE_URL = "your_database_url_here"
MAX_CONNECTIONS = 10
#仿真参数
simulation_time = 0.1  # 总仿真时间（秒）
time_step = 0.1  # 仿真步长（秒）

rx_antenna_gain = -5 # 用户天线增益
# daqi_loss = 0.1
# shanshuo_loss = 0.3
# qujihua_loss = 3
satellite_f = 5e9   #卫星频率
p1 = 44.3-30 # dBw