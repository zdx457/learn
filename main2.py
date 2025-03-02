from matplotlib import pyplot as plt
from Constellation import Constellation
from Create_Gp import CreateGp, calculate_N, find_closest_ground_point
from InterferenceAnalyzer import InterferenceAnalyzer
# from Simulation1 import Simulation1
from constant import R_e, h, P, S, region, beam_diameter, simulation_time
from groundpoint import GroundPoint
from k import Simulation

# 生成地面波位点对象列表
N = calculate_N(beam_diameter)
ground_points_data = CreateGp(N)
# print(ground_points_data)
ground_points = [GroundPoint(gp['index'], gp['x_ecef'], gp['y_ecef'], gp['z_ecef'], gp['lat'], gp['lon'], gp['alt'])
                 for gp in ground_points_data]
# 创建星座
constellation = Constellation(h, P, S)
# 筛选在划定区域内的地面点
selected_ground_points = [gp for gp in ground_points if
                               region['lat_min'] <= gp.lat <= region['lat_max'] and
                               region['lon_min'] <= gp.lon <= region['lon_max']]
interference_analyzer = InterferenceAnalyzer(selected_ground_points)
simulation=Simulation(constellation,ground_points,selected_ground_points,simulation_time, interference_analyzer)

SIN=simulation.run(region)
# SIN=simulation.run(region)
# print(SIN)


# # 创建干扰分析器
# interference_analyzer = InterferenceAnalyzer(ground_points)

# 创建仿真对象
# simulation = Simulation1(constellation, ground_points, simulation_time, interference_analyzer)

# # 计算区域中心点
# center_lat = (region['lat_min'] + region['lat_max']) / 2.0
# center_lon = (region['lon_min'] + region['lon_max']) / 2.0
# # 找到距离中心点最近的波位点
# closest_gp = find_closest_ground_point(simulation.selected_ground_points, center_lat, center_lon)
# user_gp_idx = closest_gp.index
#
# # 运行仿真
# simulation.run(region,user_gp_idx)
# visualization(simulation)
