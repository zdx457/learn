import copy
import numpy as np
from Link_Establishment import establish_links
from constant import enlarged_region, time_step, beam_diameter,f,p1
# from ephemeris import generate_ephemeris
from cluster import Cluster
from channel_h import Channel_h
from tool import calculate_distance,ecef_to_enu
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

class Simulation:
    def __init__(self, constellation, ground_points,selected_ground_points, simulation_time, interference_analyzer):
        self.ephemeris = []
        self.constellation = constellation       #星座对象
        self.ground_points = ground_points       #地面波位
        self.simulation_time = simulation_time   #仿真总时长
        self.time_step = time_step               #仿真步长
        self.all_links = []
        self.interference_analyzer=interference_analyzer
        self.cluster=None
        self.selected_ground_points = selected_ground_points
        self.selected_ground_points_dict = {gp.index:gp for gp in self.selected_ground_points}
        self.generate_ephemeris()

    def generate_ephemeris(self):
        """
        在仿真开始前，生成半小时内、每个时隙0.1s的卫星状态快照。
        """
        # 1. 先更新 constellation 到 t=0
        self.constellation.update_positions(0)
        self.constellation.update_velocity()
        total_small_steps = int(self.simulation_time / self.time_step)  # 1800 / 0.1 = 18000
        for step_idx in range(total_small_steps):
            # 深拷贝当前星座对象
            constellation_copy = copy.deepcopy(self.constellation)
            self.ephemeris.append(constellation_copy)
            # 更新星座到下一个时隙
            self.constellation.update_positions(self.time_step)
            self.constellation.update_velocity()

    def run(self, region):
        # 构建一个字典方便 ID->卫星对象
        id_to_sat = {sat.satellite_id: sat for sat in self.constellation.satellites}
        alt = 0  # 海拔为  0 m
        region1 = [[region['lat_min'], region['lon_min'], alt],
                  [region['lat_max'], region['lon_min'], alt],
                  [region['lat_max'], region['lon_max'], alt],
                  [region['lat_min'], region['lon_max'], alt],
                  ]
        self.cluster=Cluster(region1) # 区域内基站信息 m
        self.cluster.generate_all_BS()
        total_steps = int(self.simulation_time / self.time_step)
        current_global_time = 0.0

        SIN=[]
        for step in range(total_steps):
            satellites_at_small_step = self.ephemeris[step].satellites # 对应时隙星座的卫星
            current_global_time += self.time_step
            # 筛选扩大区域内卫星
            selected_satellites = [sat for sat in self.constellation.satellites if
                                   enlarged_region['lat_min'] <= sat.lat <= enlarged_region['lat_max'] and
                                   enlarged_region['lon_min'] <= sat.lon <= enlarged_region['lon_max']]
            # 存储建链关系，格式[(波位id,卫星id),()...]
            links = establish_links(self.selected_ground_points, selected_satellites)
            # 转换为satellite_links_dict,格式：[卫星i：（波位1，波位2...),卫星j：（波位1，波位2...)...]
            satellite_links_dict = {}
            for (gp_idx, sat_id) in links:
                if sat_id not in satellite_links_dict:
                    satellite_links_dict[sat_id] = []
                satellite_links_dict[sat_id].append(gp_idx)

            # 卫星下行
            BS_information = self.cluster.Base_Station  # 基站信息

            beam_link_satellite = {}
            for i in links:  # 把波位对应关系变为字典
                beam_link_satellite[i[0]] = i[1]  # {波位序号:卫星序号}

            for bs in BS_information:
                # 找基站位于的波位点序号
                distances = []
                for gp in self.selected_ground_points:
                    gp_coord_ecef = [gp.x_ecef * 1000, gp.y_ecef * 1000, gp.z_ecef * 1000]  # m
                    d = calculate_distance(bs.coordinate_ecef,gp_coord_ecef)
                    distances.append((d,gp.index))
                # 找出最小距离及其对应的序号
                min_distance, k = min(distances, key=lambda x: x[0])  # 把distances中的每个元组的第一个值作为参考
                # print(min_index)
                # 对于第k个波位点的干扰波位点
                self.interference_analyzer.compute_interference_map(k)
                interfering_indices=self.interference_analyzer.get_interfering_beams(k)
                # print(interfering_indices)

                threshold = 10 * beam_diameter * 1000  # 卫星波束半径 m
                S = 0
                I = 0
                interfering_indices.append(k)
                for gp_index in interfering_indices:#获取对应波位索引
                    gp = self.selected_ground_points_dict[gp_index] # 获取对应索引的gp对象
                    gp_coord_ecef = [gp.x_ecef * 1000, gp.y_ecef * 1000, gp.z_ecef * 1000]  # m
                    if calculate_distance(bs.coordinate_ecef, gp_coord_ecef) < threshold:
                        """
                        在阈值内的波位点
                        找对应的卫星
                        根据卫星和基站的位置算出角度，用Channel_h
                        查表得到卫星的增益
                        卫星发射功率 × 信道增益 = 对目标的干扰
                        """

                        satellite_index = beam_link_satellite[gp.index]  # 卫星索引
                        channel = Channel_h(f,  # 频率
                                            300,  # 温度
                                            -5,  # 接收端天线增益
                                            gp_coord_ecef,  # 波束中心  m
                                            np.array([i * 1000 for i in satellites_at_small_step[satellite_index].position]),  # 卫星坐标  m
                                            satellites_at_small_step[satellite_index].velocity,  # 卫星运动方向
                                            np.array(bs.coordinate_ecef)  # 用户坐标（这里用基站坐标代替） m
                                            )
                        channel_h = channel.calculate_channel_gain_from_table()  # 信道增益

                        # 卫星发射功率 × 信道增益 = 对目标的干扰
                        if gp.index != k:
                            I += p1 * 10 ** (channel_h / 10)  # 干扰
                        else:
                            S += p1 * 10 ** (channel_h / 10)  # 有用信号
                SIN.append(S/I)
                print(SIN)
                # break

                # 先绘制ENU图，并显示
                plt.figure()
                self.cluster.draw_ENU_hexagons_with_users()
                fig1 = plt.gcf()  # 获取当前图形对象
                # 获取 interfering_indices 的坐标
                interfering = [self.selected_ground_points_dict[index] for index in interfering_indices]
                interfering_ecef = [[i.x_ecef*1000, i.y_ecef*1000, i.z_ecef*1000] for i in interfering]
                interfering_enu = [ecef_to_enu(*j, reference_lat_lon_alt=self.cluster.region_center_lat_lon_alt) for j
                                   in interfering_ecef]

                x = [point[0] for point in interfering_enu]  # 提取 x 坐标列表
                y = [point[1] for point in interfering_enu]  # 提取 y 坐标列表
                plt.scatter(x, y, color='red', label='Interfering gp')
                plt.scatter(x[-1], y[-1], color='black', label='bs in gp')
                plt.scatter(bs.coordinate_enu[0],bs.coordinate_enu[1],color='green',label='bs')
                s_coord_ecef = [i * 1000 for i in satellites_at_small_step[beam_link_satellite[k]].position]
                s_coord_enu=ecef_to_enu(*s_coord_ecef, reference_lat_lon_alt=self.cluster.region_center_lat_lon_alt)
                plt.scatter(s_coord_enu[0],s_coord_enu[1],color='yellow',label='s')
                plt.plot([x[-1],s_coord_enu[0]],[y[-1],s_coord_enu[1]])
                for i, j in zip(x, y):
                    # 创建一个圆补丁
                    circle = Circle((i, j), radius=beam_diameter * 1000, fill=True, edgecolor='blue', facecolor='blue',
                                    alpha=0.3)
                    # 将圆添加到轴上
                    plt.gca().add_patch(circle)
                # 添加图例
                plt.legend()

                # 再绘制经纬度图，并显示
                plt.figure()
                self.cluster.draw_lat_lon_hexagons_with_users()
                fig2 = plt.gcf()  # 获取当前图形对象
                plt.show()
                break
        return  SIN






