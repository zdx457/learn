import copy
import numpy as np
from Link_Establishment import establish_links
from constant import enlarged_region, time_step, beam_diameter,f,p1
from constant_grounp import bs_radius,user_num_in_bs
# from ephemeris import generate_ephemeris
from cluster import Cluster
from channel_h import Channel_h
from tool import calculate_distance,ecef_to_enu,lat_lon_alt_to_ecef,ecef_to_lat_lon_alt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from calcultion_grounp_interference import calcultion_grounp_interference
from create_User import User



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
        region1 = {
            'lat_min': 27.5,
            'lat_max': 28.5,
            'lon_min': 78.5,
            'lon_max': 79.5
        }
        region1 = [[region1['lat_min'], region1['lon_min'], alt],
                  [region1['lat_max'], region1['lon_min'], alt],
                  [region1['lat_max'], region1['lon_max'], alt],
                  [region1['lat_min'], region1['lon_max'], alt],
                  ]
        self.cluster=Cluster(region1) # 区域内基站信息 m
        # self.cluster=Cluster(region1, 1000.0000009639708,2 ) # 区域内基站信息 m
        self.cluster.generate_all_BS() # 生成基站
        self.cluster.find_and_store_nearby_bs_for_all_users() # 生成用户附近基站

        total_steps = int(self.simulation_time / self.time_step)
        current_global_time = 0.0


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
            # 获得波位点{index:[enu_x,enu_y]}
            gp_center_user_enu_dict = {i.index:ecef_to_enu(i.x_ecef*1000,i.y_ecef*1000,i.z_ecef*1000,self.cluster.region_center_lat_lon_alt)[:2] for i in self.selected_ground_points}
            # 创建新用户列表
            new_user_list = []
            for index in gp_center_user_enu_dict:
                x,y = gp_center_user_enu_dict[index]
                user = User(x, y, 0, self.cluster.region_center_lat_lon_alt) #bs_r=0表示用户位置就在[x,y]处
                user.beam_index=index
                new_user_list.append(user)
            # 为新用户找到附近的基站并填充 interferece_from_bs 属性
            self.cluster.find_nearby_bs_for_user(self.cluster.hexagon_r, new_user_list, BS_information)

            # 先绘制ENU图，并显示
            plt.figure()

            self.cluster.draw_ENU_hexagons_with_users(draw_bs=False,draw_hexagon=True,draw_user=False,draw_local_coordinate_system=False)
            fig1 = plt.gcf()  # 获取当前图形对象



            # 对每个波位计算地面干扰
            SIN = {}
            for user in new_user_list:
                index=user.beam_index
                ground_point = [self.selected_ground_points_dict[index].x_ecef*1000,
                                self.selected_ground_points_dict[index].y_ecef*1000,
                                self.selected_ground_points_dict[index].z_ecef*1000] # m
                # 波位对应的卫星
                satellite= satellites_at_small_step[beam_link_satellite[index]] # 对象
                satellite_coord_ecef=[i*1000 for i in satellite.position] # m

                channel = Channel_h(f,300,-5,ground_point,np.array(satellite_coord_ecef),satellite.velocity,np.array(user.coordinate_ecef))
                channel_h = channel.calculate_channel_gain_from_table() # dB
                C = 10**((p1 + channel_h) / 10) # W
                I = calcultion_grounp_interference(f,-5,ground_point,user) # W
                C_I_dB = 10 * np.log10(C / I)
                SIN[index] = 10 * np.log10(C / I)

                # 波位区域
                region_lat_lon_alt = [[region['lat_min'], region['lon_min'], alt],
                           [region['lat_max'], region['lon_min'], alt],
                           [region['lat_max'], region['lon_max'], alt],
                           [region['lat_min'], region['lon_max'], alt],
                           ]
                region_ecef = [lat_lon_alt_to_ecef(*i) for i in region_lat_lon_alt]
                region_enu = [ecef_to_enu(*i,reference_lat_lon_alt=self.cluster.region_center_lat_lon_alt) for i in region_ecef]
                region_array = np.array(region_enu + [region_enu[0]])
                plt.plot(region_array[:, 0], region_array[:, 1], 'k--')


                # 画图
                user_enu = user.coordinate_enu
                user_x = user_enu[0]
                user_y = user_enu[1]
                plt.scatter(user_x,user_y,color='green',label='user',alpha=0.5)
                text = "{:.3f}".format(float(C_I_dB))
                plt.text(user_x, user_y, text, fontsize=10)
                # plt.text(user_x,user_y,"{:.3f}".format(float(C/I)),fontsize=10)
                # plt.text(user_x,user_y,str(I),fontsize=10)
                circle = Circle((user_x, user_y), radius=beam_diameter*1000, fill=False, edgecolor='blue', facecolor=None,
                                alpha=0.3)
                plt.gca().add_patch(circle)
                # plt.legend()

            print(SIN)
            plt.savefig('ENU_inf.png')
            plt.show()
            return  SIN






