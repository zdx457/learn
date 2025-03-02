# interference_analyzer.py

import numpy as np
from constant import R_e,  beam_diameter


class InterferenceAnalyzer:
    def __init__(self, ground_points, interference_threshold_factor=8):
        """
        初始化干扰分析器。
        参数：
        - ground_points: 地面波位点对象列表
        - interference_threshold_factor: 干扰距离阈值的倍数
        """
        self.ground_points = ground_points
        self.interference_threshold_factor = interference_threshold_factor
        self.r = beam_diameter
        self.threshold = self.interference_threshold_factor * self.r
        self.interference_map = {}

    # def compute_average_distance(self):
    #     """
    #     计算地面波位点间的平均距离的一半 r。
    #     返回：
    #     - r: 平均距离的一半
    #     """
    #     N = len(self.ground_points)
    #     A = 4 * np.pi * R_e ** 2 / N
    #     r = np.sqrt(A) / 2
    #     return r

    def compute_interference_map(self, selected_ground_index):
        """
        计算指定波位点的干扰邻居列表。

        参数：
        - selected_ground_index: 选定的波位点索引
        """
        selected_gp = next((gp for gp in self.ground_points if gp.index == selected_ground_index), None)
        if not selected_gp:
            self.interference_map[selected_ground_index] = []
            return
        selected_position = np.array([selected_gp.x_ecef, selected_gp.y_ecef, selected_gp.z_ecef])
        interfering_indices = []
        for gp in self.ground_points:
            if gp.index == selected_ground_index:
                continue
            position = np.array([gp.x_ecef, gp.y_ecef, gp.z_ecef])
            distance = np.linalg.norm(position - selected_position)
            if distance <= self.threshold:
                interfering_indices.append(gp.index)
        self.interference_map[selected_ground_index] = interfering_indices

    def get_interfering_beams(self, selected_ground_index):
        """
        获取指定波位点的干扰邻居列表。

        参数：
        - selected_ground_index: 选定的波位点索引

        返回：
        - 干扰波位点索引列表
        """
        return self.interference_map.get(selected_ground_index, [])

    def find_interfering_links(self, simulation, selected_ground_index):
        """
        对于指定的波位点，确定每个小时间隙可能产生干扰的卫星-波位点对。

        参数：
        - simulation: 仿真对象
        - selected_ground_index: 选定的波位点索引

        返回：
        - interference_links: 列表，每个元素为一个字典，包含:
          {
            'big_time_step': 大时隙编号,
            'hop': 跳波束编号,
            'satellite_id': 卫星编号,
            'ground_point_id': 波位点编号
          }
        """
        interference_links = []
        # 获取干扰波位点列表
        interfering_indices = self.get_interfering_beams(selected_ground_index)

        # 遍历所有小时间隙的结果
        small_step_results = simulation.get_small_step_results()
        for step_idx, step_res in enumerate(small_step_results):
            time = step_res['time']
            active_links = step_res['active_links']

            # 遍历所有卫星的激活链路
            for sat_id, gp_list in active_links.items():
                for gp_idx in gp_list:
                    if gp_idx in interfering_indices:
                        # 找到干扰链路
                        interference_links.append({
                            'big_time_step': step_idx // m,  # 假设step_idx对应于大时隙
                            'hop': (step_idx % m) // (m // n),  # 计算所属跳波束
                            'satellite_id': sat_id,
                            'ground_point_id': gp_idx
                        })

        return interference_links

