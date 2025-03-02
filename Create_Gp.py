import numpy as np
from constant import R_e
from math import ceil, sqrt, pi, cos, sin
from tool import ecef_to_lat_lon_alt

def calculate_N(beam_diameter, coverage_radius=R_e):
    """
    计算波位点数量 N，使得波束几乎相切。

    参数：
    - beam_diameter (float): 波束直径（单位：米）
    - coverage_radius (float): 覆盖区域半径（单位：米），默认地球半径

    返回：
    - N (int): 波位点数量
    """
    r = beam_diameter   # 波束半径
    A_hex = (3 * sqrt(3) / 2) * r**2  # 六边形单元格面积
    A_coverage = 4 * pi * coverage_radius**2  # 地球表面面积
    packing_density = 0.9069  # 六边形排列填充密度
    N = ceil((A_coverage / A_hex) * packing_density)
    return N

def CreateGp(N):
    # 地球半径（单位：公里）
    # 需要生成的地面波位数量
    # 黄金角度（以弧度表示）
    golden_angle = np.pi * (3 - np.sqrt(5))

    # 创建列表，存储每个点的ECEF坐标
    ground_stations = []

    for n in range(1, N + 1):
        # 计算z坐标
        z = (2 * n - 1) / N - 1

        # 计算半径（在单位球面上，半径为sqrt(1 - z^2)）
        radius = np.sqrt(1 - z * z)

        # 计算角度theta
        theta = golden_angle * n

        # 计算x和y坐标
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        # 将单位球面的坐标转换为地球半径的坐标
        x_ecef = R_e * x
        y_ecef = R_e * y
        z_ecef = R_e * z

        lla = ecef_to_lat_lon_alt(x_ecef,y_ecef,z_ecef)
        lat_deg = lla[0]
        lon_deg = lla[1]
        alt = lla[2]
        # # 计算经纬度
        # lon = np.arctan2(y, x)
        # hyp = np.sqrt(x * x + y * y)
        # lat = np.arctan2(z, hyp)
        #
        # # 转换为度
        # lat_deg = np.degrees(lat)
        # lon_deg = np.degrees(lon)

        # 将坐标存储到列表中
        station = {
            'index': n - 1,  # 索引从0开始
            'x_ecef': x_ecef,
            'y_ecef': y_ecef,
            'z_ecef': z_ecef,
            'lat': lat_deg,
            'lon': lon_deg,
            'alt': alt
        }
        ground_stations.append(station)

    return ground_stations

def find_closest_ground_point(selected_ground_points, center_lat, center_lon):
    """
    在选定的地面波位点中，选择离中心点最近的波位点（通过纬度和经度差异）。

    参数：
    - selected_ground_points (list of GroundPoint): 选定的地面波位点列表
    - center_lat (float): 区域中心点的纬度
    - center_lon (float): 区域中心点的经度

    返回：
    - closest_gp (GroundPoint): 离中心点最近的波位点
    """
    min_diff = float('inf')
    closest_gp = None
    for gp in selected_ground_points:
        lat_diff = abs(gp.lat - center_lat)
        lon_diff = abs(gp.lon - center_lon)
        # 可以使用欧几里得距离或其他组合方式
        total_diff = lat_diff + lon_diff  # 简单相加
        if total_diff < min_diff:
            min_diff = total_diff
            closest_gp = gp
    return closest_gp