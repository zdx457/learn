import numpy as np
import pyproj


def calculate_angle(x,y,z,origin_point,aim_point):
    """
基站计算俯仰角和方位角
    俯仰角 theta
    方位角 phi
    """

    vertice = np.array(aim_point) - np.array(origin_point)
    # 俯仰角 (0,180)
    theta = np.arccos(np.dot(vertice, z) / (
            np.linalg.norm(vertice) * np.linalg.norm(z)))  # 弧度
    # 方位角 (-90,90)
    if np.dot(y, vertice) > 0:
        vertice_Proj_toXOY = (vertice) - np.dot(vertice, z) * z  # XOY面上投影
        k = np.dot(vertice_Proj_toXOY, x)
        m1 = np.linalg.norm(vertice_Proj_toXOY) * np.linalg.norm(x)
        cos_phi = k / m1
        cos_phi = np.clip(cos_phi, -1.0, 1.0)
        phi = np.arccos(cos_phi)
    else:
        vertice_Proj_toXOY = (vertice) - np.dot(vertice, z) * z  # XOY面上投影
        k = np.dot(vertice_Proj_toXOY, x)
        m1 = np.linalg.norm(vertice_Proj_toXOY) * np.linalg.norm(x)
        cos_phi = k / m1
        cos_phi = np.clip(cos_phi, -1.0, 1.0)
        phi = - np.arccos(cos_phi)
    # 弧度转角度
    theta = np.degrees(theta)
    phi = np.degrees(phi)

    return theta, phi  # 弧度转角度


def calculate_center(region):
    region = np.array(region)
    x_coords = region[:, 0]
    y_coords = region[:, 1]
    h_coords = region[:, 2]
    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)
    min_h, max_h = min(h_coords), max(h_coords)
    x_center = (min_x + max_x) / 2
    y_center = (min_y + max_y) / 2
    h_center = (min_h + max_h) / 2
    return [x_center, y_center, h_center]


def calculate_distance(p1, p2):
    """
    计算三维空间中两点之间的距离。

    参数:
    p1: 第一个点的坐标，格式为(x1, y1, z1)
    p2: 第二个点的坐标，格式为(x2, y2, z2)

    返回:
    两点之间的距离
    """
    vector_subtraction = np.array(p1) - np.array(p2)  # 向量相减
    distance = np.linalg.norm(vector_subtraction)  # 求向量的模
    return distance


def calculate_polygon_area( vertices):
    """
    计算不规则多边形的面积，给定其顶点坐标。

    :param vertices: 顶点坐标列表，每个顶点是一个包含两个元素的列表[x, y]。
    :return: 多边形的面积。
    """
    n = len(vertices)  # 顶点的数量
    area = 0.0

    # 应用鞋带公式（Shoelace formula）计算多边形面积
    for i in range(n):
        j = (i + 1) % n  # 下一个顶点的索引（循环到最后一个顶点后回到第一个顶点）
        area += vertices[i][0] * vertices[j][1]  # 计算x1*y2
        area -= vertices[j][0] * vertices[i][1]  # 计算x2*y1

    # 返回面积的绝对值的一半
    return abs(area) / 2.0



# 创建 ECEF 到地理坐标的转换器
ecef = pyproj.Proj(proj='geocent', ellps='sphere')
lla = pyproj.Proj(proj='latlong', ellps='sphere')
# 使用 Transformer 类进行坐标转换
transformer_lla_to_ecef = pyproj.Transformer.from_proj(lla, ecef)
transformer_ecef_to_lla = pyproj.Transformer.from_proj(ecef, lla)


def ecef_to_lat_lon_alt(x, y, z):
    # 进行坐标转换
    lon, lat, alt = transformer_ecef_to_lla.transform(x, y, z)  # 输出顺序是，经度，纬度，海拔
    return [lat, lon, alt]  # 纬度，经度，海拔


# 纬度经度海拔 转 ECEF
def lat_lon_alt_to_ecef(lat, lon, alt):
    # 进行坐标转换
    x, y, z = transformer_lla_to_ecef.transform(lon, lat, alt)
    return [x, y, z]


def ecef_to_enu(x, y, z, reference_lat_lon_alt):
    lat0, lon0, alt0 = reference_lat_lon_alt
    # 将参考点的经纬度转换为 ECEF 坐标
    x0, y0, z0 = transformer_lla_to_ecef.transform(lon0, lat0, alt0)
    # 计算 ECEF 向量差
    dx = x - x0
    dy = y - y0
    dz = z - z0
    # 构建旋转矩阵
    lon_rad = np.radians(lon0)
    lat_rad = np.radians(lat0)
    R = np.array([[-np.sin(lon_rad), np.cos(lon_rad), 0],
                [-np.sin(lat_rad) * np.cos(lon_rad), -np.sin(lat_rad) * np.sin(lon_rad), np.cos(lat_rad)],
                [np.cos(lat_rad) * np.cos(lon_rad), np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)]])
    # 执行矩阵乘法
    enu = np.dot(R, np.array([dx, dy, dz]))
    E = enu[0]
    N = enu[1]
    U = enu[2]
    return [E, N, U]


def enu_to_ecef(e, n, u, reference_lat_lon_alt):
    lat0, lon0, alt0 = reference_lat_lon_alt
    # 将参考点的经纬度转换为 ECEF 坐标
    x0, y0, z0 = transformer_lla_to_ecef.transform(lon0, lat0, alt0)
    # 构建旋转矩阵
    lon_rad = np.radians(lon0)
    lat_rad = np.radians(lat0)
    R = np.array([[-np.sin(lon_rad), -np.sin(lat_rad) * np.cos(lon_rad), np.cos(lat_rad) * np.cos(lon_rad)],
                [np.cos(lon_rad), -np.sin(lat_rad) * np.sin(lon_rad), np.cos(lat_rad) * np.sin(lon_rad)],
                [0, np.cos(lat_rad), np.sin(lat_rad)]])
    # 执行矩阵乘法
    dx_dy_dz = np.dot(R, np.array([e, n, u]))
    # 计算 ECEF 坐标
    x = x0 + dx_dy_dz[0]
    y = y0 + dx_dy_dz[1]
    z = z0 + dx_dy_dz[2]
    return [x, y, z]


# 判断点是否在不规则多边形内部
def is_point_in_polygon( point, polygon):
    """
    从待判断的点向右画一条无限长的射线，然后计算这条射线与多边形边的交点数量。
    如果交点数量是奇数，则点在多边形内部；
    如果是偶数，则点在多边形外部。
    这个方法适用于凸多边形和凹多边形。
    :param point: 待判断点
    :param polygon: 多边形顶点
    :return: 是否在多边形内
    """
    if isinstance(polygon, np.ndarray):
        polygon = polygon.tolist()
    if len(polygon[0] ) >2:
        polygon = [[point[0], point[1]] for point in polygon]

    x, y = point
    n = len(polygon)
    inside = False
    for i in range(n):
        j = (i + 1) % n
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
            inside = not inside
    return inside # bool