import numpy as np
import plotly.graph_objects as go
# 定义函数：根据轨道要素计算卫星的ECI坐标
def orbital_elements_to_ECI(sat):
    # 提取轨道要素
    a = sat['a']
    e = sat['e']
    i = sat['i']
    Omega = sat['Omega']
    omega = sat['omega']
    f = sat['f']

    # 计算轨道半径
    r = a * (1 - e ** 2) / (1 + e * np.cos(f))

    # 计算在轨道平面坐标系下的位置
    x_op = r * np.cos(f)
    y_op = r * np.sin(f)
    z_op = 0.0  # 在轨道平面内，z坐标为0

    # 构建旋转矩阵
    R_omega = np.array([
        [np.cos(omega), -np.sin(omega), 0],
        [np.sin(omega), np.cos(omega), 0],
        [0, 0, 1]
    ])

    R_i = np.array([
        [1, 0, 0],
        [0, np.cos(i), -np.sin(i)],
        [0, np.sin(i), np.cos(i)]
    ])

    R_Omega = np.array([
        [np.cos(Omega), -np.sin(Omega), 0],
        [np.sin(Omega), np.cos(Omega), 0],
        [0, 0, 1]
    ])

    # 组合旋转矩阵
    R = np.dot(R_Omega, np.dot(R_i, R_omega))

    # 计算在ECI坐标系下的位置
    r_eci = np.dot(R, np.array([x_op, y_op, z_op]))

    return r_eci

# 定义函数：将ECI坐标转换为经纬度坐标
def ECI_to_LLH(r_eci):
    R_e = 6371.0
    x, y, z = r_eci
    # 计算经度
    lon = np.arctan2(y, x)
    # 计算地心纬度
    lat = np.arcsin(z / np.linalg.norm(r_eci))
    # 计算高度
    alt = np.linalg.norm(r_eci) - R_e
    # 转换为度
    lat_deg = np.degrees(lat)
    lon_deg = np.degrees(lon)
    return lat_deg, lon_deg, alt


