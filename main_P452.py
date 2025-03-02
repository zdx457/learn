# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 19:10:52 2024

@author: Jason.Mao
"""

import numpy as np  
import math
from Py452 import P452
from global_land_mask import globe
import geopandas as gpd # 一个用于处理地理空间数据的库。它扩展了pandas库的数据处理功能，使其能够处理地理数据，如矢量数据（点、线、多边形等）
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString
"""
ITU-R P.452建议书
评估在频率高于约100 MHz时地球表面上电台之间干扰的预测程序
所需输入参数：
f 频率
p 确定一个时间比例，在这个时间比例内，计算得出的基本传输损耗不会被超过。本建议书要求，0.001<p50%
d 地形剖面距离，是一个数组，di表示到第i个剖面点发射机的距离
h 地形剖面高度.是一个数组，hi表示到第i个剖面点在海平面以上的地形高度
g hi+第i个剖面点的代表性地物高度
pol 极化方式： 1 - horizontal, 2 - vertical
zone  气候区类型，1 - 沿海陆地, 2 - 内陆, 3 - 海洋
temp 空气湿度
press 干燥空气压力
Gt\Gr 发射/接收天线增益
phit_j\phit_w 发射机的经度\纬度 _j是经度_w是纬度
phir_j\phir_w 接收机的经度\纬度 _j是经度_w是纬度
htg\hrg 发射机\接收机的天线高度，可以设置为10m，一般就是几米
dct\dcr 发射机\接收机到海岸的陆地距离
"""
# 双三次插值的权重获取
def bicubic(x):
    a = -0.5
    if 0 <= abs(x) <= 1:
        return (a + 2) * abs(x) ** 3 - (a + 3) * abs(x) ** 2 + 1
    elif 1 < abs(x) <= 2:
        return a * abs(x) ** 3 - 5 * a * abs(x) ** 2 + 8 * a * abs(x) - 4 * a
    else:
        return 0
    
# 计算从发射机到接收机的每个剖点的经纬度  
def calculate_intermediate_points(phit_j, phit_w, phir_j, phir_w, num_points):  
    points = []  
    central_angle = math.acos(math.sin(phit_w) * math.sin(phir_w) + math.cos(phit_w)*math.cos(phir_w)*math.cos(phit_j - phir_j))# 计算该路径在地球中心所对应的角
    # 计算每个剖点的经纬度  
    for i in range(num_points + 1):  
        f = i / num_points  # 当前插值点在路径中的位置  
        A = math.sin((1 - f) * (central_angle)) / math.sin(central_angle)  
        B = math.sin(f * (central_angle)) / math.sin(central_angle)  
        x = A * math.cos(phit_w) * math.cos(phit_j) + B * math.cos(phir_w) * math.cos(phir_j)  
        y = A * math.cos(phit_w) * math.sin(phit_j) + B * math.cos(phir_w) * math.sin(phir_j)  
        z = A * math.sin(phit_w) + B * math.sin(phir_w)  
        # 计算经纬度  
        lat = math.atan2(z, math.sqrt(x**2 + y**2))  
        lon = math.atan2(y, x)  
        # 转换为度  
        lat = math.degrees(lat)  
        lon = math.degrees(lon)  
        points.append((lat, lon))  
    return points  

"计算某点距离海岸线的距离"
def com_coast_dis(lat,lon):
    coastline = gpd.read_file(r'coastline/ne_50m_coastline.shp')#读取海岸线文件
    coastline = coastline.to_crs(epsg=3857)
    gdf = gpd.GeoDataFrame(geometry=coastline['geometry'])
    point = Point(lon, lat)
    # 计算点到每条海岸线的距离
    distances = gdf.distance(point)
    return distances/1000

"""
传入参数依次为:频率、概率，p:0.001<p50%,一般设置为50，极化方式，
"""
def com_452(f,p,pol,phit_j, phit_w, phir_j, phir_w,Gt, Gr,htg, hrg):
    "计算路径剖面，得到d"
    phit_j = math.radians(phit_j)#将输入的角度转化为弧度
    phit_w = math.radians(phit_w)
    phir_j = math.radians(phir_j)
    phir_w = math.radians(phir_w)
    central_angle = math.acos(math.sin(phit_w) * math.sin(phir_w) + math.cos(phit_w)*math.cos(phir_w)*math.cos(phit_j - phir_j))# 计算该路径在地球中心所对应的角
    d_tot = 6371 * central_angle #发射机和接收机之间的大圆距离，单位km
    # print(d_tot)
    if d_tot <= 6:
        d_0 = 0.03
        n = d_tot / d_0
        if n <=3:
            d_0 = 0.02
    elif d_tot >= 200:
        d_0 = 1
        n = d_tot / d_0
    elif 6 < d_tot < 200:
        d_0 = 0.03 + (d_tot - 6) * (1 - 0.03) / (200 - 6)
        n = d_tot / d_0
    d = [0]
    for i in range(1,math.ceil(n)+1):
        d.append(i*d_0)
    "计算每个剖面点的经纬度"
    num_points = math.ceil(n)  
    intermediate_points = calculate_intermediate_points(phit_j, phit_w, phir_j, phir_w, num_points)  
    "根据每个剖点的经纬度得到高于海平面高度"
    H = np.loadtxt('h_above sea level/TOPO.dat')#读取高度文件
    h = [0] * len(d)
    num_rows, num_cols = H.shape   
    # 生成经纬度  
    latitudes = np.linspace((90.125), (-90.125), num_rows)        # 从 +90 到 -90  
    longitudes = np.linspace((-180.125), (180.125), num_cols)     # 从 -180 到 +180  
    # 在纬度和经度中找到最接近的点  
    for index, point in enumerate(intermediate_points):
        lat, lon = point
        # 找到包围lat的两个纬度
        index_lat = 0
        index_lon = 0
        for index_i,i in enumerate(latitudes):
            if (latitudes[index_i] > lat) and (latitudes[index_i+1] < lat):
                index_lat = index_i
                break
        for index_j,j in enumerate(longitudes):
            if (longitudes[index_j] < lon) and (longitudes[index_j+1] >lon):
                index_lon = index_j
                break
        lat1 = latitudes[index_lat]
        lat2 = latitudes[index_lat+1]
        # 找到包围lon的两个经度
        lon1 = longitudes[index_lon]
        lon2 = longitudes[index_lon+1]
        # 计算相对偏移量
        u = (lat - lat1) / (lat2 - lat1)
        v = (lon - lon1) / (lon2 - lon1)
        # 双三次插值计算
        for ii in range(-1, 3):
            for jj in range(-1, 3):
                if 0 <= index_lat + ii < num_rows and 0 <= index_lon + jj < num_cols:
                    h[index] += H[index_lat + ii, index_lon + jj] * bicubic(ii - u) * bicubic(jj - v)
    "判断每个点的气候区"
    zone = [0] * len(d)#初始化
    for index,points in enumerate(intermediate_points):#判断是内陆还是海洋
        is_land = globe.is_land(points[0], points[1])
        if is_land:
            zone[index] = 2
            coast_diatance = com_coast_dis(points[0], points[1])
            if coast_diatance.min() < 50 and h[index] <= 100:
                zone[index] = 1
        else:
            zone[index] = 3
    "发射机和接收机距离海岸线的最近距离"
    dct = com_coast_dis(phit_w, phit_j).min()
    dcr = com_coast_dis(phir_w, phir_j).min()

    "得到发射机的温空气温度和干燥气体压力"
    d = np.array(d)
    h = np.array(h)
    zone = np.array(zone)
    g=h
    press=1013
    temp=15
    
    Lb = P452.bt_loss(f, p, d, h, g, zone, htg, hrg, phit_j, phit_w, phir_j, phir_w, Gt, Gr, pol, dct, dcr, press, temp)
    return Lb
if __name__ =='__main__':
    #              频率，概率，极化方式，发送经度，发送维度，接收经度，接收纬度，发送天线增益，接收天线增益，发射天线高度，接收天线高度
    Lb = com_452(10,50,1,20, 30, 21, 31,20, 20,10, 10)
    # Lb = com_452(2,50,1,2.027573238131661, 2.027573238131661, 0, 0,-15.49724295777618, -5,30, 1.5)

    print(f"Lb = {Lb}")