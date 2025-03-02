# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:22:57 2024

@author: Jason.Mao
"""
import math
import numpy as np 
import matplotlib.pyplot as plt  

"""
NGSO天线
"""
"ITU-R.S 1528" 
#L_N=-15
def NGSO_SAT_1528_1(D,f,phi): 
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = 20*math.log10(D/lamda) + 7.7
    phi_b = math.sqrt(1200)/(D_lamda)    #波束半宽度 单位角度
    a = 2.58
    b = 6.32
    L_F = 0
    G = None  # 初始化变量 G  
    if D_lamda >= 35:
        L_N = -15
        L_B = max(15 + L_N + 0.25 * G_m, 0)  
        X = G_m + L_N + 25 * math.log10(b * phi_b)
        Y = b * phi_b * 10 ** (0.04 * (G_m + L_N - L_F))
        if 0 <= phi <= a * phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif a * phi_b < phi <= b * phi_b:  
            G = G_m + L_N
        elif b * phi_b < phi <= Y:  
            G = X - 25 * math.log10(phi) 
        elif Y < phi <= 90:  
            G = L_F  
        elif 90 < phi <= 180:  
            G = L_B
        else:  
            G = None  # 超出范围  
    elif D_lamda < 35:
        L_s = -6.75
        Y = 1.5 * phi_b
        Z = Y * 10 ** (0.04 * (G_m + L_s - L_F))  
        if 0 <= phi <= phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif phi_b < phi <= Y:  
            G = G_m - 3 * (phi / phi_b) ** 2  
        elif Y < phi <= Z:  
            G = G_m + L_s -25 * math.log10(phi/Y)
        elif Z < phi <= 180:  
            G = L_F
    return G,G_m
#L_N=-20    
def NGSO_SAT_1528_2(D,f,phi):
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = 20*math.log10(D/lamda) + 7.7
    phi_b = math.sqrt(1200)/(D_lamda)
    a = 2.58
    b = 6.32
    L_F = 0
    G = None  # 初始化变量 G  
    if D_lamda >= 35:
        L_N = -20
        L_B = max(15 + L_N + 0.25 * G_m, 0)  
        X = G_m + L_N + 25 * math.log10(b * phi_b)
        Y = b * phi_b * 10 ** (0.04 * (G_m + L_N - L_F))
        if 0 <= phi <= a * phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif a * phi_b < phi <= b * phi_b:  
            G = G_m + L_N
        elif b * phi_b < phi <= Y:  
            G = X - 25 * math.log10(phi) 
        elif Y < phi <= 90:  
            G = L_F  
        elif 90 < phi <= 180:  
            G = L_B
        else:  
            G = None  # 超出范围  
    elif D_lamda < 35:
        L_s = -6.75
        Y = 1.5 * phi_b
        Z = Y * 10 ** (0.04 * (G_m + L_s - L_F))  
        if 0 <= phi <= phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif phi_b < phi <= Y:  
            G = G_m - 3 * (phi / phi_b) ** 2  
        elif Y < phi <= Z:  
            G = G_m + L_s -25 * math.log10(phi/Y)
        elif Z < phi <= 180:  
            G = L_F
    return G,G_m    
#L_N=-25
def NGSO_SAT_1528_3(D,f,phi): 
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = 20*math.log10(D/lamda) + 7.7
    phi_b = math.sqrt(1200)/(D_lamda)
    a = 2.58
    b = 6.32
    L_F = 0
    G = None  # 初始化变量 G  
    if D_lamda >= 35:
        L_N = -25
        L_B = max(15 + L_N + 0.25 * G_m, 0)  
        X = G_m + L_N + 25 * math.log10(b * phi_b)
        Y = b * phi_b * 10 ** (0.04 * (G_m + L_N - L_F))
        if 0 <= phi <= a * phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif a * phi_b < phi <= b * phi_b:  
            G = G_m + L_N
        elif b * phi_b < phi <= Y:  
            G = X - 25 * math.log10(phi) 
        elif Y < phi <= 90:  
            G = L_F  
        elif 90 < phi <= 180:  
            G = L_B
        else:  
            G = None  # 超出范围  
    elif D_lamda < 35:
        L_s = -6.75
        Y = 1.5 * phi_b
        Z = Y * 10 ** (0.04 * (G_m + L_s - L_F))  
        if 0 <= phi <= phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif phi_b < phi <= Y:  
            G = G_m - 3 * (phi / phi_b) ** 2  
        elif Y < phi <= Z:  
            G = G_m + L_s -25 * math.log10(phi/Y)
        elif Z < phi <= 180:  
            G = L_F
    return G,G_m
#L_N=-30
def NGSO_SAT_1528_4(D,f,phi): 
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = 20*math.log10(D/lamda) + 7.7
    phi_b = math.sqrt(1200)/(D_lamda)
    a = 2.58
    b = 6.32
    L_F = 0
    G = None  # 初始化变量 G  
    if D_lamda >= 35:
        L_N = -30
        L_B = max(15 + L_N + 0.25 * G_m, 0)  
        X = G_m + L_N + 25 * math.log10(b * phi_b)
        Y = b * phi_b * 10 ** (0.04 * (G_m + L_N - L_F))
        if 0 <= phi <= a * phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif a * phi_b < phi <= b * phi_b:  
            G = G_m + L_N
        elif b * phi_b < phi <= Y:  
            G = X - 25 * math.log10(phi) 
        elif Y < phi <= 90:  
            G = L_F  
        elif 90 < phi <= 180:  
            G = L_B
        else:  
            G = None  # 超出范围  
    elif D_lamda < 35:
        L_s = -6.75
        Y = 1.5 * phi_b
        Z = Y * 10 ** (0.04 * (G_m + L_s - L_F))  
        if 0 <= phi <= phi_b:  
            G = G_m - 3 * (phi / phi_b) ** 1.5  
        elif phi_b < phi <= Y:  
            G = G_m - 3 * (phi / phi_b) ** 2  
        elif Y < phi <= Z:  
            G = G_m + L_s -25 * math.log10(phi/Y)
        elif Z < phi <= 180:  
            G = L_F
    return G,G_m
"AP30" 



""""
地球站天线
"""


"ITU-R S.465-6 NGSO地球站天线"
#发送端
def ES_465_6_tr(D,f,phi):
    eta = 0.7 #天线效率
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = math.log10(eta * (math.pi**2) * (D/lamda)**2) * 10
    phi_1 = 0.9 * 114* (D_lamda)**(-1.09)
    phi_r = 15.85 * (D_lamda)**(-0.6)
    G_1 = 32 - 25 * math.log10(phi_r)
    phi_m = 20 * (1 / D_lamda) * math.sqrt(G_m - G_1) 
    if D_lamda >= 50:
        phi_min = max(1 , 100 * (1/ D_lamda))
    elif  D_lamda < 50:
        phi_min = max(2 , 114 * D_lamda ** (-1.09))
    phi_b = 10**(42/25)
    G = None
    if D_lamda <= 54.5:
        if 0 <= phi <phi_1:
            G = G_m - 2.5 * 10**-3 * (D_lamda * phi)**2  
        elif phi_1 <= phi < phi_min:
            G = max(G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2, 32 - 25 * math.log10(phi))  
        elif phi_min <= phi < 180:
            G = max(32 - 25 * math.log10(phi), -10) 
    elif D_lamda > 54.5:
        if 0 <= phi <phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2  
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi < 180:
            G = max(32 - 25 * math.log10(phi), -10) 
    return G

#接收端
def ES_465_6_re(D,f,phi):
    eta = 0.7 #天线效率
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = math.log10(eta * (math.pi**2) * (D_lamda)**2) * 10
    phi_1 = 0.9 * 114* (D_lamda)**(-1.09)
    phi_r = 15.85 * (D_lamda)**(-0.6)
    G_1 = 32 - 25 * math.log10(phi_r)
    phi_m = 20 * (1 / D_lamda) * math.sqrt(G_m - G_1) 
    if D_lamda >= 50:
        phi_min = max(1 , 100 * (1/ D_lamda))
    elif  D_lamda < 50:
        phi_min = max(2 , 114 * D_lamda ** (-1.09))
        if phi_min > 2.5:
            phi_min = 2.5
    phi_b = 10**(42/25)
    G = None
    if D_lamda < 33.3:
        if 0 <= phi < phi_min:
            G = G_m - 2.5 * 10**-3 * (D_lamda * phi)**2  
        elif phi_min < phi <= 180:
            G = max(32 - 25 * math.log10(phi), -10) 
    elif 33.3 <= D_lamda <= 54.5:
        if 0 <= phi <phi_1:
            G = G_m - 2.5 * 10**-3 * (D_lamda * phi)**2  
        elif phi_1 <= phi < phi_min:
            G = max(G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2, 32 - 25 * math.log10(phi))  
        elif phi_min <= phi < 180:
            G = max(32 - 25 * math.log10(phi), -10) 
    elif D_lamda > 54.5:
        if 0 <= phi <phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2  
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi < 180:
            G = max(32 - 25 * math.log10(phi), -10) 
    return G

"ITU-R S.465-5 NGSO地球站天线"    
def ES_465_5(D, f, phi):
   eta = 0.7 #天线效率
   phi_b = 10 ** (42/25)
   lamda = 3 * 10**8 / f
   D_lamda = D/lamda
   G_m = math.log10(eta * (math.pi**2) * (D/lamda)**2) * 10
   if D_lamda > 100:
       G_1 = 32
   else:
       G_1 = -18 + 25 * math.log10(D_lamda)
   phi_m = 20 / D_lamda * math.sqrt(G_m - G_1)
   if D_lamda > 100:
       phi_r = 1
   else :
       phi_r = 100 / D_lamda
   if 0 <= phi < phi_m:
       G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
   elif phi_m <= phi < phi_r:
       G = G_1
   elif phi_r <= phi < phi_b:
       G = 32 -25 * math.log10(phi)
   elif phi_b <= phi <= 180:
       G = -10
   return G

"ITU-R S.1428 NGSO地球站天线"
#用于在 10.7 GHz 和 30 GHz 之间的频段内涉及非 GSO 卫星的平共处干扰评估的参考 FSS 地球站的辐射方向图
def ES_1428(D, f, phi):
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G = None
    if D_lamda < 20:
        G = None
    elif 20 <= D_lamda <=25:
        G_m = 20 * math.log10(D_lamda) + 7.7
        G_1 = 29 - 25 * math.log10(95 / D_lamda)
        phi_m = 25 / D_lamda * math.sqrt(G_m - G_1)
        if 0 < phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < 95/D_lamda:
            G = G_1
        elif 95/D_lamda <= phi < 33.1:
            G = 29 - 25 * math.log10(phi)
        elif 33.1 <= phi < 80:
            G = -9
        elif 80 <= phi <= 180:
            G = -5
    elif 25 < D_lamda <= 100:
        G_m = 20 * math.log10(D_lamda) + 7.7
        G_1 = 29 - 25 * math.log10(95 / D_lamda)
        phi_m = 25 / D_lamda * math.sqrt(G_m - G_1)
        if 0 < phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < 95/D_lamda:
            G = G_1
        elif 95/D_lamda <= phi < 33.1:
            G = 29 - 25 * math.log10(phi)
        elif 33.1 <= phi < 80:
            G = -9
        elif 80 <= phi <= 120:
            G = -4       
        elif 120 < phi <= 180:
            G = -9
    elif 100 < D_lamda:
        G_m = 20 * math.log10(D_lamda) + 8.4
        G_1 = -1 + 15 * math.log10(D_lamda)
        phi_m = 20 / D_lamda * math.sqrt(G_m - G_1)
        phi_r = 15.85 * D_lamda**(-0.6)
        if 0 < phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi < 10:
            G = 29 - 25 * math.log10(phi)
        elif 10 <= phi < 34.1:
            G = 34 - 30 * math.log10(phi)
        elif 34.1 <= phi < 80:
            G = -12
        elif 80 <= phi < 120:
            G = -7
        elif 120 <= phi <= 180:
            G = -12
    return G

"""
ITU AP30
本附录的规定适用于第3区11.7 GHz至12.2 GHz频段、第1区11.7 GHz至12.5 GHz频段、第2区12.2 GHz至12.7 GHz频段的卫星广播业务
以及第1、2和3区分配给这些频段的其他业务，只要它们与这些频段的卫星广播业务的关系有关。
"""           
# AP30 539页天线方向图 565接收天线  571发射天线，704接收天线，713发射天线

"ITU AP30_B 地球站天线 "   
def ES_AP_30_B(D,f,phi):
    # 使用ITU-R S.580-6对其进行拓展，5.5 m for the 6/4 GHz band; 2.7 m for the 13/10-11 GHz band. (WRC-07) 
    lamda = 3 * 10**8 / f
    D_lamda = D / lamda
    eta = 0.7
    G_m = 10 * math.log10(eta * (math.pi * D_lamda)**2)
    phi_b = 10**(42/25)
    if D_lamda >= 100:
        phi_r = 15.85 * (D_lamda)**(-0.6)
        G_1 = -1 + 15 * math.log10(D_lamda)
    elif D_lamda < 100:
        phi_r = 100/D_lamda
        if D_lamda < 50:
            G_1 = 2 + 15 * math.log10(D_lamda)
        elif 50 <= D_lamda < 100:
            G_1 = -21 + 25 * math.log10(D_lamda)
    phi_m = 20 / D_lamda * math.sqrt(G_m - G_1)
    if  D_lamda >= 50:
        if 0 <= phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi <= 19.95:
            G = 29 - 25 * math.log10(phi)
        elif 19.95 < phi < phi_b:
            G = min(-3.5 , 32 - 25 * math.log10(phi))
        elif phi_b <= phi <=180:
            G = -10
    elif D_lamda < 50:
        if 0 <= phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi < phi_b:
            G = 52 - 10 * math.log10(D_lamda) - 25 * math.log10(phi)
        elif phi_b <= phi <=180:
            G = 10 - 10 * math.log10(D_lamda) 
    return G  
"ITU AP30 "
                
"ITU AP8 地球站天线"        
def ES_AP_8(D,f,phi):
    lamda = 3 * 10**8 / f
    D_lamda = D/lamda
    G_m = 20 * math.log10(D_lamda) + 7.7
    G_1 = 2 + 15 * math.log10(D_lamda)
    phi_b = 48
    phi_m = 20 / D_lamda * math.sqrt(G_m - G_1)
    if D_lamda >= 100:
        phi_r = 15.85 * D_lamda**(-0.6)
        if 0 <= phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi < phi_b:
            G = 32 - 25 * math.log10(phi)
        elif phi_b <= phi <=180:
            G = -10
    elif D_lamda < 100:
        phi_r = 100/D_lamda
        if 0 <= phi < phi_m:
            G = G_m - 2.5 * 10**(-3) * (D_lamda * phi)**2
        elif phi_m <= phi < phi_r:
            G = G_1
        elif phi_r <= phi < phi_b:
            G = 52 - 10 * math.log10(D_lamda) - 25 * math.log10(phi)
        elif phi_b <= phi <=180:
            G = 10 - 10 * math.log10(D_lamda)
    return G


"""
基站和用户设备天线
输入：天线的几何图形：仰角：theta 0——180°；方位角：phi -180——180°
"""
def IMT_2101(theta,phi):
    """
    单元方向图
    """
    # 单个单元的水平3dB带宽/度
    phi_3dB = 65
    # 单个单元的垂直3dB带宽/度
    theta_3dB = 65
    # 前后比/dB
    Am = 30
    SLAv = 30
    # 元件增益/dBi
    Gmax = 5
    # 水平辐射方向图
    A_eh_phi = -min(12 * (phi/phi_3dB)**2 , Am)
    # 垂直辐射方向图
    A_ey_theta = -min(12 * ((theta-90)/theta_3dB)**2 , SLAv)
    # 单个元件方向图
    A_E = Gmax - min(-(math.floor(A_eh_phi + A_ey_theta)), Am)
    """
    复合天线方向图
    """  
    # 天线阵列配置（行*列）
    NH = 8 # 水平
    NV = 8 # 垂直
    # 辐射单元间隔
    d_lamda_H = 0.5
    d_lamda_V = 0.5
    # 下倾角
    theta_etilt = 10
    phi_escan = 10
    # 叠加矢量计算
    v_nm = np.zeros((NV, NH), dtype=complex)
    w_nm = np.zeros((NV, NH), dtype=complex)
    for n in range(1,NV+1):
        for m in range(1,NH+1):          
            v_nm[n-1,m-1]  = np.exp( 1j * 2*math.pi * ( (n-1) * d_lamda_V * math.cos(math.radians(theta)) + (m-1) * d_lamda_H * math.sin(math.radians(theta)) * math.sin(math.radians(phi))))
    # 加权计算
    for n in range(1,NV+1):
        for m in range(1,NH+1):
            w_nm[n-1,m-1] = 1/math.sqrt(NH*NV) * np.exp(1j  * 2*math.pi * ((n-1) * d_lamda_V * math.sin(math.radians(theta_etilt)) - (m-1) * d_lamda_H * math.cos(math.radians(theta_etilt)) * math.sin(math.radians(phi_escan))))
    A = A_E +10 * np.log10(np.abs(np.sum(w_nm * v_nm))**2)
    return A         
