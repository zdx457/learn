import math
import numpy as np

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