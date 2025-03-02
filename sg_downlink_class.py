# 星地间 下行链路干扰
import math
import numpy as np
from scipy.special import erfc
from antenna.sate_1528_15 import *
from antenna.ue_2101 import *


class SSDownlink:
    def __init__(self, angel1, angel2, D, f, T, bandwidth, bandwidth2, data_rate, theta, phi, pt_list,
                 num_satellites, gr_max, angel_1, angel_2, distance_to_satellite_list, wavelength_list, overlap_list,
                 tx_power, loss, distance_to_satellite, overlap):
        # 单条干扰链路参数
        self.angel1 = angel1  # 卫星1-地面站2与卫星1-地面1(即天线最大增益方向)间离轴角
        self.angel2 = angel2  # 卫星1-地面站2与卫星2-地面2线间离轴角
        self.D = D  # 天线直径
        self.f = f  # 通信卫星的频率，单位为Hz
        self.k = 1.38e-23  # 玻尔兹曼常数，单位：J/K
        self.T = T  # 为接收机噪声温度  单位：K
        self.bandwidth = bandwidth  # 发射卫星系统的带宽，单位：赫兹
        self.bandwidth2 = bandwidth2  # bandwidth2为地面接收机的工作带宽
        self.data_rate = data_rate  # 信号传输速率,单位：比特每秒
        self.tx_antenna_gain = NGSO_SAT_1528_1(D, f, angel1)
        self.rx_antenna_gain = -5  # 手机终端接收增益 dBi
        self.theta = theta  # 基站用户仰角
        self.phi = phi  # 基站用户方向角
        # n条干扰链路参数
        self.pt_list = pt_list  # 干扰卫星的发射功率 瓦特
        self.num_satellites = num_satellites  # 同频工作卫星个数
        self.gr_max = gr_max  # 地面站天线的峰值增益 （单位：dB）
        self.angel_1 = angel_1  # 卫星i-地面站2与卫星i-地面i线间离轴角
        self.angel_2 = angel_2  # 卫星i-地面站2与卫星2-地面2线间离轴角
        self.distance_to_satellite_list = distance_to_satellite_list  # 每颗卫星到接收机的距离列表
        self.wavelength_list = wavelength_list  # 每颗卫星的波长列表
        self.overlap_list = overlap_list  # 每颗干扰卫星与接收频段的重叠带宽列表
        # 单个干扰卫星参数
        self.tx_power = tx_power  # 通信链路 发射端天线功率，单位：瓦特
        self.loss = loss  # 其他损耗，单位：dBi
        self.distance_to_satellite = distance_to_satellite  # 卫星与地面站之间的距离（单位：千米）
        self.overlap = overlap  # 干扰卫星与接收频段的重叠带宽
        self.acs = 5  # dB
        self.aclr = 25  # dB

    def calculate_received_power(self):
        loss_pg = 32.45 + 20 * math.log10(self.distance_to_satellite) + 20 * math.log10(self.f / 1000000)
        received_power = 10 * math.log10(
            self.tx_power) + self.tx_antenna_gain + self.rx_antenna_gain - loss_pg - self.loss
        return received_power

    def calculate_eb_n0(self):
        k = 1.38e-23  # 玻尔兹曼常数，单位：J/K
        noise_power_density = k * self.T * self.bandwidth
        eb_n0 = 10 ** ((self.calculate_received_power() - 10 * math.log10(self.data_rate) - 10 * math.log10(
            noise_power_density)) / 10)
        return eb_n0

    def calculate_ber_2ask(self):
        return 0.5 * erfc(np.sqrt(self.calculate_eb_n0()))

    # 计算频谱效率
    def calculate_spectral_efficiency(self):
        spectral_efficiency = self.data_rate / self.bandwidth
        return spectral_efficiency

    # 计算BER for BPSK
    def calculate_ber_bpsk(self):
        ber_bpsk = 0.5 * erfc(np.sqrt(self.calculate_eb_n0()))
        return ber_bpsk

    # 计算BER for QPSK
    def calculate_ber_qpsk(self):
        ber_qpsk = 0.5 * erfc(np.sqrt(self.calculate_eb_n0()))
        return ber_qpsk

    # 计算BER for 16-QAM
    def calculate_ber_16qam(self):
        M = 16
        ber_16qam = (3 / 2) * (1 - (1 / np.sqrt(M))) * erfc(
            np.sqrt((3 * np.log2(M) / (M - 1)) * self.calculate_eb_n0()))
        return ber_16qam

    def calculate_N(self):
        N = 10 * math.log10(self.k * self.T * self.bandwidth2)
        return N

    def calculate_I_N(self):
        num_terms = self.num_satellites
        term1 = 0
        fsl = 32.45 + 20 * math.log10(self.distance_to_satellite) + 20 * math.log10(self.f / 1000000)

        for i in range(num_terms):
            gt_list = NGSO_SAT_1528_1(self.D, self.f, self.angel_1[i])
            rx_antenna_gain = self.rx_antenna_gain  # 设置固定值
            gt_list = 10 ** (gt_list / 10)  # 转换成线性值
            rx_antenna_gain = 10 ** (rx_antenna_gain / 10)  # 转换成线性值
            tr_power = self.pt_list[i] * self.overlap_list[i] / self.bandwidth  # 干扰功率
            wavelength = self.wavelength_list[i]  # 下行链路工作载波的波长
            acir = 1 / (1 / self.aclr + 1 / self.acs)  # dB aclr是发送端卫星-邻道泄漏比  acs是接收端手机-邻道选择性
            acir = 10 ** (acir / 10)  # 线性
            term = tr_power * gt_list * rx_antenna_gain * (
                    wavelength**2 / (acir * (4 * math.pi * 1000 * self.distance_to_satellite_list[i])) ** 2)
            term1 += term

        I = 10 * math.log10(term1) # 集总
        I = I - fsl
        I_N = I - 10 * math.log10(self.k * self.T * self.bandwidth2)

        return I_N, I

    def calculate_C_N(self):
        distance_to_satellite = 20 * math.log10(self.distance_to_satellite * 1000)
        loss_pg = 32.45 + 20 * math.log10(self.distance_to_satellite) + 20 * math.log10(self.f / 1000000)
        tx_power_db = 10 * math.log10(self.tx_power)
        tx_antenna_gain = IMT_2101(self.theta, self.phi)  # 仰角 方位角
        rx_antenna_gain = self.rx_antenna_gain  # 设置固定值
        wavelength_db = 20 * math.log10(self.wavelength_list[0])  # 下行链路工作载波的波长
        C = tx_power_db + tx_antenna_gain + rx_antenna_gain + wavelength_db - distance_to_satellite - loss_pg   # 基站发射载波功率
        N = 10 * math.log10(self.k * self.T * self.bandwidth2)
        C_N = C - N
        return C_N, C

    def calculate_epfd(self):
        num_terms = self.num_satellites
        epfd = 0
        term1 = 0
        gt_list = [0] * num_terms

        for i in range(num_terms):
            gt_list = NGSO_SAT_1528_1(self.D, self.f, self.angel_1[i])
            rx_antenna_gain = self.rx_antenna_gain

            gt_list = 10 ** (gt_list / 10)
            rx_antenna_gain = 10 ** (rx_antenna_gain / 10)
            term = (self.pt_list[i] * gt_list * rx_antenna_gain) / (self.overlap_list[i] *
                                                                    1000 * self.distance_to_satellite_list[
                                                                        i] ** 2 * self.gr_max * 4 * math.pi)
            term1 += term
            epfd = 10 * math.log10(term1)

        return epfd

    def calculate_pfd(self):
        pfd = 0
        num_terms = self.num_satellites

        for i in range(num_terms):
            gt_list = NGSO_SAT_1528_1(self.D, self.f, self.angel_1[i])
            pfd = 10 * math.log10(self.pt_list[i]) + gt_list - 10 * math.log10(
                4 * math.pi * self.distance_to_satellite_list[i] ** 2) + 10 * math.log10(self.bandwidth)

        return pfd

    def calculate_t(self):
        t = 0
        num_terms = self.num_satellites


        for i in range(num_terms):
            gt_list = NGSO_SAT_1528_1(self.D, self.f, self.angel_1[i])
            rx_antenna_gain = self.rx_antenna_gain
            # gt_list = gt_list_[0]
            gt_list = 10 ** (gt_list / 10)
            rx_antenna_gain = 10 ** (rx_antenna_gain / 10)
            term = (self.pt_list[i] * gt_list * rx_antenna_gain) * (3e8 / self.f) ** 2 / (
                    self.k * (4 * math.pi * 1000 * self.distance_to_satellite_list[i]) ** 2 * self.T)
            t += term
            t = 10 * math.log10(t)

        return t

    # 定义第i颗卫星特定用户的容量channel_capacity
    # distance_to_satellite_list 为卫星与用户之间的距离
    # bandwidth为LEO卫星带宽
    def calculate_channel_capacity(self, I):
        loss_sh = 5
        loss_atm = 5
        gt = 10 ** (self.tx_antenna_gain / 10)
        gr = 10 ** (self.rx_antenna_gain / 10)
        loss_sh = 10 ** (loss_sh / 10)
        loss_atm = 10 ** (loss_atm / 10)
        wavelength = 3e8 / self.f
        h = gt * gr / (((4 * math.pi * 1000 * self.distance_to_satellite) / wavelength) ** 2 * loss_sh * loss_atm)
        I = 10 ** (I / 10)
        sinr = self.tx_power * h / (I + self.k * self.T * self.bandwidth2)
        channel_capacity = self.bandwidth * np.log2(1 + sinr)
        return channel_capacity

    def calculate_C_I(self, C, I):
        C_I = C - I
        return C_I

    def calculate_C_NI(self, C, I):
        N = self.k * self.T * self.bandwidth2
        C = 10 ** (C / 10)
        I = 10 ** (I / 10)
        C_NI = C / (N + I)
        C_NI = 10 ** (C_NI / 10)
        return C_NI

    def calculate_i_N(self, i):
        i_N = i - 10 * math.log10(self.k * self.T * self.bandwidth2)
        return i_N

    # 计算单链路干扰 i 单位：dB
    def calculate_i(self):
        tr_power = self.pt_list[0] * self.overlap / self.bandwidth
        tr_power_db = 10 * math.log10(tr_power)
        wavelength = 3e8 / self.f
        i = tr_power_db + self.tx_antenna_gain + self.rx_antenna_gain + 20 * math.log10(
            wavelength / (1000 * self.distance_to_satellite))
        return i


# 实例化类
ss_downlink = SSDownlink(15, 5, 2.7, 5e9, 300, 5e6, 5e6, 3e6, [1000, 900, 800], 3, 10, [10, 15, 10],
                         [15, 10, 10], [400, 500, 600], [0.06, 0.05, 0.06], [5000, 4000, 4500], 20.0, 0, 600, 4000)

# 调用类方法
received_power = ss_downlink.calculate_received_power()
eb_n0 = ss_downlink.calculate_eb_n0()
ber_2ask = ss_downlink.calculate_ber_2ask()
spectral_efficiency = ss_downlink.calculate_spectral_efficiency()
ber_bpsk = ss_downlink.calculate_ber_bpsk()
ber_qpsk = ss_downlink.calculate_ber_qpsk()
ber_16qam = ss_downlink.calculate_ber_16qam()
I_N, I = ss_downlink.calculate_I_N()
i = ss_downlink.calculate_i()
i_N = ss_downlink.calculate_i_N(i)
C_N, C = ss_downlink.calculate_C_N()
C_I = ss_downlink.calculate_C_I(C, I)
C_i = ss_downlink.calculate_C_I(C, i)
N = ss_downlink.calculate_N()
C_NI = ss_downlink.calculate_C_NI(C, I)
C_Ni = ss_downlink.calculate_C_NI(C, i)
epfd = ss_downlink.calculate_epfd()
pfd = ss_downlink.calculate_pfd()
t = ss_downlink.calculate_t()
channel_capacity = ss_downlink.calculate_channel_capacity(I)

print("接收端信号功率为:", received_power, "dB")
print("Eb/N0 值为:", eb_n0)
print("BER for 2ASK:", ber_2ask)
print(f"BER for BPSK: {ber_bpsk:.2f}")
print(f"BER for QPSK: {ber_qpsk:.2f}")
print(f"BER for 16-QAM: {ber_16qam:.2f}")

print("I/N (集总干扰与噪声比): {:.2f} dBW".format(I_N))
print("i/N (干扰与噪声比): {:.2f} dB".format(i_N))
print("I:", I, "dBW")
print("i:", i, "dBW")
print("N:", N, "dB")
print("C:", C, "dBW")
print("C/N (载波与噪声比): {:.2f} dB".format(C_N))
print("C/I (载波与集总干扰比): {:.2f} dBW".format(C_I))
print("C/i (载波与干扰比): {:.2f} dB".format(C_i))
print("C/(N+I) (载波与噪声+集总干扰比): {:.2f} dBW".format(C_NI))
print("C/(N+i) (载波与噪声+干扰比): {:.2f} dB".format(C_Ni))

print("EPFD (等效载噪比密度): {:.2f} ".format(epfd))
print("PFD : {:.2f} dBW/m2".format(pfd))
print("delta_T/T: {:.2f} ".format(t))
print("频谱效率: {:.2f} ".format(spectral_efficiency))
print("信道容量: {:.2f} bps".format(channel_capacity))
