# 星地间 下行链路干扰
import math
import numpy as np
from scipy.special import erfc
from antenna.sate_1528_15 import *
from antenna.ue_2101 import *


class SGDownlink:
    def __init__(self, angel1, angel2, D, f, T, bandwidth, bandwidth2, data_rate, theta, phi,
                 tx_power, tu_power, loss, distance_to_satellite, overlap):
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
        self.tx_antenna_gain = NGSO_SAT_1528_1(D, f, angel1)  # 干扰发射卫星天线增益
        self.rx_antenna_gain = -5  # 手机终端接收增益 dBi
        self.theta = theta  # 基站用户仰角
        self.phi = phi  # 基站用户方向角
        # 单个干扰卫星参数
        self.tx_power = tx_power  # 干扰链路 发射端天线功率，单位：瓦特
        self.tu_power = tu_power  # 基站发射功率
        self.loss = loss  # 其他损耗，单位：dBi
        self.distance_to_satellite = distance_to_satellite  # 卫星与地面站之间的距离（单位：千米）
        self.overlap = overlap  # 干扰卫星与接收频段的重叠带宽

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

    def calculate_h(self):
        wavelength = 3e8 / self.f
        gt = NGSO_SAT_1528_1(self.D, self.f, self.angel1)
        gt = 10 ** (gt / 10)  # 转换成线性值
        rx_antenna_gain = 10 ** (self.rx_antenna_gain / 10)  # 转换成线性值
        h = gt * rx_antenna_gain * (
                (wavelength ** 2) / (4 * math.pi * 1000 * self.distance_to_satellite) ** 2)
        h = 10 * math.log10(h)
        return h

    # 计算单链路干扰 I 单位：dB
    def calculate_I(self, h):
        tr_power = self.tx_power * self.overlap / self.bandwidth
        tr_power_db = 10 * math.log10(tr_power)
        fsl = 32.45 + 20 * math.log10(self.distance_to_satellite) + 20 * math.log10(self.f / 1000000)
        I = tr_power_db + h - fsl #dB
        return I

    def calculate_I_N(self, I, N):
        I_N = I - N
        return I_N

    def calculate_C_N(self, N):
        distance_to_satellite = 20 * math.log10(self.distance_to_satellite * 1000)
        loss_pg = 32.45 + 20 * math.log10(self.distance_to_satellite) + 20 * math.log10(self.f / 1000000)  # 自由空间损耗
        tx_power_db = 10 * math.log10(self.tu_power)
        tx_antenna_gain = IMT_2101(self.theta, self.phi)  # 仰角 方位角
        rx_antenna_gain = self.rx_antenna_gain  # 设置固定值
        wavelength_db = 20 * math.log10(3e8 / self.f)  # 下行链路工作载波的波长
        C = tx_antenna_gain + rx_antenna_gain + tx_power_db + wavelength_db - distance_to_satellite - loss_pg  # 基站
        C_N = C - N
        return C_N, C

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

    def calculate_epfd(self):
        epfd = 0
        gt = NGSO_SAT_1528_1(self.D, self.f, self.angel1)  # 干扰卫星发射增益
        gt = 10 ** (gt / 10)
        term = (self.tx_power * gt) / (self.overlap *  # 接收增益按固定值算，等于最大增益，故公式里抵消
                                       1000 * self.distance_to_satellite ** 2 * 4 * math.pi)
        epfd = 10 * math.log10(term)

        return epfd

    def calculate_pfd(self):
        gt = NGSO_SAT_1528_1(self.D, self.f, self.angel1)  # 干扰卫星发射增益
        pfd = 10 * math.log10(self.tx_power) + gt - 10 * math.log10(
            4 * math.pi * self.distance_to_satellite ** 2) + 10 * math.log10(self.bandwidth)

        return pfd

    def calculate_t(self):
        gt = NGSO_SAT_1528_1(self.D, self.f, self.angel1)
        rx_antenna_gain = self.rx_antenna_gain
        gt = 10 ** (gt / 10)
        rx_antenna_gain = 10 ** (rx_antenna_gain / 10)
        term = (self.tx_power * gt * rx_antenna_gain) * (3e8 / self.f) ** 2 / (
                self.k * (4 * math.pi * 1000 * self.distance_to_satellite) ** 2 * self.T)
        t = 10 * math.log10(term)

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


if __name__ == '__main__':
    # 实例化类
    ss_downlink = SGDownlink(15, 5, 2.7, 5e9, 300, 5e6, 5e6, 3e6, 10, 5, 10,
                             20, 0, 600, 4000)

    # 调用类方法
    received_power = ss_downlink.calculate_received_power()
    eb_n0 = ss_downlink.calculate_eb_n0()
    ber_2ask = ss_downlink.calculate_ber_2ask()
    spectral_efficiency = ss_downlink.calculate_spectral_efficiency()
    ber_bpsk = ss_downlink.calculate_ber_bpsk()
    ber_qpsk = ss_downlink.calculate_ber_qpsk()
    ber_16qam = ss_downlink.calculate_ber_16qam()
    h = ss_downlink.calculate_h()
    I = ss_downlink.calculate_I(h)
    N = ss_downlink.calculate_N()
    I_N = ss_downlink.calculate_I_N(I, N)
    C_N, C = ss_downlink.calculate_C_N(N)
    C_I = ss_downlink.calculate_C_I(C, I)

    C_NI = ss_downlink.calculate_C_NI(C, I)
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
    print("I:", I, "dBW")
    print("N:", N, "dB")
    print("C:", C, "dBW")
    print("C/N (载波与噪声比): {:.2f} dB".format(C_N))
    print("C/I (载波与集总干扰比): {:.2f} dBW".format(C_I))
    print("C/(N+I) (载波与噪声+集总干扰比): {:.2f} dBW".format(C_NI))

    print("EPFD (等效载噪比密度): {:.2f} ".format(epfd))
    print("PFD : {:.2f} dBW/m2".format(pfd))
    print("delta_T/T: {:.2f} ".format(t))
    print("频谱效率: {:.2f} ".format(spectral_efficiency))
    print("信道容量: {:.2f} bps".format(channel_capacity))
