import pandas as pd
from antenna_satellite import *
import numpy as np

from tool import calculate_angle


class Channel_h:
    def __init__(self, f, T, rx_antenna_gain, beam_center, coordinate_s, direction_vertices, coordinate_user):     #
        self.f = f  # 通信卫星的频率，单位为Hz
        self.k = 1.38e-23  # 玻尔兹曼常数，单位：J/K
        self.T = T  # 为接收机噪声温度  单位：K
        self.rx_antenna_gain=rx_antenna_gain
        #损耗参数
        self.daqi_loss=0.1  #dB
        self.shanshuo_loss = 0.3  # dB
        self.qujihua_loss = 3  # dB
        #坐标参数
        self.coordinate_s=np.array(coordinate_s)  #单位m
        self.direction_vertices=np.array(direction_vertices)  #运动方向向量
        self.coordinate_user=np.array(coordinate_user)
        self.beam_center=beam_center

        self.s_to_u=  self.coordinate_user - self.coordinate_s # 卫星指向用户
        self.s_to_b=  self.beam_center - self.coordinate_s # 卫星指向波束中心
        self.distance=np.linalg.norm(self.s_to_u) # m 用户和卫星的距离
        self.gain_array = None

    def calculate_channel_gain(self):
        """
        计算信道增益 G
        """
        D = 4.535 # 天线口径 m

        angle_off_axis =  self.calculate_angle_v1_v2(self.s_to_b,self.s_to_u)#离轴角
        tx_antenna_gain = NGSO_SAT_1528_1(D, self.f, angle_off_axis)[0]
        FSPL = 32.45 + 20 * math.log10(self.distance / 1e3) + 20 * math.log10(self.f / 1e6)
        G = tx_antenna_gain + self.rx_antenna_gain - FSPL - self.qujihua_loss
        return G,angle_off_axis

    def calculate_channel_gain_from_table(self):
        Theta, Phi = self.calculate_angle_theta_and_phi(self.beam_center, self.coordinate_s, self.direction_vertices, self.coordinate_user)
        tx_antenna_gain = self.find_gain('2G-3D.csv', Theta, Phi)
        FSPL = 32.45 + 20 * math.log10(self.distance / 1e3) + 20 * math.log10(self.f / 1e6)
        G = tx_antenna_gain + self.rx_antenna_gain - FSPL - self.daqi_loss - self.shanshuo_loss - self.qujihua_loss
        return G


    def calculate_angle_theta_and_phi(self, beam_center, coordinate_s, direction_vertices, coordinate_user,):
        """
        俯仰角 theta
        方位角 phi
        """
        v = direction_vertices # 卫星运动方向
        x = (beam_center - coordinate_s) / np.linalg.norm(beam_center - coordinate_s) # x方向
        z = (v-np.dot(v,x)*x)/np.linalg.norm(v-np.dot(v,x)*x) # z方向
        y = np.cross(x,z)
        vertice = coordinate_user - coordinate_s  # u-s

        theta,phi= calculate_angle(x,y,z,coordinate_s,coordinate_user)
        return theta, phi # 角度

    def find_gain(self, path, theta, phi):
        if self.gain_array is None:
            try:
                # 尝试从 .npy 文件加载数组
                self.gain_array = np.load('gain_array.npy')
            except FileNotFoundError:

                # 如果文件不存在，从 CSV 文件生成数组
                df = pd.read_csv(path)
                phi_vals = df['Phi[deg]'].values
                theta_vals = df['Theta[deg]'].values
                gain = df['dB(RealizedGainRHCP)'].values
                phi_vals = (phi_vals * 2 + 360)
                theta_vals = (theta_vals * 2 + 360)
                phi_vals = np.clip(phi_vals, 0, 720)
                theta_vals = np.clip(theta_vals, 0, 720)
                self.gain_array = np.zeros((721, 721))
                self.gain_array[phi_vals.astype(int), theta_vals.astype(int)] = gain
                # 保存数组到 .npy 文件
                np.save('gain_array.npy', self.gain_array)


        phi_new = int(phi * 2 + 360)
        theta_new = int(theta * 2 + 360)
        gain = self.gain_array[phi_new, theta_new]
        return gain

    def calculate_angle_v1_v2(self, vector1, vector2):
        # 将输入转换为 numpy 数组
        vector1 = np.array(vector1)
        vector2 = np.array(vector2)
        # 计算向量的点积
        dot_product = np.dot(vector1, vector2)
        # 计算向量的模
        norm_vector1 = np.linalg.norm(vector1)
        norm_vector2 = np.linalg.norm(vector2)
        # 计算夹角的余弦值
        cos_angle = dot_product / (norm_vector1 * norm_vector2)
        # 修正可能溢出的值
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        # 计算夹角（弧度制）
        angle = np.arccos(cos_angle)
        # 将弧度转换为角度
        angle = np.degrees(angle)
        return angle



if __name__=='__main__':
    # P=10**(44.3/10)/1000          # 单波束功率44.3dBm
    # # 路损计算
    # f = 2 * 1e9
    # h = 600 * 1e3
    # N = 2  # 波束层数
    # cell_radius = 12.5 * 1e3  # 波束半径m
    # c = np.array([0, 100, 200])
    # s1 = np.array([0, 0, 600 * 1e3])
    # vertices = np.array([0, 1, 0])
    # u1 = np.array([0,0,0])
    #
    # channel = Channel_h(f,300,5,c,s1,vertices,u1)
    # G = channel.calculate_channel_gain()

    c = Channel_h(None,None,None,0,0,None,0)

    print(c.find_gain('2G-3D.csv',90,0))
    print(c.find_gain('2G-3D.csv', -89,-179))
    print(c.find_gain('2G-3D.csv', 89, 1))
    print(c.find_gain('2G-3D.csv', 90, 0))