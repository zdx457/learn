import numpy as np
import pandas as pd


def create_gain_array_from_csv(file_path):
    # 读取 csv 文件
    df = pd.read_csv(file_path)
    # 提取 phi, theta 和 gain 列的数据
    phi = df['Phi[deg]'].values
    theta = df['Theta[deg]'].values
    gain = df['dB(RealizedGainRHCP)'].values
    # 将 phi 和 theta 的值乘以 2 并加上 360
    phi = (phi * 2 + 360)
    theta = (theta * 2 + 360)
    # 对 phi 和 theta 进行范围处理
    phi = np.clip(phi, 0, 720)
    theta = np.clip(theta, 0, 720)
    # 创建一个 721x721 的二维数组，初始化为 0
    gain_array = np.zeros((721, 721))
    # 使用高级索引将 gain 值存储在相应位置
    gain_array[theta.astype(int), phi.astype(int)] = gain
    return gain_array


def process_gain_array(file_path, row_index, col_index):
    # 创建增益数组
    gain_array = create_gain_array_from_csv(file_path)
    # 输出指定位置的值
    specified_value = gain_array[row_index, col_index]
    # 将 numpy 数组转换为 DataFrame
    df_gain = pd.DataFrame(gain_array)
    # 将 DataFrame 存储为 csv 文件，设置 header=False 去掉列标题
    df_gain.to_csv('new_gain_array.csv', index=False, header=False)
    return specified_value


if __name__ == "__main__":
    file_path = '2G-3D.csv'  # 替换为你的 csv 文件的实际路径
    row_index = 540  # phi*2+360
    col_index = 360  # seita*2+360
    specified_value = process_gain_array(file_path, row_index, col_index)
    print(f"指定位置 ({row_index}, {col_index}) 的值为: {specified_value}")