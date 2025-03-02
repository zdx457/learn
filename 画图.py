import pandas as pd
import matplotlib.pyplot as plt

# 设置图片清晰度
# plt.rcParams['figure.dpi'] = 100

# 设置全局字体大小
plt.rcParams.update({'font.size': 16})

# 读取 Excel 文件
file_path = '倾角对应距离.xlsx'
df = pd.read_excel(file_path)

# 提取倾角和距离数据
angles = df['angle']
distances = df['distance']

# 绘制倾角和距离的二维图
plt.figure(figsize=(10, 6))
plt.plot(angles, distances, c='r',marker='.', linestyle='-')
# 设置标题，添加字体设置以避免可能的字体问题
plt.rcParams['font.sans-serif'] = ['SimSun']  # 设置中文字体为黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
plt.title('不同倾角所对应距离的二维图')
plt.xlabel('仰角( °)')
plt.xticks(rotation=45)
plt.ylabel('隔离距离 (m)')
plt.grid(True)

# 显示图形
plt.show()
