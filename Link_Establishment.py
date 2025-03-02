
import numpy as np
from constant import R_e

def establish_links(selected_ground_points, selected_satellites):
    """
    建立地面波位点与卫星之间的链路。

    参数：
    - selected_ground_points: 在原始区域内的地面波位点列表
    - selected_satellites: 在扩大区域内的卫星列表

    返回：
    - 链路对列表，每个元素为 (ground_index, satellite_id)
    """
    links = []
    for ground in selected_ground_points:
        ground_pos = np.array([ground.x_ecef, ground.y_ecef, ground.z_ecef])
        max_alpha = None
        selected_satellite_id = None

        for sat in selected_satellites:
            sat_pos = sat.position

            # 计算 d 和 r
            d = np.linalg.norm(sat_pos - ground_pos)
            r = np.linalg.norm(sat_pos)

            # 计算 α₁
            numerator = R_e**2 + d**2 - r**2
            denominator = 2 * R_e * d
            cos_alpha = numerator / denominator

            # 防止数值误差导致 cos_alpha 超出 [-1, 1]
            cos_alpha = np.clip(cos_alpha, -1.0, 1.0)

            alpha = np.degrees(np.arccos(cos_alpha))

            if alpha >= 90:
                if max_alpha is None or alpha > max_alpha:
                    max_alpha = alpha
                    selected_satellite_id = sat.satellite_id

        if selected_satellite_id is not None:
            link = (ground.index, selected_satellite_id)
            links.append(link)

    return links
