# ground_point.py

class GroundPoint:
    def __init__(self, index, x_ecef, y_ecef, z_ecef, lat, lon, alt):
        """
        初始化地面波位点对象。

        参数：
        - index: 波位点索引
        - x_ecef, y_ecef, z_ecef: 波位点的 ECEF 坐标（km）
        - lat: 纬度（度）
        - lon: 经度（度）
        """
        self.index = index
        self.x_ecef = x_ecef
        self.y_ecef = y_ecef
        self.z_ecef = z_ecef
        self.lat = lat
        self.lon = lon
        self.alt = alt
        # self.linked_satellite_id = None  # 当前时间步建立链路的卫星编号

    # def set_link(self, satellite_id):
    #     """
    #     设置波位点的链路信息。
    #
    #     参数：
    #     - satellite_id: 链接的卫星编号
    #     """
    #     self.linked_satellite_id = satellite_id
    #
    # def clear_link(self):
    #     """
    #     清除波位点的链路信息。
    #     """
    #     self.linked_satellite_id = None

