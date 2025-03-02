# satellite.py
import numpy as np
from Walker_convert import orbital_elements_to_ECI, ECI_to_LLH

class Satellite:
    def __init__(self, satellite_id, orbital_elements, wavelength, pt):
        self.satellite_id = satellite_id
        self.orbital_elements = orbital_elements  # 字典形式，包含 a, e, i, Omega, omega, f
        self.position = None  # ECEF 坐标(x,y,z)
        self.wavelength = wavelength
        self.pt = pt
        self.lat = None
        self.lon = None
        self.alt = None
        self.velocity = None

    def calculate_velocity(self, time):
        n = self.orbital_elements['n']  # 平均运动
        delta_M = n * time
        f = (self.orbital_elements['f'] + delta_M) % (2 * np.pi)
        orbital_elements_1 = {
            'a': self.orbital_elements['a'],
            'e': self.orbital_elements['e'],
            'i': self.orbital_elements['i'],
            'Omega': self.orbital_elements['Omega'],
            'omega': self.orbital_elements['omega'],
            'f': f,
            'n': self.orbital_elements['n']
        }
        r_ecef = orbital_elements_to_ECI(orbital_elements_1)
        self.velocity = r_ecef - self.position

    def update_position(self, time):
        # 更新真近点角
        n = self.orbital_elements['n']  # 平均运动
        delta_M = n * time
        f = (self.orbital_elements['f'] + delta_M) % (2 * np.pi)
        self.orbital_elements['f'] = f

        # 计算 ECI 坐标
        r_eci = orbital_elements_to_ECI(self.orbital_elements)

        # 在不考虑地球自转的情况下，ECI 坐标等同于 ECEF 坐标
        r_ecef = r_eci

        self.position = r_ecef
        self.lat, self.lon, self.alt = ECI_to_LLH(r_ecef)
