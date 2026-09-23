"""The atmospheres the sky solver computes, and their optical constants.

Numba-free so that other code (engine parameter exports, tests) can read the same
definitions the solver uses. solver.py imports everything here.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

# Binned ASTM G173 extraterrestrial irradiance and 233 K ozone cross sections,
# transcribed from Eric Bruneton's BSD-3-Clause 2017/2018 reference demo.
# Source details and complete upstream notice are in SOURCES.md/LICENSES.txt.
SW=np.arange(360.,831.,10.)
SF=np.array([1.11776,1.14259,1.01249,1.14716,1.72765,1.73054,1.6887,1.61253,1.91198,2.03474,2.02042,2.02212,1.93377,1.95809,1.91686,1.8298,1.8685,1.8931,1.85149,1.8504,1.8341,1.8345,1.8147,1.78158,1.7533,1.6965,1.68194,1.64654,1.6048,1.52143,1.55622,1.5113,1.474,1.4482,1.41018,1.36775,1.34188,1.31429,1.28303,1.26758,1.2367,1.2082,1.18737,1.14683,1.12362,1.1058,1.07124,1.04992])
O3=np.array([1.18e-27,2.182e-28,2.818e-28,6.636e-28,1.527e-27,2.763e-27,5.52e-27,8.451e-27,1.582e-26,2.316e-26,3.669e-26,4.924e-26,7.752e-26,9.016e-26,1.48e-25,1.602e-25,2.139e-25,2.755e-25,3.091e-25,3.5e-25,4.266e-25,4.672e-25,4.398e-25,4.701e-25,5.019e-25,4.305e-25,3.74e-25,3.215e-25,2.662e-25,2.238e-25,1.852e-25,1.473e-25,1.209e-25,9.423e-26,7.455e-26,6.566e-26,5.105e-26,4.15e-26,4.228e-26,3.237e-26,2.451e-26,2.801e-26,2.534e-26,1.624e-26,1.465e-26,2.078e-26,1.383e-26,7.105e-27])

@dataclass(frozen=True)
class Atmosphere:
    name:str
    radius_m:float
    scale_height_m:float
    density_scale:float
    ozone_du:float=300.
    ground_albedo:float=.1
    top_scale_heights:float=10.
    sun_radius_deg:float=32./120.
    visible_filter:float=1.
    @property
    def top(self): return self.radius_m+self.scale_height_m*self.top_scale_heights
    @property
    def ozone_scale(self): return self.scale_height_m/8000.

EARTH=Atmosphere('earth',6371000.,8000.,1.)
# Optical profile proxy: terrestrial effective scale height stretched by g_E/g_M.
# At the same profile temperature, 1.2x surface number density gives 7.26x column.
MOON=Atmosphere('moon',1737400.,8000.*9.80665/1.62,1.2)
MOON_ZERO=Atmosphere('moon_no_ozone',MOON.radius_m,MOON.scale_height_m,1.2,0.)


def rayleigh(lam_nm, atm: Atmosphere):
    """Molecular scattering coefficient at the surface (1/m), as in the solver."""
    return 1.24062e-6 * (np.asarray(lam_nm, dtype=float) / 1000.) ** -4 * atm.density_scale


def ozone_absorption(lam_nm, atm: Atmosphere):
    """Ozone absorption coefficient at the ozone layer's peak (1/m), as in the solver."""
    return np.interp(lam_nm, SW, O3) * (atm.ozone_du * 2.687e20 / (15000 * atm.ozone_scale))
