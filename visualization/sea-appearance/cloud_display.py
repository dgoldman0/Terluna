"""Declared angular filtering of noisy cloud radiance for display references.

This is a rendering approximation, not a new transport calculation. Filter a
regular one-degree elevation/azimuth grid so irregular sampling does not turn
individual photon outliers into coloured stripes. Raw groups remain retained.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter


def delta_interpolator(elevation, azimuth, cloudy, clear, sigma_deg=0.):
    delta = cloudy-clear
    def interpolator(e, a, values):
        periodic = np.r_[a[-1]-360,a,a[0]+360]
        field = np.concatenate([values[:,-1:],values,values[:,:1]],axis=1)
        return RegularGridInterpolator((e,periodic),field,bounds_error=False,fill_value=None)
    if sigma_deg <= 0:
        return interpolator(elevation,azimuth,gaussian_filter(delta,(.6,.6,0)))
    # No finer-than-model morphology is implied by this display grid.
    ee = np.linspace(elevation[0],elevation[-1],round(elevation[-1]-elevation[0])+1)
    aa = np.arange(360.)
    e,a = np.meshgrid(ee,aa,indexing='ij')
    field = interpolator(elevation,azimuth,delta)(np.stack([e,a],-1))
    field = gaussian_filter(field,(sigma_deg/(ee[1]-ee[0]),sigma_deg,0),mode=('nearest','wrap','nearest'))
    return interpolator(ee,aa,field)
