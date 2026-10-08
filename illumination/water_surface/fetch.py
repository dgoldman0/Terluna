"""Conditional fetch-to-wave-age closure from Elfouhaily et al. (1997), eq. 37.

X = g F / U10**2. This is an Earth empirical closure, here evaluated at supplied
gravity; it is not a terrain-resolving wave forecast or a duration-growth model.
Source: https://doi.org/10.1029/97JC00467, p. 15788.
"""
import numpy as np


def inverse_wave_age(fetch_m, wind_10m_m_s, gravity):
    values = np.asarray([fetch_m, wind_10m_m_s, gravity], dtype=float)
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError('Fetch, wind speed and gravity must be finite and positive')
    x = gravity*fetch_m/wind_10m_m_s**2
    return float(.84*np.tanh((x/22000.)**.4)**(-.75))
