"""An open sea's slopes from a friction velocity alone, where no wave run resolves its longer waves.

The sea appearance study (research/studies/sea_appearance/slopes.py) adds a wave run's resolved slopes to the
short-wave regime of the unified spectrum of Elfouhaily et al. (1997) above the run's highest wavenumber
(short_waves.py). Away from the wave runs' coasts there is no run, so open_sea_slopes stands the spectrum's long-wave
regime (31) of a fully developed sea in for it, up to a chosen wavenumber, and keeps the short-wave regime (40) above,
under the 10 m wind of the neutral surface layer with Charnock's roughness at the given gravity. Below the short waves'
threshold only the long waves tilt the water. In light winds the long-wave regime carries more slope than a wave run
gives: at the nearside coasts' median winds about a quarter more.
"""
from __future__ import annotations

from illumination.water_surface.short_waves import (FULLY_DEVELOPED, _peak, neutral_wind, short_wave_slopes,
                                                    slope_variances)


def neutral_friction_velocity(wind, gravity, charnock, roughness_floor, von_karman, height=10.0):
    """The friction velocity whose neutral surface layer gives the wind at the height, by bisection."""
    low, high = 1e-5, 5.0
    for _ in range(80):
        middle = (low + high) / 2
        at = neutral_wind(middle, gravity, charnock, roughness_floor, von_karman, height)
        low, high = (middle, high) if at < wind else (low, middle)
    return (low + high) / 2


def open_sea_slopes(u_star, gravity, split_wavenumber, charnock, roughness_floor, von_karman):
    """Mean square slope (the trace of the slope covariance) of an open sea under a friction velocity: the long waves
    of a fully developed sea up to split_wavenumber (rad/m) and the short waves above. Returns the total, its two parts
    and the 10 m wind."""
    u10 = float(neutral_wind(u_star, gravity, charnock, roughness_floor, von_karman))
    k_p, _ = _peak(u10, FULLY_DEVELOPED, gravity)
    long_waves = (slope_variances(u10, u_star, FULLY_DEVELOPED, gravity, k_to=split_wavenumber)['total']
                  if k_p / 10 < split_wavenumber else 0.0)
    short = short_wave_slopes(u_star, gravity, split_wavenumber)['total']
    return dict(total=long_waves + short, long_waves=long_waves, short_waves=short, wind_10m=u10)
