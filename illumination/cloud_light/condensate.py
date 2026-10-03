"""Visible liquid-cloud extinction for prescribed effective-radius scenarios."""
import numpy as np


def liquid_optical_depth(liquid_water_path_kg_m2, effective_radius_m, liquid_density_kg_m3):
    """Geometric-optics Qext=2, uniform prescribed effective radius.

    This converts cloud liquid alone. Ice habits, precipitation, spectral Mie
    structure and the scattering source function require separate treatment.
    """
    path = np.asarray(liquid_water_path_kg_m2, float)
    if (not np.isfinite(path).all() or np.any(path < 0) or
            not np.isfinite([effective_radius_m, liquid_density_kg_m3]).all() or
            min(effective_radius_m, liquid_density_kg_m3) <= 0):
        raise ValueError('Nonnegative water path and positive radius and density required')
    return 3 * path / (2 * liquid_density_kg_m3 * effective_radius_m)
