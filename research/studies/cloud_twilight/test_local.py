import json

import numpy as np
import pytest

from illumination.cloud_light.condensate import liquid_optical_depth
from climate.crm.cloud_columns import sha256
from shared.constants import SYNODIC_MONTH_DAYS
from .local import HERE, ROOT


def test_liquid_extinction_matches_monodisperse_cross_sections():
    # Independent sphere mass and Qext=2 projected area.
    n_per_m2, radius, density = 1e10, 12e-6, 997.
    mass = n_per_m2 * 4 / 3 * np.pi * radius ** 3 * density
    tau = n_per_m2 * 2 * np.pi * radius ** 2
    assert liquid_optical_depth(mass, radius, density) == pytest.approx(tau)
    with pytest.raises(ValueError):
        liquid_optical_depth(-1, radius, density)


def test_saved_local_clouds_preserve_spinup_boundary_and_mass_accounting():
    p = json.loads((HERE / 'results/local_clouds.json').read_text())
    assert p['source']['sampled_span_days'][0] >= SYNODIC_MONTH_DAYS
    assert p['source']['sampled_span_days'][1] <= 2 * SYNODIC_MONTH_DAYS
    assert p['checks']['maximum_mean_density_relative_error'] < 2e-5
    assert p['checks']['maximum_height_band_mass_closure_kg_m2'] < 1e-10
    assert sum(b['columns'] for b in p['local_time_bins']) == p['source']['columns']
    assert 0 <= p['cloud_column_fraction'] <= 1
    for key, values in p['species_path_kg_m2'].items():
        assert sum(b['species_water_path_kg_m2'][key]['mean'] for b in p['height_bands']) == pytest.approx(values['mean'])
    for path, expected in p['producer']['files'].items():
        assert sha256(ROOT / path) == expected
