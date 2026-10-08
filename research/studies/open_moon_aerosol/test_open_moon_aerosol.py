"""Checks of the aerosol study: the particle machinery against single sizes and κ-Köhler theory, the night's decay, the
cases' order, fog, the sulfur arithmetic and the regions' areas."""
import math

import numpy as np
import pytest

from atmosphere.electricity import conductivity as cd
from research.studies.open_moon_aerosol import budget as bu
from research.studies.open_moon_aerosol import inputs as inp
from research.studies.open_moon_aerosol import particles as pa


def test_a_narrow_mode_takes_up_ions_as_one_size_does_and_holds_its_mass():
    m = pa.Mode(1.0e8, 1.0e-7, 1.01, 0.3)
    assert pa.ion_sink([m], 101325.0, 293.15, 0.0) == pytest.approx(cd.aerosol_attachment(0.05, 1e8, 101325.0, 293.15),
                                                                     rel=1e-3)
    wide = pa.Mode(2.0e8, 1.5e-7, 1.8, 0.2)
    assert pa.number_for_mass(wide.mass_kg_m3(), wide.d_m, wide.sigma) == pytest.approx(wide.number_m3)


def test_cloud_nuclei_follow_kappa_koehler_theory():
    t, s = 293.15, 0.003
    dc = pa.critical_diameter(0.3, s, t)
    a = 4 * cd.SURFACE_TENSION * cd.WATER_MOLAR_MASS / (cd.GAS_CONSTANT * t * cd.WATER_DENSITY)
    assert 0.3 == pytest.approx(4 * a ** 3 / (27 * dc ** 3 * math.log(1 + s) ** 2))           # Petters and Kreidenweis Eq. 10
    assert pa.ccn([pa.Mode(1e8, 5 * dc, 1.2, 0.3)], s, t) == pytest.approx(1e8, rel=1e-6)
    assert pa.ccn([pa.Mode(1e8, dc / 5, 1.2, 0.3)], s, t) < 1.0
    assert pa.ccn([pa.Mode(1e8, dc, 1.5, 0.3)], s, t) == pytest.approx(0.5e8)


def test_the_night_decays_toward_its_end_and_a_layer_fills_from_its_flux():
    assert bu.night_mean(1.0) == 1.0
    assert bu.night_mean(math.exp(-2.0)) == pytest.approx((1 - math.exp(-2.0)) / 2.0)
    assert pa.night_layer(0.0, 1.0, 200.0, 0.0, 10.0) == pytest.approx(10.0 * 3600.0 / 200.0 / 2.0)
    steady = pa.night_layer(0.0, 1.0, 200.0, 1e-3, 1000.0)
    assert steady == pytest.approx(1.0 / (200.0 * 1e-3), rel=1e-3)


def test_the_cases_run_from_the_fewest_particles_to_the_most_and_fog_takes_the_ions():
    for region in ('seas', 'wet_land', 'fog_desert', 'polar_dry_land'):
        for phase in ('day', 'night'):
            n = [sum(m.number_m3 for m in bu.modes(region, case, phase, 300.0)) for case in bu.CASES]
            assert n[0] < n[1] < n[2], (region, phase)
    modes = bu.modes('wet_land', 'central', 'day', 400.0)
    clear = pa.conductivity(modes, 0.04e6, 121000.0, 298.0, 0.8)
    mu = cd.mobility_air(cd.MOBILITY_STANDARD[0], 121000.0, 298.0)
    fog = cd.droplet_attachment(298.0, mu, 60e6, 8e-6) + pa.ion_sink(modes, 121000.0, 298.0, 0.99)
    assert cd.conductivity_of(0.04e6, 121000.0, 298.0, fog)[0] < 0.1 * clear


def test_the_seas_sulfur_reaches_the_plains_at_parts_per_trillion():
    # 1.2 umol/m2/day over 28 % of the Moon, carried onto 62 %, 70 % as sulfur dioxide, 1.5 days in a 5-km layer
    flux = 1.2e-6 * 0.7 * 0.28 / 0.62
    assert bu.sulfur_dioxide_ppt(1.2, 0.28, 0.62, 5000.0) == pytest.approx(flux * 1.5 / 5000.0 / 49.0 * 1e12)
    assert bu.sulfur_dioxide_ppt(3.9, 0.28, 0.62, 5000.0) < 700.0 / 50.0       # far below the savanna's 0.7 ppb


def test_the_regions_cover_the_moon():
    if not inp.CLIMATOLOGY.exists():
        pytest.skip('the GCM climatology product is missing (python -m climate.gcm.climatology A28_dim5_moon:20-29)')
    assert sum(inp.areas().values()) == pytest.approx(1.0)
    w = inp.regime('fog_desert')
    assert w['night']['relative_humidity'] > 0.9 and w['day']['rain_mm_day'] < 0.5
