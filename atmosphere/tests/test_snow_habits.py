"""Size conventions, projected area, humidity and explicit mechanical load bounds."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import flakes as sf
from atmosphere.snow import habits


def test_dimension_conversion_matches_hogan_instead_of_relabeling_mean_dimension():
    old = habits.GEOMETRIES['legacy_bf95']
    corrected = habits.GEOMETRIES['bf95_dmax']
    dmax = np.array([1e-4, 1e-3, 1e-2, 0.6])
    # Hogan's reported Dmax ~= 1.25 Dmean, with a rounded SI coefficient.
    assert corrected.mass(dmax) == pytest.approx(old.mass(dmax / 1.25), rel=0.002)
    assert corrected.size(old.mass(dmax)) / dmax == pytest.approx(1.25, rel=0.002)
    assert corrected.size(corrected.mass(dmax)) == pytest.approx(dmax, rel=1e-12)
    assert corrected.mass(0.15) == pytest.approx(0.0003291234, rel=1e-6)


def test_legacy_geometry_preserves_the_recorded_baseline():
    g = habits.GEOMETRIES['legacy_bf95']
    d = np.array([0.0002, 0.01, 0.2])
    assert g.mass(d) == pytest.approx(sf.flake_mass(d), rel=1e-12)
    assert g.area(d) == pytest.approx(0.35 * math.pi / 4 * d ** 2)


def test_ellipse_area_is_converted_to_circle_area():
    out = habits.aggregate_area_geometry('plates', 1024)
    # Multiplying the ellipse's area by q must match multiplying the Dmax
    # circle by the corrected area ratio. Using q directly would fail this.
    dmax = 0.01
    ellipse_area = math.pi / 4 * dmax * (dmax * out['horizontal_aspect_ratio'])
    circle_area = math.pi / 4 * dmax ** 2
    assert out['circle_area_ratio'] * circle_area == pytest.approx(out['ellipse_area_ratio'] * ellipse_area)
    assert out['circle_area_ratio'] == pytest.approx(0.3128188564, rel=1e-9)
    assert out['circle_area_ratio'] < out['ellipse_area_ratio']


def test_geometry_regressions_reject_unsupported_monomer_extrapolation():
    for invalid in (1, 2049, 4.5):
        with pytest.raises(ValueError):
            habits.aggregate_area_geometry('needles', invalid)
    with pytest.raises(ValueError):
        habits.aggregate_area_geometry('graupel', 64)


def test_ice_growth_can_occur_below_water_saturation():
    assert habits.humidity_regime(-15, 0.9)['supersaturation_ice'] > 0
    assert habits.humidity_regime(-15, 0.8)['regime'] == 'ice sublimation'
    assert habits.humidity_regime(-15, 1.0)['supersaturation_ice'] == pytest.approx(0.15741745, rel=1e-6)


def test_temperature_anchors_follow_column_lapse_without_assigning_habit_fractions():
    from atmosphere.radiative_convective.thermodynamics import DryAir
    from shared.constants import MOON_SURFACE_GRAVITY
    column = sf.Column('test', MOON_SURFACE_GRAVITY, DryAir(101325), 29000, 66000, 0.001)
    rows = habits.habit_ladder(column)
    assert rows[0]['height_m'] == pytest.approx(44000)
    assert rows[-1]['height_m'] == pytest.approx(31000)
    assert rows[0]['pressure_pa'] < rows[-1]['pressure_pa']
    assert all('fraction' not in row for row in rows)


def test_reference_load_uses_velocity_difference_and_required_strength_is_explicit():
    inputs = dict(diameter_m=0.6, air_density_kg_m3=0.7, area_ratio=0.35)
    zero = habits.differential_load(**inputs, velocity_difference_m_s=0)
    assert zero['reference_differential_force_n'] == 0
    at_point_one = habits.differential_load(**inputs, velocity_difference_m_s=0.1)
    assert at_point_one['reference_differential_force_n'] == pytest.approx(0.0003463605901)
    at_point_two = habits.differential_load(**inputs, velocity_difference_m_s=0.2)
    assert at_point_two['reference_differential_force_n'] == pytest.approx(
        4 * at_point_one['reference_differential_force_n'])
    assert habits.critical_velocity_difference(at_point_one['reference_differential_force_n'], **inputs) == pytest.approx(0.1)


def test_geometries_have_resolvable_source_ids_and_no_claimed_strength():
    ledger = json.loads((Path(habits.__file__).with_name('habits_sources.json')).read_text())
    sources = {source['id'] for source in ledger['sources']}
    for g in habits.GEOMETRIES.values():
        assert set(g.source_ids) <= sources
        assert 'not validated' in g.evidence
        assert 'strength' not in g.to_dict()
