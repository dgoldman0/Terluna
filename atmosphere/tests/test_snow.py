"""The snow module's physics: saturation, fall speed, growth, melting and the winter column."""
import math

import pytest

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from atmosphere.radiative_convective.thermodynamics import DryAir
from atmosphere.snow import flakes as sf


def air_at_freezing(pressure=100000.0):
    column = sf.Column('test', STANDARD_GRAVITY, DryAir(101325.0), 0.0, pressure, 6.5e-3)
    return column, column.state(0.0)


def test_saturation_over_ice_and_liquid_meet_at_the_triple_point_and_part_below():
    assert sf.ice_saturation_pressure(273.16) == pytest.approx(611.657, rel=2e-3)
    assert sf.liquid_saturation_pressure(273.16) == pytest.approx(611.657, rel=2e-3)
    ratio = sf.liquid_saturation_pressure(258.15) / sf.ice_saturation_pressure(258.15)
    assert 1.14 < ratio < 1.17                                          # the dendritic zone's supersaturation


def test_a_centimetre_flake_falls_as_heymsfield_and_westbrook_give_and_slower_at_lunar_gravity():
    _, s = air_at_freezing()
    d = 0.01
    m = sf.flake_mass(d)
    x = 8.0 * m * STANDARD_GRAVITY * s['rho'] / (math.pi * math.sqrt(0.35) * s['eta'] ** 2)
    re = 16.0 * (math.sqrt(1.0 + 4.0 * math.sqrt(x) / (64.0 * math.sqrt(0.35))) - 1.0) ** 2
    earth, _ = sf.fall_speed(d, m, STANDARD_GRAVITY, s)
    assert earth == pytest.approx(s['eta'] * re / (s['rho'] * d), rel=1e-12)
    assert 1.2 < earth < 1.5                                            # large aggregates fall at 1-1.5 m/s
    moon, _ = sf.fall_speed(d, m, MOON_SURFACE_GRAVITY, s)
    assert 0.35 < moon / earth < 0.5


def test_ice_grows_at_liquid_saturation_below_freezing_and_holds_at_ice_saturation():
    column = sf.Column('test', STANDARD_GRAVITY, DryAir(101325.0), 0.0, 80000.0, 6.5e-3)
    s = column.state(column.height_of(258.15))
    assert sf.deposition_rate(1e-3, 10.0, s) > 0.0
    s_ice = dict(s, e=sf.ice_saturation_pressure(s['t']))
    assert sf.deposition_rate(1e-3, 10.0, s_ice) == pytest.approx(0.0, abs=1e-18)


def test_melting_starts_only_where_the_air_gives_the_flake_heat():
    column, s = air_at_freezing()
    assert sf.melting_rate(5e-3, 50.0, s) == pytest.approx(0.0, abs=1e-15)          # saturated air at 0 C
    warm = column.state(-500.0, relative_humidity=1.0)                               # 3.25 K warmer, saturated
    assert sf.melting_rate(5e-3, 50.0, warm) > 0.0
    dry = column.state(-100.0, relative_humidity=0.5)                                # its evaporation cools it more
    assert sf.melting_rate(5e-3, 50.0, dry) == 0.0


def test_the_column_is_hydrostatic_through_its_freezing_level():
    column = sf.Column('test', MOON_SURFACE_GRAVITY, DryAir(101325.0), 29000.0, 66000.0, 0.92e-3)
    assert column.temperature(29000.0) == pytest.approx(sf.T0)
    assert column.height_of(sf.T0 - 9.2) == pytest.approx(39000.0)
    assert column.state(30000.0)['p'] < 66000.0 < column.state(28000.0)['p']


def test_flakes_melt_deeper_below_the_freezing_level_in_drier_air():
    column = sf.Column('test', STANDARD_GRAVITY, DryAir(101325.0), 1500.0, 85000.0, 6.5e-3)
    moist, dry = sf.melt(column, 5.0, 0.9), sf.melt(column, 5.0, 0.7)
    assert 0 < moist['half_melted_below_m'] < moist['melted_below_m'] < dry['melted_below_m'] < 1500.0


def test_the_aggregation_models_laws_match_the_single_flake_laws():
    import numpy as np
    from atmosphere.snow import aggregation as ag
    column = sf.Column('test', MOON_SURFACE_GRAVITY, DryAir(101325.0), 29000.0, 66000.0, 0.92e-3)
    s = column.state(35000.0)
    d = np.array([3e-4, 2e-3, 1e-2])
    m = sf.flake_mass(d)
    u, re = ag.speeds(d, m, MOON_SURFACE_GRAVITY, s)
    for i in range(d.size):
        one, one_re = sf.fall_speed(d[i], m[i], MOON_SURFACE_GRAVITY, s)
        assert u[i] == pytest.approx(one, rel=1e-12) and re[i] == pytest.approx(one_re, rel=1e-12)
        single = sf.deposition_rate(d[i], re[i], s)
        assert ag.deposition(d[i:i + 1], re[i:i + 1], s)[0] == pytest.approx(single, rel=1e-12)


def test_the_breakup_law_follows_phillips_yano_and_khain_for_snow():
    import numpy as np
    from atmosphere.snow import aggregation as ag
    d = 1e-3
    small = sf.flake_mass(d)
    energy = 1e-9
    scale = 1.58e7 * (1.0 + 1.33e-4 / d ** 1.5) * math.pi * d * d
    expected = scale * (1.0 - math.exp(-(7.08e6 * 3.5e-3 * energy / scale) ** 0.5))
    assert ag.fragments(np.array([small]), np.array([energy]))[0] == pytest.approx(min(expected, 100.0), rel=1e-6)
    assert ag.fragments(np.array([small]), np.array([0.0]))[0] == 0.0
    assert ag.fragments(np.array([small]), np.array([1.0]))[0] == 100.0


def test_meetings_keep_the_snows_mass_whether_they_stick_or_break():
    import numpy as np
    from atmosphere.snow import aggregation as ag
    rng = np.random.default_rng(1)
    d = np.exp(rng.uniform(math.log(2e-4), math.log(5e-3), 400))
    m = sf.flake_mass(d)
    xi = np.round(rng.uniform(1e3, 1e5, 400))
    u = 0.1 + 0.5 * (d / d.max()) ** 0.5
    m_frag = sf.flake_mass(1e-4)
    for temperature, breakup in ((273.15, False), (253.15, True)):
        mm, xx = m.copy(), xi.copy()
        before = float(np.sum(mm * xx))
        made = ag.collide(mm, xx, u, d, np.arange(d.size), 100.0, temperature, 30.0, rng, breakup, m_frag)
        assert float(np.sum(mm * xx)) + made * m_frag == pytest.approx(before, rel=1e-9)
        assert (mm > 0).all() and (xx >= 0).all()


def test_a_short_column_brings_snow_to_the_freezing_level():
    from atmosphere.snow import aggregation as ag
    column = sf.Column('test', STANDARD_GRAVITY, DryAir(101325.0), 1500.0, 85000.0, 6.5e-3)
    out = ag.summary(ag.column(column, -5.0, 1e4, 1.0, bins=6, inject=4, fill=1.0, gather=0.5))
    assert out['snowfall_mm_h'] > 0.0 and out['median_by_mass_mm'] >= 0.1
