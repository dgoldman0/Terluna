"""The single-drop solver: Earth's measured fall speeds, the physical limits, the recovered solver's outputs, and its
numerical bookkeeping."""
import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from atmosphere.radiative_convective.thermodynamics import DryAir, earthlike_air
from atmosphere.rain import drops

RAIN = Path(drops.__file__).resolve().parent
VENDOR = RAIN / 'vendor' / 'rainprops'
EARTH_AIR = DryAir(101325.0)
# The recovered solver's own inputs (11 October 2026): the thunder product's lunar mixture and its Earth mixture.
RECOVERED_MOON_AIR = DryAir(121590.0, {'N2': 0.815269668774136, 'O2': 0.1745785015215067,
                                       'Ar': 0.009751829704357395, 'CO2': 0.0004})
RECOVERED_EARTH_AIR = DryAir(101325.0, {'N2': 0.78084, 'O2': 0.20946, 'Ar': 0.00934, 'CO2': 0.00036})


def column(gravity=STANDARD_GRAVITY, pressure=101325.0, temperature=293.15, humidity=0.5, air=EARTH_AIR, **kw):
    return drops.Column(gravity_m_s2=gravity, surface_pressure_pa=pressure, surface_temperature_k=temperature,
                        surface_relative_humidity=humidity, air=air, **kw)


def moon(temperature=293.15, humidity=0.8, air=None, **kw):
    return column(MOON_SURFACE_GRAVITY, 121590.0, temperature, humidity, air or earthlike_air(121590.0, 400.0), **kw)


def test_vendored_rainprops_matches_its_provenance():
    provenance = json.loads((VENDOR / 'PROVENANCE.json').read_text())
    assert provenance['license'] == 'MIT'
    assert provenance['commit'] == 'e2df9c2102112526f05109051dfffeeb5c318967'
    for entry in provenance['files']:
        data = (VENDOR / entry['path']).read_bytes()
        assert len(data) == entry['bytes'], entry['path']
        assert hashlib.sha256(data).hexdigest() == entry['sha256'], entry['path']
    assert 'MIT License' in (VENDOR / 'LICENSE').read_text()


def test_gunn_kinzer_earth_fall_speeds():
    """34 measured fall speeds, 0.1-5.8 mm, in air at 20 C, 1 atm and 50% humidity."""
    radius, measured = np.loadtxt(VENDOR / 'data' / 'gunn1949_r_m_v_mpers.csv', delimiter=',').T
    earth = column()
    calculated = np.array([drops.terminal_speed(2 * r, earth) for r in radius])
    assert len(calculated) == 34
    assert np.sqrt(np.mean((calculated - measured)**2)) < 0.205
    assert np.max(np.abs(calculated / measured - 1)) < 0.075
    assert np.all(np.diff(calculated) > 0)


@pytest.mark.parametrize('world, recovered', [
    ('earth', [0.7119406497, 2.0237887842, 3.7766117115, 6.2743207389, 7.7794554003, 8.6029302799, 9.0059685552,
               9.1640102693]),
    ('moon', [0.1575519060, 0.5560916911, 1.1489825144, 2.1126800694, 2.8838463249, 3.5001127756, 3.9829318068,
              4.3540459128]),
    ('moon_air_earth_gravity', [0.6850381373, 1.9072598670, 3.5306803578, 5.8245582241, 7.1809617409,
                                7.9066144923, 8.2515721292, 8.3783521213]),
])
def test_port_reproduces_the_recovered_terminal_speeds(world, recovered):
    """The recovered run's terminal-speeds.json (20 C, 80%): within 0.1%, the port's surface tension and saturation
    pressure being the only changes."""
    c = {'earth': column(humidity=0.8, air=RECOVERED_EARTH_AIR),
         'moon': moon(air=RECOVERED_MOON_AIR),
         'moon_air_earth_gravity': column(STANDARD_GRAVITY, 121590.0, 293.15, 0.8, RECOVERED_MOON_AIR)}[world]
    for mm, v in zip((0.2, 0.5, 1, 2, 3, 4, 5, 6), recovered):
        assert drops.terminal_speed(mm * 1e-3, c) == pytest.approx(v, rel=1e-3)


def test_port_reproduces_the_recovered_representative_drops():
    """The recovered representative-drops.json: cloud drops from the central lunar column and drips from rest."""
    central = moon(303.106, 0.765, RECOVERED_MOON_AIR)
    assert central.cloud_base_m == pytest.approx(3629.772, rel=1e-4)
    assert drops.fall(1e-3, central, sample_heights_m=[central.cloud_base_m, 0])['status'] == 'evaporated'
    for mm, diameter_mm, speed, time in ((1.5, 0.902656, 1.048099, 2438.409), (2.0, 1.620508, 1.790074, 1775.958),
                                         (5.0, 4.856406, 3.975390, 892.904)):
        r = drops.fall(mm * 1e-3, central, sample_heights_m=[central.cloud_base_m, 0])
        assert r['survived']
        assert r['surface_diameter_m'] * 1e3 == pytest.approx(diameter_mm, rel=1e-4)
        assert r['surface_speed_m_s'] == pytest.approx(speed, rel=1e-3)
        assert r['fall_time_s'] == pytest.approx(time, rel=5e-4)
    for height, mm, speed in ((0.1, 1.0, 0.514884), (1.0, 3.0, 1.595214), (3.0, 3.0, 2.336950)):
        r = drops.fall(mm * 1e-3, replace(central, cloud_base_m=height), initial_speed_m_s=0.0,
                       sample_heights_m=[height, 0])
        assert r['surface_speed_m_s'] == pytest.approx(speed, rel=1e-3)


def test_stokes_limit():
    """Tiny drops fall at (rho_w - rho_a) g D^2 / 18 mu, so as g and as D^2."""
    for c in (column(), moon()):
        state = c.state(0.0)
        for d in (2e-6, 5e-6, 1e-5):
            stokes = ((drops.WATER_DENSITY - state['density_kg_m3']) * c.gravity_m_s2 * d**2
                      / (18 * state['viscosity_pa_s']))
            assert drops.terminal_speed(d, c) == pytest.approx(stokes, rel=5e-3)
        assert drops.terminal_speed(1e-5, c) / drops.terminal_speed(5e-6, c) == pytest.approx(4.0, rel=5e-3)
    same_air = column(STANDARD_GRAVITY, 121590.0, 293.15, 0.8, earthlike_air(121590.0, 400.0))
    ratio = drops.terminal_speed(5e-6, same_air) / drops.terminal_speed(5e-6, moon())
    assert ratio == pytest.approx(STANDARD_GRAVITY / MOON_SURFACE_GRAVITY, rel=5e-3)


def test_inertial_regime_scalings():
    """A rigid sphere's speed at large Reynolds number follows sqrt(g) and nearly rho_a^(-1/2); deformation, which
    lunar gravity weakens, slows large drops less than that."""
    lunar = moon()
    same_air = replace(lunar, gravity_m_s2=STANDARD_GRAVITY, lapse_k_m=None, cloud_base_m=None)
    d = 6e-3
    rigid = drops.terminal_speed(d, same_air, deform=False) / drops.terminal_speed(d, lunar, deform=False)
    assert rigid == pytest.approx(math.sqrt(STANDARD_GRAVITY / MOON_SURFACE_GRAVITY), rel=0.02)
    deformable = drops.terminal_speed(d, same_air) / drops.terminal_speed(d, lunar)
    assert deformable < 0.8 * rigid
    dense = replace(lunar, surface_pressure_pa=2 * 121590.0, lapse_k_m=None, cloud_base_m=None)
    exponent = math.log(drops.terminal_speed(d, dense, deform=False) / drops.terminal_speed(d, lunar, deform=False))
    assert -0.5 < exponent / math.log(2.0) < -0.4


def test_slowing_by_gravity_spans_the_regimes():
    """With gravity alone changed in the Moon's air at 20 C and 80%, the slowing falls with size through CM1's single
    factor; the review's rerun of the recovered solver gave 4.35, 3.43, 3.07, 2.76, 2.49 and 1.93."""
    lunar = moon(air=RECOVERED_MOON_AIR)
    same_air = replace(lunar, gravity_m_s2=STANDARD_GRAVITY, lapse_k_m=None, cloud_base_m=None)
    slowing = [drops.terminal_speed(mm * 1e-3, same_air) / drops.terminal_speed(mm * 1e-3, lunar)
               for mm in (0.2, 0.5, 1, 2, 3, 6)]
    assert np.round(slowing, 2).tolist() == [4.35, 3.43, 3.07, 2.76, 2.49, 1.93]
    cm1 = (9.81 / MOON_SURFACE_GRAVITY)**0.6
    assert slowing[2] > cm1 > slowing[3]


def test_stability_screen_scales_as_inverse_root_gravity():
    earth, lunar = column(STANDARD_GRAVITY, 121590.0, 293.15, 0.8, RECOVERED_MOON_AIR), moon(air=RECOVERED_MOON_AIR)
    ratio = drops.stable_diameter(lunar) / drops.stable_diameter(earth)
    assert ratio == pytest.approx(math.sqrt(STANDARD_GRAVITY / MOON_SURFACE_GRAVITY), rel=1e-6)
    assert round(ratio, 2) == 2.46
    assert drops.stable_diameter(column()) * 1e3 == pytest.approx(6.68, abs=0.01)


def test_shape_rounds_with_size_and_gravity():
    earth, lunar = column(), moon()
    for c in (earth, lunar):
        state = c.state(0.0)
        ratios = [drops.axis_ratio(mm * 1e-3, c.gravity_m_s2, state['density_kg_m3'], 293.15)
                  for mm in (0.05, 1, 2, 4, 6)]
        assert ratios[0] > 0.9999
        assert np.all(np.diff(ratios) < 0)
    state = lunar.state(0.0)
    assert (drops.axis_ratio(6e-3, lunar.gravity_m_s2, state['density_kg_m3'], 293.15)
            > drops.axis_ratio(6e-3, earth.gravity_m_s2, earth.state(0.0)['density_kg_m3'], 293.15) + 0.2)


@pytest.mark.parametrize('c', [column(), moon()], ids=['earth', 'moon'])
def test_drips_stay_below_free_fall_and_terminal_speed(c):
    for height in (0.1, 1.0):
        r = drops.fall(3e-3, replace(c, cloud_base_m=height), initial_speed_m_s=0.0, sample_heights_m=[height, 0])
        assert r['survived']
        assert 0 < r['surface_speed_m_s'] < math.sqrt(2 * c.gravity_m_s2 * height)
        assert r['surface_speed_m_s'] < r['surface_terminal_speed_m_s']


def test_saturated_air_keeps_the_drop():
    """Saturated isothermal air neither takes water from a drop nor cools it; over 100 m the column's uniform vapour
    fraction leaves the air aloft 1.2% short of saturation, and the drop's small loss grows with the shortfall."""
    short = drops.fall(2e-3, column(humidity=1.0, lapse_k_m=0.0, cloud_base_m=0.1), initial_speed_m_s=0.0,
                       sample_heights_m=[0.1, 0.0])
    assert abs(short['evaporated_mass_fraction']) < 1e-8
    assert short['final_drop_temperature_k'] == pytest.approx(293.15, abs=1e-5)
    tall = drops.fall(2e-3, column(humidity=1.0, lapse_k_m=0.0, cloud_base_m=100.0), sample_heights_m=[100.0, 0.0])
    drier = drops.fall(2e-3, column(humidity=0.9, lapse_k_m=0.0, cloud_base_m=100.0), sample_heights_m=[100.0, 0.0])
    assert 0 < tall['evaporated_mass_fraction'] < 1e-3
    assert drier['evaporated_mass_fraction'] > 5 * tall['evaporated_mass_fraction']


def test_updraft_lengthens_the_fall_and_keeps_the_slip():
    """Ground speed drops by the updraft; the speed relative to the air, which sets drag, does not."""
    still = column(humidity=1.0, lapse_k_m=0.0, cloud_base_m=100.0)
    rising = replace(still, updraft_m_s=1.0)
    a = drops.fall(2e-3, still, sample_heights_m=[100.0, 0.0])
    b = drops.fall(2e-3, rising, sample_heights_m=[100.0, 0.0])
    assert b['fall_time_s'] > a['fall_time_s'] * 1.1
    assert b['surface_air_relative_speed_m_s'] == pytest.approx(a['surface_air_relative_speed_m_s'], rel=1e-3)
    assert a['surface_speed_m_s'] - b['surface_speed_m_s'] == pytest.approx(1.0, abs=1e-3)


def test_mass_closes_and_tighter_integration_agrees():
    central = moon(303.106, 0.765)
    for mm in (2.0, 5.0):
        ordinary = drops.fall(mm * 1e-3, central, sample_heights_m=[central.cloud_base_m, 0])
        tight = drops.fall(mm * 1e-3, central, sample_heights_m=[central.cloud_base_m, 0], rtol=2e-8)
        assert ordinary['status'] == tight['status'] == 'surface'
        assert abs(ordinary['evaporated_mass_fraction'] - tight['evaporated_mass_fraction']) < 1e-4
        kept = ordinary['surface_mass_kg'] / ordinary['initial_mass_kg']
        assert kept + ordinary['evaporated_mass_fraction'] == pytest.approx(1.0, abs=1e-12)
        track = ordinary['trajectory']
        assert track['altitude_m'][0] == central.cloud_base_m and track['altitude_m'][-1] == 0.0


def test_domain_is_enforced():
    with pytest.raises(ValueError, match='warm-water'):
        column(temperature=320.0)
    with pytest.raises(ValueError, match='below freezing'):
        column(temperature=275.0, humidity=0.5, cloud_base_m=5000.0)
    with pytest.raises(ValueError):
        drops.terminal_speed(0.05, column())
    with pytest.raises(ValueError):
        drops.liquid_saturation_pressure(270.0)
