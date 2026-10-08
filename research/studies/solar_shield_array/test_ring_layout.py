import json
from types import SimpleNamespace

import numpy as np
import pytest

from shared import constants as K
from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import ring_layout as rl


def flat_calibration(rate=0., e_middle=.11, e_edge=.14):
    """A calibration with a uniform turning and a forced eccentricity rising with the tilt, as the flights find."""
    radii, tilts = [17000., 23000.], [-30., 0., 30.]
    ecc = [[e_edge, e_middle, e_edge]]*2
    return rl.Calibration(radii, tilts, [[rate]*3]*2, ecc)


def test_the_two_sided_layout_keeps_its_clearance_pair_by_pair():
    cal = flat_calibration()
    heights = np.arange(-600., 600.+1e-9, 9.286)
    radii = rl.two_sided(cal, heights, 19000., 1.)
    tilts = np.degrees(np.arcsin(heights/radii))
    ecc = cal.at(radii, tilts)[1]
    check = rl.clear(radii, heights, ecc, 1.)
    assert check['holds'] and check['shingle_km'] >= 1.-1e-6 and check['nodes_km'] >= 1.-1e-6
    # The radius grows toward both edges, the sides alternating.
    middle = np.argmin(np.abs(heights))
    assert radii[middle] == radii.min()
    assert radii[0] > radii[middle] and radii[-1] > radii[middle]


def test_the_staircase_keeps_its_strips_clear_at_apolune_and_perilune():
    cal = flat_calibration()
    heights = np.arange(-600., 600.+1e-9, 9.286)
    radii = rl.staircase(cal, heights, 19000., 1., True)
    assert np.all(np.diff(radii[np.argsort(heights)]) > 0)
    ecc = cal.at(radii, np.degrees(np.arcsin(heights/radii)))[1]
    assert rl.clear(radii, heights, ecc, 1.)['holds']


def test_rings_at_one_radius_have_no_clearance():
    cal = flat_calibration()
    heights = np.arange(-100., 100.+1e-9, 9.286)
    radii = np.full(len(heights), 20000.)
    ecc = cal.at(radii, np.degrees(np.arcsin(heights/radii)))[1]
    assert not rl.clear(radii, heights, ecc, 1.)['holds']


def test_the_calibration_holds_its_lookups_to_the_flown_ranges():
    cal = rl.Calibration([18000., 22000.], [-26., 26.], [[-.01, -.01], [.01, .01]], [[.1, .1], [.1, .1]])
    inside, _ = cal.at(20000., 0.)
    beyond, _ = cal.at(30000., 0.)
    assert np.isclose(inside[0], 0.) and np.isclose(beyond[0], .01)
    assert not cal.covers(np.array([30000.]), np.array([0.]))


def test_the_turning_is_the_trend_of_the_pole_against_the_sun():
    # The Sun turns at 0.9856 degrees a day about a fixed orbit pole; a ring tilted 20 degrees toward it turns 0.003
    # degrees a day faster.
    sun_rate, extra, tilt = np.radians(.9856), np.radians(.003), np.radians(20.)
    days = np.arange(0., 60., 1/24)
    times = days*K.JULIAN_DAY

    def at(t):
        d = t/K.JULIAN_DAY
        sun = 1.5e11*np.array([np.cos(sun_rate*d), np.sin(sun_rate*d), 0.])
        return dict(positions=dict(sun=sun, earth=np.array([-3.84e8, 0., 0.])),
                    moon_v=np.array([0., 1e3, 0.]), earth_v=np.zeros(3))
    env = SimpleNamespace(at=at)
    phase = (sun_rate+extra)*days
    pole = np.c_[np.sin(tilt)*np.cos(phase), np.sin(tilt)*np.sin(phase), -np.cos(tilt)*np.ones_like(phase)]
    # Any position and velocity whose cross product points along the pole.
    a = np.c_[-np.sin(phase), np.cos(phase), np.zeros_like(phase)]
    r = 2e7*a
    v = 1e3*np.cross(pole, a)
    traj = np.concatenate([r, v], axis=1)[:, None, :]
    out = rl.turning(env, times, traj, every=1)[0]
    assert np.isclose(out['rate_deg_per_day'], .003, atol=1e-6)
    assert np.isclose(out['lean_deg'], 20., atol=1e-6)


def test_the_matched_layout_grows_toward_both_edges_and_keeps_its_clearance_when_held():
    # A turning that needs more radius the further a ring is tilted: the matched layout puts each ring near that radius.
    radii_grid, tilts = [17000., 23000.], [-30., 0., 30.]
    rate = [[-.05-.1*abs(t)/30. for t in tilts], [.25-.1*abs(t)/30. for t in tilts]]
    cal = rl.Calibration(radii_grid, tilts, rate, [[.14, .11, .14]]*2)
    heights = np.arange(-900., 900.+1e-9, 9.286)
    ideal = rl.ideal_radii(cal, heights)
    assert ideal[np.argmin(np.abs(heights))] < ideal[0] and ideal[np.argmin(np.abs(heights))] < ideal[-1]
    radii = rl.matched(cal, heights, ideal, 0., 1., 1., 'held')
    ecc = cal.at(radii, np.degrees(np.arcsin(heights/radii)))[1]
    assert rl.clear(radii, heights, ecc, 1., held=True)['holds']
    assert radii[np.argmin(np.abs(heights))] < radii[0] and radii[np.argmin(np.abs(heights))] < radii[-1]


def test_a_stack_on_steeply_rising_forced_eccentricities_runs_away():
    # Forced eccentricity rising fast with radius: each strip further out swings in at perilune, so the steps grow.
    cal = rl.Calibration([17000., 23000.], [-30., 0., 30.], [[0.]*3]*2, [[.10, .10, .10], [.40, .40, .40]])
    heights = np.arange(-400., 400.+1e-9, 9.286)
    assert rl.staircase(cal, heights, 19000., 1., True, 'forced', cap=21000.) is None
    assert rl.staircase(cal, heights, 19000., 1., True, 'held', cap=21000.) is not None


def test_held_eccentricities_keep_the_perilune_clearance_and_turning_follows_them():
    heights = np.array([-9.286, 0., 9.286])
    radii = np.array([20000., 19998., 20002.])
    forced = np.array([.13, .11, .13])
    kept = rl.held_profile(radii, heights, forced, 1.)
    peri = radii*(1-kept)/(1+kept)
    assert peri[2]-peri[1] >= 1.-1e-9 and peri[0]-peri[1] >= 1.-1e-9
    assert np.all(kept <= forced+1e-12)
    change = rl.held_turning_change(forced, kept)
    assert np.all(change <= 1e-12) and change.min() < 0      # less eccentricity, a slower turn


def test_the_flight_code_version_is_bumped_with_the_flight_code():
    # Saved flights are keyed to FLIGHT_CODE; change it whenever flights, turning or frozen_centre change.
    import hashlib
    import inspect
    sources = inspect.getsource(rl.flights)+inspect.getsource(rl.turning)+inspect.getsource(rl.frozen_centre)
    assert hashlib.sha256(sources.encode()).hexdigest()[:16] == '9b8cbc85b322d84c', 'bump ring_layout.FLIGHT_CODE'


@pytest.mark.skipif(not rl.OUT.exists(), reason='ring_layout.json not written')
def test_product_is_current_and_its_layouts_hold_their_clearance():
    product = json.loads(rl.OUT.read_text())
    assert product['schema'] == 'terluna.research.ring-layout/1' and not product['smoke']
    producer = product['producer']
    assert sorted(producer['files']) == sorted(rl.FILES)
    for name, value in producer['files'].items():
        assert rl.digest(rl.ROOT/name) == value, name
    assert producer['inputs']['results/plane_motion.json'] == rl.digest(rl.PLANES_OUT)
    assert producer['inputs']['results/frozen_rings.json'] == rl.digest(rl.FROZEN_OUT)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == rl.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])
    for clearance, layouts in product['layouts'].items():
        assert not layouts['one_radius']['clearance']['holds']
        for mode in ('forced', 'held'):
            for name, layout in layouts[mode].items():
                if 'runs_past_km' in layout:
                    continue
                assert layout['clearance']['holds'], (clearance, mode, name)
                # Layouts that fail may stack past the calibrated radii, which the product flags; the matched
                # layouts, the design, stay within them.
                assert layout['within_calibration'] or not (mode == 'held' and name.startswith('matched')), \
                    (clearance, mode, name)
        assert 'runs_past_km' not in layouts['held']['matched']
        coupled = layouts['held']['matched_coupled']
        assert coupled['within_calibration'] and coupled['clearance']['holds']
    assert product['flown']['matched']['steering_over_reach_max'] < 1
