"""Tests of the ring fleet's glow screen (glow.py).

    nice -n 10 python -m pytest research/studies/solar_shield_array/test_glow.py -q    # about three CPU minutes

The product test rebuilds results/glow.json from its inputs and needs the protection domain's optical inputs, the WHI
spectrum, the middle atmosphere's columns and the correlated-k tables; without them it is skipped and says why.
"""
import json
import math
import re

import numpy as np
import pytest

from shared import constants as K
from . import glow as g

MD = g.HERE/'glow.md'


def inputs_missing():
    missing = [str(p.relative_to(g.ROOT)) for p in g.OPTICAL if not p.exists()]
    try:
        from atmosphere.middle_atmosphere import fetch_inputs
        whi = fetch_inputs.path('whi2008_ref_solar_irradiance_ver2.dat')
        if not whi.exists():
            missing.append(str(whi))
        for name in g.COLUMNS:
            for path in (g.MIDDLE/'profiles'/f'{name}.csv', g.MIDDLE/'cases'/f'{name}.json'):
                if not path.exists():
                    missing.append(str(path.relative_to(g.ROOT)))
        from atmosphere.radiative_convective import ck
        if not (ck.CACHE/ck.table_key()).is_dir():
            missing.append('correlated-k tables')
    except Exception as error:                       # the line inputs the table key hashes are absent
        missing.append(f'correlated-k inputs ({error.__class__.__name__})')
    return missing


MISSING = inputs_missing()
needs_inputs = pytest.mark.skipif(bool(MISSING), reason=f'inputs not restored: {MISSING[:3]}')


@pytest.fixture(scope='module')
def product():
    return json.loads(g.OUT.read_text())


@pytest.fixture(scope='module')
def heat():
    return json.loads(g.HEAT.read_text())


@pytest.fixture(scope='module')
def layout():
    return json.loads(g.LAYOUT.read_text())


# ---------------------------------------------------------------------------------------------- optics and stack

def test_oblique_recursion_is_the_protection_domains_at_normal_incidence():
    from protection import model
    w = np.linspace(.25, 20., 400)
    layers = [np.full(w.shape, 1.46+1e-3j), np.full(w.shape, 2.4+.05j), np.full(w.shape, 1.2+0j)]
    thick = [.3, .1, 4.]
    r0, t0 = model.stack_rt(w, layers, thick)
    r1, t1 = g.stack_rt(w, layers, thick, 1.)
    assert np.max(np.abs(r0-r1)) < 1e-12 and np.max(np.abs(t0-t1)) < 1e-12
    # Lossless layers conserve energy at every incidence, both polarizations averaged.
    lossless = [np.full_like(w, 1.46), np.full_like(w, 2.4), np.full_like(w, 1.46)]
    for c in (1., .6, .2, .05):
        r, t = g.stack_rt(w, lossless, thick, c)
        assert np.max(np.abs(r+t-1)) < 1e-12


def test_blackbody_bands_close_and_match_quadrature():
    from scipy.integrate import quad
    edges = g.thermal_bins()
    for T in (90., 175., 300., 480.):
        total = g.band_radiance(np.array([T]), edges).sum()
        assert total == pytest.approx(K.STEFAN_BOLTZMANN*T**4/np.pi, rel=1e-12)
    for u in (.2, .5, 1.3, 4., 12.):
        exact = 15/np.pi**4*quad(lambda x: x**3/np.expm1(x), 0, u)[0]
        assert g.blackbody_fraction(np.array([u]))[0] == pytest.approx(exact, rel=1e-10, abs=1e-15)


def test_grey_layers_reproduce_the_provisioning_two_layer_formula():
    from provisioning import screen
    rng = np.random.default_rng(12)
    for _ in range(20):
        e_out, e_mo, e_mi = rng.uniform(.03, .9, 3)
        share, (outer, inner) = screen.two_layers(e_out, e_mo, e_mi)
        plates = g.Plates.stack([g.grey_plate(e_out, e_mi), g.grey_plate(e_out, e_mo)])
        res = g.solve(plates, [40., 40.])
        assert res['inward']/80. == pytest.approx(share, rel=1e-9)
        assert K.STEFAN_BOLTZMANN*res['T'][1]**4/40. == pytest.approx(outer, rel=1e-8)
        assert K.STEFAN_BOLTZMANN*res['T'][0]**4/40. == pytest.approx(inner, rel=1e-8)
        assert abs(res['residual']) < 1e-10*80.


def test_two_layer_estimate_is_rebuilt_from_the_same_inputs(heat):
    from provisioning import screen
    films, window = heat['tiles']['films'], heat['fleet']['window_ring_share']
    expected = {k: ({kk: vv for kk, vv in v.items() if kk != 'one_layer_share'} if isinstance(v, dict) else v)
                for k, v in screen.heat_mirror(films, window).items()}
    assert g.two_layer_check(heat) == expected
    stored = json.loads(g.FIRST.read_text())['heat']['heat_mirror']
    assert g.two_layer_check(heat) == {k: ({kk: vv for kk, vv in v.items() if kk != 'one_layer_share'}
                                           if isinstance(v, dict) else v) for k, v in stored.items()}


def test_first_screen_glow_is_rebuilt_from_the_same_inputs(heat, layout):
    from provisioning import screen
    rebuilt = g.first_screen(heat, layout)
    assert [round(x, 2) for x in rebuilt] == screen.heat(heat, layout)['fleet_infrared_w_m2'] == [5.66, 9.92]
    assert [round(x, 2) for x in rebuilt] == json.loads(g.FIRST.read_text())['heat']['fleet_infrared_w_m2']


@needs_inputs
def test_spectral_stacks_conserve_energy_and_coatings_stay_opaque():
    th = g.thermal_plate(.1, 4.)
    plate = dict(rho_out=th['rho_out'], rho_in=th['rho_in'], tau=th['tau'], edges_um=th['edges_um'])
    for n in (1, 3, 9):
        res = g.solve(g.Plates.repeat(plate, n), [50./n]*n)
        assert abs(res['residual']) < 1e-11*50.
        assert abs(res['F']) < 1e-9
    mirror = g.coated(plate, inner=.05, outer=.85)
    assert np.all(mirror['tau'] == 0) and np.allclose(mirror['rho_in'], .95) and np.allclose(mirror['rho_out'], .15)
    res = g.solve(g.Plates.stack([mirror, plate, plate]), [20., 20., 20.])
    assert abs(res['residual']) < 1e-11*60.
    # The film's transmittance is the same from both faces (reciprocity).
    assert th['tau_asymmetry'] < 1e-12


# --------------------------------------------------------------------------------------------------- geometry

def test_a_moon_facing_mirror_returns_the_moon_onto_itself_as_discs_overlap():
    R = K.MOON_RADIUS/1e3
    for r in (17000., 20000.):
        for delta in (0., 1.5, 3., 4., 6.):
            d = math.radians(delta)
            P = np.array([[[0., r, 0.]]])
            n = np.array([[[math.sin(d), math.cos(d), 0.]]])
            geo = g.specular_geometry({'P': P}, n)
            landed = (geo['ir_day']+geo['ir_night'])[0, 0]/geo['ir_in'][0, 0]
            x = d/math.asin(R/r)
            overlap = (2/math.pi)*(math.acos(x)-x*math.sqrt(1-x*x)) if x < 1 else 0.
            assert landed == pytest.approx(overlap, abs=.02)
            # The tile receives the Moon's light as a disc of radiance: pi (R/r)^2 per unit radiance, by cos(tilt).
            assert geo['ir_in'][0, 0] == pytest.approx(np.pi*(R/r)**2*math.cos(d), rel=1e-3)


def test_glow_lands_by_day_and_night_as_the_sub_point_lies():
    cpsi, share = g.landing_share()
    assert share[-1] == pytest.approx(1.) and share[0] == pytest.approx(0., abs=1e-12)
    assert share[cpsi.size//2] == pytest.approx(.5, abs=1e-9)
    assert np.all(np.diff(share) >= -1e-12)


def test_the_layout_radius_is_the_sunward_crossings(layout, heat):
    ledger = json.loads(g.LEDGER.read_text())
    rings = g.fleet_rings(layout, heat, ledger)
    S = g.fleet_samples(rings, 360, ring_index=[0, 929, 1858])
    sunward = np.argmin(np.abs(S['phase']), axis=1)
    r_sun = S['r'][np.arange(3), sunward]
    assert r_sun == pytest.approx(rings['r_apo'][[0, 929, 1858]], rel=1e-4)
    # The rings cross the aperture at the stack's heights.
    z = S['P'][np.arange(3), sunward, 2]
    assert z == pytest.approx(rings['heights'][[0, 929, 1858]], abs=1.)


# ----------------------------------------------------------------------------------------------------- product

@needs_inputs
def test_product_regenerates_exactly_from_inputs_with_matching_hashes(product):
    for name, value in product['inputs'].items():
        assert g.digest(g.ROOT/name) == value, name
    for name, value in product['producer']['files'].items():
        assert g.digest(g.ROOT/name) == value, name
    from atmosphere.radiative_convective import ck
    assert product['producer']['correlated_k_table'] == ck.table_key()
    assert all(product['checks'].values()), product['checks']
    rebuilt = g.build_product(g.compute(verbose=False))
    assert json.dumps(rebuilt, indent=1)+'\n' == g.OUT.read_text()


def test_the_screens_own_sums_close(product):
    r = product['results']
    for key, row in r['answers'].items():
        for k in ('global', 'day', 'night'):
            ir = row['glow_W_m2'][k]+row['infrared_return_W_m2'][k]
            assert row['infrared_total_W_m2'][k] == pytest.approx(ir, rel=1e-5), key
            assert row['heat_total_W_m2'][k] == pytest.approx(ir+row['sunlight_return_W_m2'][k], rel=1e-5), key
        # The global mean is the mean of the two hemispheres.
        for part in ('glow_W_m2', 'infrared_return_W_m2', 'sunlight_return_W_m2'):
            assert row[part]['global'] == pytest.approx(.5*(row[part]['day']+row[part]['night']), rel=1e-4), key
        assert row['absorbed_PW'] == pytest.approx(sum(row['absorbed_by_half_PW'].values()), rel=1e-5)
        assert row['absorbed_PW'] == pytest.approx(sum(row['absorbed_by_ring_PW'].values()), rel=1e-5)
        d = row['dimming']
        a = g.DESIGN['planetary_albedo']
        expect = 100*(row['infrared_total_W_m2']['global']/((1-a)*r['design_sunlight_W_m2'])
                      + row['sunlight_return_W_m2']['global']/r['design_sunlight_W_m2'])
        assert d['percent'] == pytest.approx(expect, rel=1e-5)
        # The product keeps six significant figures, so the difference closes to the larger term's rounding.
        assert d['to_add_percent'] == pytest.approx(d['percent']-d['coating_provides_percent'],
                                                    abs=1e-5*max(1., abs(d['percent'])))


def md_table(title):
    """The rows of the markdown table under a heading in glow.md, as lists of cells."""
    text = MD.read_text()
    block = text.split(title, 1)[1]
    rows = []
    for line in block.splitlines()[1:]:
        if line.startswith('|'):
            rows.append([c.strip() for c in line.strip('|').split('|')])
        elif rows:
            break
    return rows[2:]


def number(cell):
    return float(re.sub(r'[^0-9.\-]', '', cell.replace('−', '-')))


def test_glow_md_states_the_products_numbers_and_sums(product):
    answers = product['results']['answers']
    rows = md_table('## The answers screened')
    assert len(rows) >= 20
    for cells in rows:
        key = re.search(r'`([a-z0-9_]+)`', cells[0]).group(1)
        a = answers[key]
        glow, ir, light, total, dim = (number(c) for c in cells[1:6])
        assert glow == pytest.approx(a['glow_W_m2']['global'], abs=.051), key
        assert ir == pytest.approx(a['infrared_return_W_m2']['global'], abs=.051), key
        assert light == pytest.approx(a['sunlight_return_W_m2']['global'], abs=.051), key
        assert total == pytest.approx(a['heat_total_W_m2']['global'], abs=.051), key
        assert total == pytest.approx(glow+ir+light, abs=.151), key
        assert dim == pytest.approx(a['dimming']['percent'], abs=.051), key
    night = md_table('## By lunar day and night')
    for cells in night:
        key = re.search(r'`([a-z0-9_]+)`', cells[0]).group(1)
        a = answers[key]
        values = [number(c) for c in cells[1:7]]
        assert values[0] == pytest.approx(a['glow_W_m2']['day'], abs=.051), key
        assert values[1] == pytest.approx(a['glow_W_m2']['night'], abs=.051), key
        assert values[2] == pytest.approx(a['infrared_return_W_m2']['day'], abs=.051), key
        assert values[3] == pytest.approx(a['infrared_return_W_m2']['night'], abs=.051), key
        assert values[4] == pytest.approx(a['heat_total_W_m2']['night'], abs=.051), key
        assert values[5] == pytest.approx(a['night_lux_from_returns']['mean'], abs=.051), key
