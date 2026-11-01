"""The storm-water study: Green and Ampt's solution under steady rain at both gravities, the gravity scalings (the
conductivity and suction terms, the sorptivity they leave unchanged, the 36-fold time for gravity to take over, the
plan's first-hour uptake and the drops' kinetic energy against atmosphere/rain's solved fall speeds), the seal against
Rawls et al.'s crust factor, mass conservation on synthetic storms on every surface, the storm libraries and the
Poisson mixing, and the product's binding to its code and inputs. The product itself regenerates with run.py."""
import json
import math

import numpy as np
import pytest

from research.studies.storm_water import extract as EX
from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import run as R
from research.studies.storm_water import storms as ST

MOON = IN.MOON


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def ponded_time(p, f, dtheta):
    """Green and Ampt's time (h) to take in f mm when ponded from the start on a bare soil."""
    a = p['psi'] * dtheta
    return (f - a * math.log(1.0 + f / a)) / p['ks']


def steady_uptake(p, rate, hours, dtheta):
    """Green and Ampt under steady rain on a bare soil: ponding at F_p = K psi dtheta / (i - K), then the implicit
    solution, solved by bisection for F at the given time."""
    a = p['psi'] * dtheta
    fp = p['ks'] * a / (rate - p['ks'])
    tp = fp / rate
    if hours <= tp:
        return rate * hours
    lo, hi = fp, rate * hours
    for _ in range(200):
        f = 0.5 * (lo + hi)
        t = tp + (f - fp - a * math.log((f + a) / (fp + a))) / p['ks']
        lo, hi = (f, hi) if t < hours else (lo, f)
    return 0.5 * (lo + hi)


@pytest.mark.parametrize('gravity', [1.0, MOON])
@pytest.mark.parametrize('soil', ['sandy_loam', 'loam'])
def test_steady_rain_matches_green_ampt(gravity, soil):
    p = IN.parameters(soil, gravity=gravity, retention_mm=0.0)
    dth = np.array([0.35])
    for rate in (20.0, 60.0, 150.0):
        st = IN.State(dth.shape)
        for _ in range(120):
            IN.step(p, st, dth, rate, 1.0 / 60.0)
        exact = steady_uptake(p, rate, 2.0, 0.35)
        assert st.f[0] == pytest.approx(exact, rel=5e-3)
        assert st.f[0] + st.runoff[0] == pytest.approx(rate * 2.0, rel=1e-12)


def test_gravity_scalings():
    for name in IN.SURFACES:
        e, m = IN.parameters(name, gravity=1.0, pressure_pa=121590.0), IN.parameters(name, pressure_pa=121590.0)
        assert m['ks'] / e['ks'] == pytest.approx(MOON)                    # k rho g / mu
        assert m['psi'] / e['psi'] == pytest.approx(1.0 / MOON)            # capillary pressure over rho g
        assert m['r'] / e['r'] == pytest.approx(1.0 / MOON) if e['zc'] else m['r'] == 0.0
        assert m['ks'] * m['psi'] == pytest.approx(e['ks'] * e['psi'])     # so the sorptivity is unchanged
        sorp_e, sorp_m = math.sqrt(2 * e['ks'] * e['psi'] * 0.35), math.sqrt(2 * m['ks'] * m['psi'] * 0.35)
        assert sorp_m == pytest.approx(sorp_e)
        assert (sorp_m / m['ks']) ** 2 / (sorp_e / e['ks']) ** 2 == pytest.approx(MOON ** -2)
    assert MOON ** -2 == pytest.approx(36.4, abs=0.2)


def test_first_hour_uptake_as_the_plan_derives():
    """The plan: a dry sandy loam takes 30 mm in the first hour on the Moon against 33-36 on Earth, a loam 15 against
    16-17, at Rawls' parameters; dry here is a deficit of 0.35."""
    def first_hour(soil, gravity, dtheta):
        p = IN.parameters(soil, gravity=gravity)
        lo, hi = 0.0, 200.0
        for _ in range(100):
            f = 0.5 * (lo + hi)
            lo, hi = (f, hi) if ponded_time(p, f, dtheta) < 1.0 else (lo, f)
        return f
    assert first_hour('sandy_loam', MOON, 0.35) == pytest.approx(30.0, abs=1.0)
    assert 33.0 <= first_hour('sandy_loam', 1.0, 0.35) <= 37.0
    assert first_hour('loam', MOON, 0.35) == pytest.approx(15.0, abs=1.0)
    assert 16.0 <= first_hour('loam', 1.0, 0.35) <= 17.5


@pytest.mark.parametrize('gravity', [1.0, MOON])
def test_seal_reproduces_rawls_crust_factor(gravity):
    """With the front's suction set aside, the two-layer resistance gives Rawls et al.'s (1990) effective conductivity
    SC Ks / (1 + psi_i / L) at Earth's gravity (their eq. 5, with psi_i, SC and Zc for a loam), times g at any other."""
    p = IN.parameters('sealed_loam', gravity=gravity)
    ks, sc, psi_i, zc = IN.SOILS['loam']['ks_mm_h'], IN.SEAL['subcrust_factor'], IN.SEAL['matric_drop_mm'], IN.SEAL['thickness_mm']
    for depth in (20.0, 100.0, 500.0):
        ke = depth / (p['r'] + (depth - zc) / p['k_sub'])
        assert ke == pytest.approx(gravity * sc * ks / (1.0 + psi_i / depth), rel=1e-12)


def test_drop_energy_against_the_rain_solver():
    """The plan's v ~ sqrt(g / rho_air) at 1.2 atm gives 0.138 of Earth's energy per mm; atmosphere/rain's solved
    speeds in the resort's air give 0.11-0.17 at 2-4 mm drops, with 3 mm closest."""
    k = float(IN.drop_energy_factor(121590.0))
    assert k == pytest.approx(MOON / 1.2, rel=1e-3) and k == pytest.approx(0.138, abs=0.001)
    d = json.loads((R.ROOT / 'atmosphere/rain/results/rain_drops.json').read_text())
    f = d['fall_speed']['resort']
    assert d['airs']['resort']['moon']['pressure_pa'] == pytest.approx(121590.0)
    ratio = {dm: (vm / ve) ** 2 for dm, vm, ve in zip(d['diameters_mm'], f['moon']['terminal_speed_m_s'],
                                                       f['earth']['terminal_speed_m_s']) if ve}
    assert 0.10 < ratio[2.0] < k < ratio[4.0] < 0.18
    assert ratio[3.0] == pytest.approx(k, rel=0.05)


@pytest.mark.parametrize('name', list(IN.SURFACES))
@pytest.mark.parametrize('gravity', [1.0, MOON])
def test_mass_conserves_on_synthetic_storms(name, gravity):
    rng = np.random.default_rng(7)
    p = IN.parameters(name, gravity=gravity, pressure_pa=121590.0)
    dth = IN.deficit(p, np.array([0.0, 0.2, 0.6]))
    st = IN.State(dth.shape)
    rain = 0.0
    for storm in range(4):
        for _ in range(int(rng.integers(5, 90))):
            r = float(rng.exponential(30.0))
            IN.step(p, st, dth, r, 1.0 / 60.0)
            rain += r / 60.0
        IN.end_storm(p, st)
    assert np.allclose(st.f + st.runoff, rain, rtol=1e-12)
    assert (st.runoff >= 0).all() and (st.store == 0).all()


def test_forming_seal_lies_between_the_bare_and_sealed_loam():
    lib = dict(storms=[np.full(40, 60.0)], weights=[1.0], step_min=1.0)
    run = {}
    for label, p in (('bare', IN.parameters('loam')), ('sealed', IN.parameters('sealed_loam')),
                     ('moon_drops', IN.parameters('sealed_loam_forming', pressure_pa=121590.0)),
                     ('earth_drops', IN.parameters('sealed_loam_forming', energy_factor=1.0))):
        run[label] = float(IN.single_storms(p, lib, np.array([0.1]))[0][0, 0])
    assert run['bare'] < run['moon_drops'] < run['earth_drops'] < run['sealed']


def test_walnut_gulch_family_against_measured_storms():
    lib = ST.walnut_gulch()
    fit = lib['fit']
    for row in fit['osborn']:                       # the fitted counts and minutes above 50-200 mm/h
        assert 0.6 < row['events'] / row['measured_events'] < 1.5
        assert 0.6 < row['minutes'] / row['measured_minutes'] < 1.6
    for row in fit['mendez']:                       # held out: the annual maxima at 5-60 minutes
        assert 0.85 < row['two_year'] / row['measured_two_year'] < 1.3
        assert 0.85 < row['ten_year'] / row['measured_ten_year'] < 1.3
    d = ST.depths(lib)
    assert d.mean() == pytest.approx(lib['parameters']['mean_depth_mm'], rel=1e-3)


def test_storm_lives_reader():
    lib = ST.from_storm_lives(R.STORMS, R.ROOT)
    raw = json.loads((R.ROOT / R.STORMS).read_text())
    assert len(lib['storms']) == sum(1 for t in raw['tracks'] if t['whole'] and sum(t['series']['rain_mm_h']) > 0)
    t = next(t for t in raw['tracks'] if t['whole'])
    assert ST.depths(lib)[0] == pytest.approx(sum(t['series']['rain_mm_h']) * raw['output_s'] / 3600.0)
    with pytest.raises(ValueError):
        ST.from_storm_lives('atmosphere/rain/results/rain_drops.json', R.ROOT)


def test_poisson_mixing():
    lam = np.array([0.0, 0.3, 2.0, 12.0])
    s = np.full(4, 0.1)
    line = np.arange(31)[:, None] * np.ones((1, R.SATURATIONS.size)) * 2.5
    assert np.allclose(R.poisson_mix(line, lam, s), 2.5 * lam)          # beyond the table it extends linearly
    assert np.allclose(R.poisson_mix(np.ones((31, R.SATURATIONS.size)), lam, s), 1.0)


def test_bucket_reads_as_dry_topsoil():
    soil = IN.SOILS['sandy_loam']
    s = IN.initial_saturation(soil, np.array([0.0, 0.0125, 0.5]))
    residual = soil['porosity'] - soil['effective_porosity']
    assert s[0] == pytest.approx((soil['wilting_point'] - residual) / soil['effective_porosity'])
    assert s[2] == pytest.approx((soil['field_capacity'] - residual) / soil['effective_porosity'])
    assert IN.initial_saturation(soil, 0.0125, 'air_dry') == 0.0


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, sha in product['producer']['files'].items():
        assert R.sha256(path) == sha, path
    for path, sha in product['producer']['inputs'].items():
        assert R.sha256(path) == sha, path
    for key in ('evidence', 'reading_rule'):
        assert len(product[key]) > 200
    from shared.provenance import constants_changed
    assert not constants_changed(product['producer']['constants'])


@pytest.mark.skipif(not EX.OUT.is_file(), reason='the design-years extract is ignored and regenerated by extract.py')
def test_ignored_inputs_match(product):
    for path, sha in product['producer']['ignored_inputs'].items():
        assert EX.sha256(R.ROOT / path) == sha, path


def test_product_cells_and_classes(product):
    cells = product['cells']
    assert len(cells) == product['basis']['cells']
    assert {c['unep'] for c in cells} == set(R.CLASSES)
    assert product['basis']['land_shares']['semi-arid'] == pytest.approx(0.251, abs=0.002)
    for c in cells:
        for n in R.SURFACES:
            for b in R.BRACKETS:
                for m in R.MODES:
                    v = c['runoff_ratio'][n][b][m]['yearly']
                    assert v is None or 0.0 <= v <= 1.0
    for n in R.SURFACES:
        for b in R.BRACKETS:
            cl = product['classes'][n][b]
            for name in R.CLASSES:
                assert cl['independent']['classes'][name]['runoff_ratio'] <= cl['continuing']['classes'][name]['runoff_ratio'] + 1e-9
    assert product['walnut_gulch_check']['rule'].startswith('Set before the chain was run')


@pytest.mark.skipif(not EX.OUT.is_file(), reason='the design-years extract is ignored and regenerated by extract.py')
def test_composite_reproduces_the_classes():
    """The extract's tenth-of-day composites give back the aridity each cell was classed by, as design_years.py forms it."""
    x, _ = EX.load()
    by, mean = EX.composite(x)
    pet = EX.demand(mean['tas'], mean['ps'], mean['rss'] + mean['rls'])
    assert np.allclose(mean['pr'] / pet, x['aridity'], rtol=1e-5)
    assert all(lo <= a < hi for a, c in zip(x['aridity'], x['cls']) for name, lo, hi in EX.UNEP if name == c)
