"""The coast-water study: the belts' steady flow against its analytic solution and geography's screen, the night's
balance against Earth's measured dew and FAO-56's own equation, the 16 px/deg and 100 m belt rules on a synthetic
shore, the product's binding to its code and inputs, and the figures the README states. The full regeneration, about
seven minutes, runs with TERLUNA_COAST_WATER_REGENERATE=1 where the ignored inputs are present."""
import json
import math
import os
import re

import numpy as np
import pytest

from research.studies.coast_water import belts as B
from research.studies.coast_water import dew as D
from research.studies.coast_water import run as R
from research.studies.coast_water import shores as S

YEAR_S = R.YEAR_S
README = R.HERE / 'README.md'


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def ignored_present():
    return all((R.ROOT / p).is_file() for p in R.IGNORED)


# --------------------------------------------------------------------------------------------------------- belts
def test_belt_width_limits_and_geography_screen():
    """A steep shore's belt is reach over slope, a flat one sqrt(2TR/E): geography's 0.56-0.79 km at 5e-4 m2/s and a
    metre a year with 10-20 m to spare, 1.8-2.5 km at ten times the transmissivity."""
    e = 1.0 / YEAR_S
    assert B.width(0.1, 0.0, 15.0) == pytest.approx(150.0)
    assert B.width(0.0, e / (2 * 5e-4), 10.0) == pytest.approx(B.flat_width(5e-4, e, 10.0))
    assert [round(float(B.flat_width(5e-4, e, r)) / 1000, 2) for r in (10.0, 20.0)] == [0.56, 0.79]
    assert [round(float(B.flat_width(5e-3, e, r)) / 1000, 1) for r in (10.0, 20.0)] == [1.8, 2.5]
    # the island screen: a wooded island a kilometre across lowers its middle about 8 m at the central value
    assert float(B.water_table(500.0, 500.0, e / (2 * 5e-4))) == pytest.approx(7.9, abs=0.1)


def test_belt_edge_solves_the_steady_flow():
    """Finite differences of T h'' = E under the belt, h = 0 at the shore and no flow at the edge, put the water table
    c W^2 below the shore at the edge, and the edge on a steady slope where the ground and that fall reach the reach."""
    t, e, s, reach = 5e-4, 1.2 / YEAR_S, 0.004, 20.0
    c = e / (2 * t)
    w = float(B.width(s, c, reach))
    n = 2001
    x = np.linspace(0.0, w, n)
    dx = x[1] - x[0]
    a = np.zeros((n, n))
    rhs = np.full(n, e / t * dx * dx)
    a[0, 0], rhs[0] = 1.0, 0.0
    for i in range(1, n - 1):
        a[i, i - 1:i + 2] = (1.0, -2.0, 1.0)
    a[-1, -2:], rhs[-1] = (2.0, -2.0), e / t * dx * dx             # the mirror node: no flow at the edge
    h = np.linalg.solve(a, rhs)
    assert -h[-1] == pytest.approx(c * w * w, rel=1e-4)
    assert s * w - h[-1] == pytest.approx(reach, rel=1e-4)
    assert np.allclose(-h, B.water_table(x, w, c), rtol=1e-3, atol=1e-6)


def ramp(spacing_km, rise_per_km, n_land=40, n_sea=10, rows=60):
    """A straight north-south shore with the sea to the west and land rising eastward at a steady slope."""
    x = (np.arange(n_land + n_sea) - n_sea + 0.5) * spacing_km
    height = np.where(x > 0, x * rise_per_km, -10.0)
    return np.tile(height, (rows, 1)), x


@pytest.mark.parametrize('slope', [0.002, 0.01, 0.05])
def test_fine_rule_on_a_straight_shore(slope):
    """belts.py's rules on a Cartesian grid (run.fine_belts, the 100 m check's): on a straight shore along the grid,
    rising steadily, the belt per km of shore at 100 m spacing is the analytic width within a tenth (its first 100 m
    strip counts pi/4, the mean over directions), and at 1.9 km spacing the belt lies within the first pixel, where
    the strip is pi/4 of the analytic width on such a shore."""
    t, e, reach = 5.17e-4, 1.1 / YEAR_S, 15.0
    c = e / (2 * t)
    w = float(B.width(slope, c, reach))
    assert w < 1895.0
    for spacing, expected, tolerance in ((0.1, w, 0.1), (1.895, B.STAIRCASE * w, 0.02)):
        h, _ = ramp(spacing, slope * 1000.0, n_land=int(12 / spacing) + 5)
        area, _, _ = R.fine_belts(h, spacing * 1000.0, c, reach)
        rows = slice(5, -5)
        per_km = area[rows].sum() / (h[rows].shape[0] * spacing)
        assert per_km * 1000 == pytest.approx(expected, rel=tolerance), (spacing, per_km * 1000, expected)


# ---------------------------------------------------------------------------------------------------- the night
def test_earth_night_dew_within_measured_dew():
    """Ten clear hours under Earth's 87 W/m2 deficit at 90-96% humidity give the canopy 0.2-0.3 mm, within the
    measured 0.13 +/- 0.10 mm a night typical and the lysimeters' 0.5-0.6 mm most (Ritter et al. 2019)."""
    for rh in (0.90, 0.96):
        canopy, _ = D.earth_night(rh, 87.0, 10.0)
        assert 0.03 < float(canopy) < 0.6
    dry, _ = D.earth_night(0.70, 87.0, 10.0)
    assert float(dry) == 0.0


def test_canopy_balance_closes_and_dew_rises_with_humidity():
    t, p, u2 = np.full(4, 293.15), np.full(4, 121590.0), np.full(4, 0.6)
    rh = np.array([0.75, 0.90, 0.96, 0.99])
    e = rh * D.saturation_pressure(t)
    q = D.EPSILON * e / (p - (1 - D.EPSILON) * e)
    sky = D.STEFAN_BOLTZMANN * t ** 4 - 25.0
    dew, tc = D.canopy_dew(t, q, p, sky, u2)
    assert dew[0] == 0.0 and np.all(np.diff(dew) > 0)
    heat, vapour = D.leaf_conductances(t, p, u2, D.LEAF_WIDTHS_M[1], t - tc)
    ra = D.RA_REFERENCE / D.wind_2m(u2)
    rh_, rv = ra + 1 / (2 * D.CANOPY_LAI * heat), ra + 1 / (2 * D.CANOPY_LAI * vapour)
    rho = D.density(t, q, p)
    residual = ((1 - D.G_NIGHT) * D.LEAF_EMISSIVITY * (sky - D.STEFAN_BOLTZMANN * tc ** 4)
                + rho * D.CP * (t - tc) / rh_ + D.LV * np.maximum(rho * (q - D.saturation_humidity(tc, p)), 0) / rv)
    assert np.all(np.abs(residual) < 1e-6)


def test_demand_matches_fao56_daily_equation():
    """With the daily reference's 70 s/m and no soil heat, the general equation gives FAO-56's equation 6 at sea level
    within its constants' rounding."""
    t, p, u2, rn = 293.15, 101325.0, 2.0, 150.0
    rh = 0.6
    e = rh * D.saturation_pressure(t)
    q = D.EPSILON * e / (p - (1 - D.EPSILON) * e)
    sky = rn + D.STEFAN_BOLTZMANN * t ** 4
    ours = float(D.demand(np.array(0.0), np.array(sky), np.array(t), np.array(q), np.array(p), np.array(u2),
                          np.array(True), resistance=(70.0, 70.0), soil=(0.0, 0.0))) * 86400.0
    tc = t - 273.15
    es = 0.6108 * math.exp(17.27 * tc / (tc + 237.3))
    delta = 4098.0 * es / (tc + 237.3) ** 2
    gamma = 0.665e-3 * p / 1000.0
    fao = (0.408 * delta * rn * 0.0864 + gamma * 900.0 / (tc + 273.0) * u2 * (es - rh * es)) / \
        (delta + gamma * (1 + 0.34 * u2))
    assert ours == pytest.approx(fao, rel=0.03)


def test_dew_point_inverts_saturation():
    t, p = np.array([285.0, 300.0]), np.array([121590.0, 121590.0])
    q = D.saturation_humidity(t, p)
    assert np.allclose(D.dew_point(q, p), t, atol=1e-6)


# -------------------------------------------------------------------------------------------- the product's binding
def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, sha in product['producer']['files'].items():
        assert R.digest(path) == sha, path
    for path, sha in product['producer']['inputs'].items():
        assert R.digest(path) == sha, path
    for key in ('evidence', 'reading_rule'):
        assert len(product[key]) > 200


@pytest.mark.skipif(not ignored_present(), reason='the ignored inputs are on the research drive')
def test_ignored_inputs_match(product):
    for path, sha in product['producer']['ignored_inputs'].items():
        assert R.digest(path) == sha, path


def test_basis_is_climate_targets(product):
    targets = json.loads((R.ROOT / R.TARGETS).read_text())
    row = next(r for r in targets['resort_band']['rows'] if abs(r['ground_k'] - 301.0) < 1e-9)
    assert product['basis']['share_of_the_way'] == pytest.approx(row['share_of_the_way'], rel=1e-6)
    assert product['basis']['share_of_the_way'] == pytest.approx(
        (301.0 - 295.0230045582492) / (306.5234791232695 - 295.0230045582492), rel=1e-6)


def test_share_matched_masks_give_climate_targets_shares(product):
    d = product['deserts']['share_matched']
    for count in ('every_cell', 'model_land'):
        assert d[count]['share'] == pytest.approx(d[count]['climate_targets_share'], abs=0.001)


def test_water_table_grids_match_geography(product):
    table = json.loads((R.ROOT / R.TABLE).read_text())['cases']['central']['states']['mature']
    check = product['water_table']['check']
    assert check['land_km2'] == pytest.approx(table['land_km2'], rel=1e-4)
    assert check['median_depth_m'] == pytest.approx(table['median_depth_m'], abs=0.01)


@pytest.mark.skipif(os.environ.get('TERLUNA_COAST_WATER_REGENERATE') != '1' or not ignored_present(),
                    reason='regeneration takes about seven minutes; set TERLUNA_COAST_WATER_REGENERATE=1')
def test_regenerates(product):
    again = R.rounded(R.results())
    for key in ('basis', 'deserts', 'demand', 'belts', 'water_table', 'dew'):
        assert json.dumps(again[key], sort_keys=True) == json.dumps(product[key], sort_keys=True), key


# ------------------------------------------------------------------------------------------- the README's figures
NUMBER = re.compile(r'[−+-]?\d+(?:,\d{3})*(?:\.\d+)?')


def flat(text):
    return re.sub(r'\s+', ' ', text)


def numbers(text):
    return [float(m.replace('−', '-').replace(',', '')) for m in NUMBER.findall(text)]


def table(text, header):
    """The rows of the markdown table whose header row starts with a cell, keyed by their first cell."""
    rows, found = {}, False
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith('|'):
            if found:
                break
            continue
        cells = [c.strip() for c in line.strip('|').split('|')]
        if not found:
            found = cells[0] == header
        elif not set(cells[0]) <= set('-'):
            rows[cells[0]] = cells[1:]
    assert found and rows, header
    return rows


def same(written, value, decimals):
    assert round(written, decimals) == round(value, decimals), (written, value)


def test_readme_belts_table(product):
    rows = table(README.read_text(), 'Desert shore')
    d = product['belts']['deserts']['rain_mask.every_cell']
    lo, hi = d['central_15m']['by_class'], d['central_25m']['by_class']
    names = {'Sea coasts': 'sea_coast', 'Island coasts': 'island_coast',
             'Common-level lake shores': 'common_level_lake_shore', 'Open lake shores': 'open_lake_shore',
             'Closed lake shores': 'closed_lake_shore'}
    for label, cls in names.items():
        shore, width, per_km, belt, supply = (numbers(c) for c in rows[label])
        same(shore[0], round(lo[cls]['shore_km'], -2), 0)
        w15, w25 = lo[cls]['width_m'], hi[cls]['width_m']
        assert width == [round(w15['median']), round(w25['median']), round(w15['p10']), round(w25['p10']),
                         round(w15['p90']), round(w25['p90'])], (label, width)
        assert per_km == [round(1000 * v[cls]['belt_km2'] / v[cls]['shore_km']) for v in (lo, hi)], (label, per_km)
        assert belt == [round(lo[cls]['belt_km2'], -2), round(hi[cls]['belt_km2'], -2)], (label, belt)
        assert supply == [round(v[cls]['supply_m3_yr_per_km'] / 1e6, 2) for v in (lo, hi)]
    shore, _, _, belt, uptake = (numbers(c) for c in rows['All'])
    same(shore[0], round(d['central_15m']['shore_km'], -3), 0)
    assert belt == [round(d['central_15m']['belt_km2'], -2), round(d['central_25m']['belt_km2'], -2)]
    assert uptake == [round(d['central_15m']['uptake_km3_yr']), round(d['central_25m']['uptake_km3_yr'])]


def test_readme_slice_table(product):
    rows = table(README.read_text(), 'Slice')
    designs = product['slice_plan']['cost']['designs']
    names = {'2 km throughout': 'uniform_2_km', '500 m at the shore, stretched': 'stretched_500_m',
             '1 km throughout': 'uniform_1_km', '250 m at the shore, stretched': 'stretched_250_m'}
    for label, key in names.items():
        columns, step, hours = (numbers(c) for c in rows[label])
        v = designs[key]
        assert columns == [v['columns']] and step == [round(v['step_s'])], label
        assert hours == [round(h, 1) for h in v['wall_hours']], (label, hours, v['wall_hours'])


def found(pattern, text):
    m = re.search(pattern, text)
    assert m, pattern
    return [float(x.replace(',', '')) for x in m.groups()]


def test_readme_deserts(product):
    text = flat(README.read_text())
    d = product['deserts']
    pct = lambda name, count, n=1: round(100 * d[name][count]['share'], n)
    assert found(r"moved linearly in the ground's temperature leaves ([\d.]+)% of the land under 0\.5 mm a day "
                 r"counting every cell's land, and ([\d.]+)% on the model's own land cells, ([\d.]+) million km² at "
                 r"16 px/deg; "
                 r"moved by a constant factor for each kelvin it leaves ([\d.]+)% \(([\d.]+)%\), ([\d.]+) million km², "
                 r"close to climate_targets' ([\d.]+)% \(([\d.]+)%\)", text) == [
        pct('rain_mask', 'every_cell'), pct('rain_mask', 'model_land'),
        round(d['rain_mask']['every_cell']['land_km2_16ppd'] / 1e6, 2),
        pct('log_rain_mask', 'every_cell'), pct('log_rain_mask', 'model_land'),
        round(d['log_rain_mask']['every_cell']['land_km2_16ppd'] / 1e6, 2),
        round(100 * d['share_matched']['every_cell']['climate_targets_share'], 1),
        round(100 * d['share_matched']['model_land']['climate_targets_share'], 1)]
    against = d['against_interpolated_shares']
    lin, fac = against['additive']['bands'], against['log_linear']['bands']
    assert found(r'Of the ([\d.]+) points between the linear form and the interpolated shares, 30–45° carries '
                 r'([\d.]+): it is ([\d.]+)% desert in the design case and its rain doubles from ([\d.]+) to '
                 r'([\d.]+) mm a day', text) == [
        round(against['additive']['points'], 1), round(lin['30-45']['points'], 1),
        round(100 * lin['30-45']['design_desert'], 1), round(lin['30-45']['rain_mm_day']['design'], 2),
        round(lin['30-45']['rain_mm_day']['test'], 2)]
    assert found(r'leaves ([\d.]+)% of it desert at 301 K against ([\d.]+)% by the shares \(([\d.]+)% by the factor '
                 r'form\)\. 45–60° carries ([\d.]+): its rain rises from ([\d.]+) to ([\d.]+) mm a day, and it is '
                 r'([\d.]+)% desert against ([\d.]+)% \(([\d.]+)%\)\. 15–30° adds ([\d.]+) and 60–75° takes back '
                 r'([\d.]+)', text) == [
        round(100 * lin['30-45']['mask'], 1), round(100 * lin['30-45']['interpolated'], 1),
        round(100 * fac['30-45']['mask'], 1), round(lin['45-60']['points'], 1),
        round(lin['45-60']['rain_mm_day']['design'], 2), round(lin['45-60']['rain_mm_day']['test'], 2),
        round(100 * lin['45-60']['mask'], 1), round(100 * lin['45-60']['interpolated'], 1),
        round(100 * fac['45-60']['mask'], 1), round(lin['15-30']['points'], 1), round(-lin['60-75']['points'], 1)]
    k = d['points_per_kelvin']
    assert found(r'falling ([\d.]+) points a kelvin over 301–302 K where the factor form falls ([\d.]+) and the '
                 r'shares ([\d.]+)', text) == [round(-k['additive'], 1), round(-k['log_linear'], 1),
                                                round(-k['interpolated_shares'], 1)]
    assert found(r"the linear form's land under ([\d.]+) mm a day", text) == [
        round(d['share_matched']['every_cell']['threshold_mm_day'], 2)]


def test_readme_reach(product):
    text = flat(README.read_text())
    b = product['belts']['deserts']
    share = lambda region, case: 100 * b[region][case]['share_within_reach']
    cases = ('central_15m', 'central_25m', 'slow_15m', 'sandy_25m')
    assert found(r"The belts cover ([\d.]+)–([\d.]+)% of the linear mask's desert land at the central crust and "
                 r"([\d.]+)–([\d.]+)% across the crust's range, counting every cell's land; ([\d.]+)–([\d.]+)% "
                 r"\(([\d.]+)–([\d.]+)%\) on the model's own land cells; ([\d.]+)–([\d.]+)% \(([\d.]+)–([\d.]+)%\) on "
                 r"the factor mask; ([\d.]+)–([\d.]+)% \(([\d.]+)–([\d.]+)%\) on the share-matched deserts", text) == [
        round(share(region, c), 2) for region in ('rain_mask.every_cell', 'rain_mask.model_land',
                                                  'log_rain_mask.every_cell', 'share_matched.every_cell')
        for c in cases]
    check = product['belts']['resolution_check']
    ratios = check['ratio_range']
    lows, highs = [v[0] for v in ratios.values()], [v[1] for v in ratios.values()]
    central = ratios['central_15m'] + ratios['central_25m']
    assert found(r'gives ([\d.]+)–([\d.]+) times the belts resolved there \(([\d.]+)–([\d.]+) at the central crust\), '
                 r'above 1 in every case', text) == [round(min(lows), 2), round(max(highs), 2), round(min(central), 2),
                                                     round(max(central), 2)]
    assert min(lows) > 1.0
    # the corrected range, case by case: each case's share over its own largest and smallest ratio
    corrected = product['belts']['within_reach_corrected']
    for region in ('rain_mask.every_cell', 'log_rain_mask.every_cell'):
        by_case = {c: [share(region, c) / ratios[c][1], share(region, c) / ratios[c][0]] for c in ratios}
        assert [100 * v for v in corrected[region]['all_cases']] == pytest.approx(
            [min(v[0] for v in by_case.values()), max(v[1] for v in by_case.values())], rel=1e-5)
        mid = [v for c, v in by_case.items() if c.startswith('central')]
        assert [100 * v for v in corrected[region]['central']] == pytest.approx(
            [min(v[0] for v in mid), max(v[1] for v in mid)], rel=1e-5)
    lin, fac = corrected['rain_mask.every_cell'], corrected['log_rain_mask.every_cell']
    assert found(r'puts within reach \*\*([\d.]+)–([\d.]+)% of the desert land\*\* \(([\d.]+)–([\d.]+)% at the central '
                 r'crust\), and ([\d.]+)–([\d.]+)% on the factor mask', text) == [
        round(100 * lin['all_cases'][0], 1), round(100 * lin['all_cases'][1], 1), round(100 * lin['central'][0], 2),
        round(100 * lin['central'][1], 2), round(100 * fac['all_cases'][0], 2), round(100 * fac['all_cases'][1], 2)]
    halving = {n: [v['ratio_200m_over_100m'] for v in c['cases'].values()] for n, c in check['coasts'].items()}
    others = [x for n, v in halving.items() if n != 'w_procellarum' for x in v]
    assert found(r'whose 200 m belts are ([\d.]+)–([\d.]+) times its 100 m ones; on the other three coasts 200 m gives '
                 r'([\d.]+)–([\d.]+) times', text) == [round(min(halving['w_procellarum']), 2),
                                                       round(max(halving['w_procellarum']), 2), round(min(others), 2),
                                                       round(max(others), 2)]
    bands = b['rain_mask.every_cell']
    assert found(r'By band at the central crust: ([\d.]+)–([\d.]+)% at 30–45°, ([\d.]+)–([\d.]+)% at 45–60°, '
                 r'([\d.]+)–([\d.]+)% at 60–75° and ([\d.]+)–([\d.]+)% at 75–90°', text) == [
        round(100 * bands[c]['bands'][k]['share_within_reach'], 2 if k == '60-75' else 1)
        for k in ('30-45', '45-60', '60-75', '75-90') for c in ('central_15m', 'central_25m')]


def test_readme_night(product):
    text = flat(README.read_text())
    n = product['dew']['at_301_k']['central']
    classes = ('shore_land_in_water_cells', 'coastal_land_cells', 'inland_land_cells')
    assert found(r'the air at 2 m (\d+)% humid over the shore land of its water cells, (\d+)% on its coastal land '
                 r'cells and (\d+)% inland', text) == [round(100 * n[c]['desert']['rh']) for c in classes]
    assert found(r"water cells' air is about as humid as the land's \((\d+)% against (\d+)% and (\d+)%\)", text) == [
        round(100 * n[c]['desert']['rh']) for c in classes]
    onset = product['dew']['onset']
    assert found(r'gather dew above ([\d.]+)% humidity inland, ([\d.]+)% on the coastal land cells and ([\d.]+)% under',
                 text) == [round(100 * onset[c]['onset_rh'], 1) for c in ('inland_land_cells', 'coastal_land_cells',
                                                                          'shore_land_in_water_cells')]
    fog = product['dew']['desert_fog_at_301_k']
    assert found(r'its water is at most ([\d.]+) mm a lunar night with Earth\'s least deposition onto short '
                 r'vegetation, 2\.1 cm/s, and ([\d.]+) mm with its most, 8 cm/s, (\d+)% of the deserts\' rain',
                 text) == [
        round(fog['water_mm_lunar_night'][1], 2), round(fog['water_mm_lunar_night_at_earths_most'], 2),
        round(100 * fog['share_of_desert_rain_at_earths_most'])]
    reach = product['demand']['capillary_reach_m']['capillary_rise_m']
    assert found(r'sandy loam lifts the water from ([\d.]+) m and the finer soil from ([\d.]+) m', text) == [
        round(v, 1) for v in reach]
    dew = max(n[c]['desert']['canopy_mm_lunar_day'] for c in classes)
    rain = [n[c]['desert']['rain_mm_lunar_day'] for c in classes]
    cool = product['dew']['at_301_k']['cool']
    assert found(r'it gathers ([\d.]+) mm a lunar night against ([\d.]+)–([\d.]+) mm of rain a lunar day, and '
                 r'([\d.]+) mm with the air', text) == [round(dew, 3), round(min(rain), 1), round(max(rain), 1),
                                                        round(max(cool[c]['desert']['canopy_mm_lunar_day']
                                                                  for c in classes), 3)]
    rings = product['dew']['cm1']['rings']
    assert found(r'a canopy gathers ([\d.]+)–([\d.]+) mm a lunar night, ([\d.]+) mm on the equator\'s', text) == [
        round(rings['45n_wet']['canopy_mm_lunar_night'], 1), round(rings['45s_wet']['canopy_mm_lunar_night'], 1),
        round(rings['equator']['canopy_mm_lunar_night'], 1)]
    jura = product['dew']['cm1']['jura_box']
    assert found(r'night land breeze of ([\d.]+)–([\d.]+) m/s offshore', text) == sorted(
        round(-v['night_mean'], 2) for v in jura['onshore_wind_m_s'].values())
    earth = {(x['deficit_w_m2'], x['rh']): x['canopy_mm'] for x in product['dew']['earth_check']}
    assert found(r'the canopy gathers ([\d.]+) mm at 90% humidity and ([\d.]+) mm at 96%, between', text) == [
        round(earth[(87.0, 0.9)], 2), round(earth[(87.0, 0.96)], 2)]


def test_readme_closed_lakes(product):
    text = flat(README.read_text())
    c = product['belts']['closed_lakes']['rain_mask']
    lo = c['central_15m']
    assert found(r'The (\d+) closed desert lakes all lie at (\d+)–(\d+)°\. Round the median one the belt draws '
                 r'([\d.]+)–([\d.]+) times the lake\'s own evaporation at the central crust \(([\d.]+)–([\d.]+) across '
                 r'the range\)', text) == [
        lo['lakes'], math.floor(lo['lake_latitudes_deg'][0]), math.ceil(lo['lake_latitudes_deg'][1]),
        round(lo['draw_over_evaporation']['median'], 1), round(c['central_25m']['draw_over_evaporation']['median'], 1),
        round(min(v['draw_over_evaporation']['median'] for v in c.values()), 1),
        round(max(v['draw_over_evaporation']['median'] for v in c.values()), 1)]
    assert found(r'its row\'s sea cells\' at 301 K, ([\d.]+)–([\d.]+) mm a day under their marine air', text) == [
        round(v, 1) for v in lo['row_evaporation_mm_day']]
    assert found(r'the median lake is ([\d.]+) km², (\w+) pixels', text.replace('two pixels', '2 pixels')) == [
        round(lo['lake_km2']['median'], 2), lo['lake_pixels_median']]
    assert found(r'the belts draw ([\d.]+)–([\d.]+) km³ a year against the lakes\' ([\d.]+) km³', text) == [
        round(lo['draw_km3_yr'], 2), round(c['central_25m']['draw_km3_yr'], 2), round(lo['evaporation_km3_yr'], 1)]
    # the lakes' draw is the closed lake shores' own uptake, lake by lake
    d = product['belts']['deserts']['rain_mask.every_cell']['central_15m']['by_class']['closed_lake_shore']
    assert lo['draw_km3_yr'] == pytest.approx(d['uptake_km3_yr'], rel=1e-4)


def test_readme_islands(product):
    text = flat(README.read_text())
    isl = product['belts']['islands']
    lo, hi = isl['cases']['central_15m'], isl['cases']['central_25m']
    small = lambda r: r['by_size']['0-10_km2']
    assert isl['wooded_share'] == 0.95
    assert found(r'Of the Moon\'s ([\d,]+) islands, holding ([\d,]+) km² of land, those under 10 km² have (\d+)–(\d+)% '
                 r'of their land within reach at the central crust, and (\d+)–(\d+) of the ([\d,]+) lie wholly within '
                 r'it, their strips covering 95% of their land or more; islands of 10–100 km² (\d+)–(\d+)%', text) == [
        isl['count'], round(isl['land_km2'], -2), round(100 * small(lo['all'])['within_reach_share']),
        round(100 * small(hi['all'])['within_reach_share']), small(lo['all'])['wooded_count'],
        small(hi['all'])['wooded_count'], small(lo['all'])['count'],
        round(100 * lo['all']['by_size']['10-100_km2']['within_reach_share']),
        round(100 * hi['all']['by_size']['10-100_km2']['within_reach_share'])]
    s = isl['small']
    assert found(r'The small islands average ([\d.]+) km², one to (\w+) pixels', text.replace('to four', 'to 4')) == [
        round(s['mean_km2'], 1), round(s['pixels']['p90'])] and round(s['pixels']['p10']) == 1
    assert found(r'stay within the reach, (\d+)–(\d+) of the small islands are green throughout', text) == [
        small(lo['all'])['disc_rule_count'], small(hi['all'])['disc_rule_count']]
    assert found(r'In the deserts, (\d+)–(\d+) of (\d+) small islands are green throughout by the strips and '
                 r'(\d+)–(\d+) by the discs, and island land as a whole is ([\d.]+)–([\d.]+)% within reach', text) == [
        small(lo['desert'])['wooded_count'], small(hi['desert'])['wooded_count'], small(lo['desert'])['count'],
        small(lo['desert'])['disc_rule_count'], small(hi['desert'])['disc_rule_count'],
        round(100 * lo['desert']['within_reach_share'], 1), round(100 * hi['desert']['within_reach_share'], 1)]


def test_readme_water_table(product):
    text = flat(README.read_text())
    w = product['water_table']['rain_mask']['every_cell']
    assert found(r'stands at the ground on ([\d.]+)% of the rain mask\'s desert land and within 25 m of it on '
                 r'([\d.]+)%, a median (\d+) m down', text) == [round(100 * w['at_ground'], 2),
                                                                round(100 * w['within_25_m'], 1),
                                                                round(w['median_depth_m'])]
    assert found(r'only (\d+)% of it in cells beside a sea or lake', text) == [
        round(100 * w['within_25_m_beside_water']['one_cell'])]
    young = [w['young_1000_yr'], w['young_10000_yr']]
    assert found(r'the water table stands at the ground on ([\d.]+)% of the desert land and within 25 m on ([\d.]+)%, '
                 r'a median ([\d,]+) m down', text) == [
        round(100 * young[0]['at_ground'], 2), round(100 * young[0]['within_25_m'], 1),
        round(young[0]['median_depth_m'], -1)]
    assert round(100 * young[1]['at_ground'], 2) == round(100 * young[0]['at_ground'], 2)
    assert round(100 * young[1]['within_25_m'], 1) == round(100 * young[0]['within_25_m'], 1)


# ------------------------------------------------------------------------------------------------- the shores' geometry
def test_labels_merge_across_the_seam_and_the_poles():
    lab = np.zeros((8, 12), int)
    lab[2, 0], lab[2, 11] = 1, 2                    # one body cut by the seam
    lab[6, 0], lab[5, 11] = 3, 4                    # touching diagonally across the seam
    lab[0, 3], lab[0, 9] = 5, 6                     # across the north pole: half a turn apart in the polar row
    lab[7, 2], lab[7, 8] = 8, 9                     # and across the south pole
    lab[4, 5] = 7                                   # alone
    out = S.merge_labels(lab)
    assert out[2, 0] == out[2, 11] and out[6, 0] == out[5, 11] and out[0, 3] == out[0, 9] and out[7, 2] == out[7, 8]
    assert len(set(out[out > 0].tolist())) == 5 and out[0, 0] == 0


def test_banded_distances_against_brute_force():
    """The banded transforms give each pixel the distance to its nearest water pixel in its band's metric, the band's
    own east-west spacing, joined across the seam."""
    rng = np.random.default_rng(7)
    rows, cols, dy = 40, 64, 1.895
    lat = 60.0 - np.arange(rows) * 0.5
    water = rng.random((rows, cols)) < 0.02
    water[:, 0] = False
    water[5, 0] = True
    dist, near = S.distances(water, lat, dy, band_rows=8, margin_rows=rows)
    wr, wc = np.nonzero(water)
    for r in range(rows):
        band = slice((r // 8) * 8, min(rows, (r // 8) * 8 + 8))
        dx = dy * math.cos(math.radians(lat[band].mean()))
        dcol = np.abs(np.arange(cols)[:, None] - wc[None, :])
        dcol = np.minimum(dcol, cols - dcol)
        brute = np.sqrt(((r - wr)[None, :] * dy) ** 2 + (dcol * dx) ** 2).min(axis=1)
        assert np.allclose(dist[r], brute, rtol=1e-6, atol=1e-6), r
        nr, nc = np.divmod(near[r], cols)
        dn = np.minimum(np.abs(np.arange(cols) - nc), cols - np.abs(np.arange(cols) - nc))
        assert np.allclose(np.sqrt(((r - nr) * dy) ** 2 + (dn * dx) ** 2), brute, rtol=1e-6, atol=1e-6)
    # across the seam: the land in the last column, one column from the water in the first
    assert dist[5, cols - 1] == pytest.approx(dy * math.cos(math.radians(lat[0:8].mean())), rel=1e-6)


def test_shore_edges_lengths_and_depths():
    """A water pixel's four neighbours: those beside it to the east and west share an edge of dy and span their row's
    dx across it, those above and below share the row edge's length and span dy."""
    rows, cols, dy = 4, 6, 2.0
    lat = np.array([61.5, 60.5, 59.5, 58.5])
    dx = dy * np.cos(np.radians(lat))
    dx_edge = dy * np.cos(np.radians(np.array([62.0, 61.0, 60.0, 59.0, 58.0])))
    wc = np.zeros((rows, cols), np.int8)
    wc[1, 3] = S.SEA
    lv = np.where(wc > 0, 0.0, np.nan).astype(np.float32)
    h = np.where(wc > 0, -10.0, 5.0).astype(np.float32)
    keep = wc == 0
    cidx, edge, across, rise, lake = S.shore_edges(keep, wc, lv, np.zeros((rows, cols), np.int32), h, dy, dx, dx_edge)
    got = {tuple(divmod(int(i), cols)): (edge[k, 0], across[k, 0], rise[k, 0]) for k, i in enumerate(cidx)}
    assert set(got) == {(1, 2), (1, 4), (0, 3), (2, 3)}
    assert got[(1, 2)] == pytest.approx((dy, dx[1], 5.0)) and got[(1, 4)] == pytest.approx((dy, dx[1], 5.0))
    assert got[(0, 3)] == pytest.approx((dx_edge[1], dy, 5.0)) and got[(2, 3)] == pytest.approx((dx_edge[2], dy, 5.0))
    assert np.all(edge[:, 1:] == 0) and np.all(lake == 0)


@pytest.mark.parametrize('transmissivity, slope', [(5.17e-4, 0.002), (5.17e-3, 0.0005)])
def test_pipeline_on_a_straight_shore(transmissivity, slope):
    """The 16 px/deg pipeline on a synthetic strip of land at the equator between two straight shores along the grid,
    rising steadily from each: banded distances, the shore's edges and depths and run.belts_for give each shore pi/4
    of the analytic width in its first pixel, the staircase factor on a shore of one direction, and the analytic
    width beyond it in the pixels inland."""
    rows, cols, ppd = 48, 96, S.PPD
    lat = (rows / 2 - 0.5 - np.arange(rows)) / ppd
    dy = 2.0 * math.pi * R.MOON_RADIUS / 1000.0 / 360.0 / ppd
    dx = dy * np.cos(np.radians(lat))
    dx_edge = dy * np.cos(np.radians(rows / 2 / ppd - np.arange(rows + 1) / ppd))
    row_area = dx * dy
    col = np.arange(cols)
    sea = (col >= 30) & (col < 50)
    from_shore = np.minimum((col - 49) % cols, (30 - col) % cols)        # columns from the nearest water column
    height = np.where(sea, -10.0, slope * (from_shore - 0.5)[None, :] * dx[:, None] * 1000.0)
    water = np.broadcast_to(sea, (rows, cols)).copy()
    wc = np.where(water, S.SEA, 0).astype(np.int8)
    lv = np.where(water, 0.0, np.nan).astype(np.float32)
    dist, near = S.distances(water, lat, dy)
    keep = ~water & (dist <= S.NEAR_KM)
    cidx, edge, across, rise, lake = S.shore_edges(keep, wc, lv, np.zeros((rows, cols), np.int32),
                                                   height.astype(np.float32), dy, dx, dx_edge)
    coastal = np.zeros(rows * cols, bool)
    coastal[cidx] = True
    inner = np.flatnonzero(keep.ravel() & ~coastal)
    nr, nc = np.divmod(near.ravel()[inner].astype(np.int64), cols)
    ir, ic = np.divmod(inner, cols)
    step_c = (nc - ic) % cols
    step_c = np.where(step_c > cols // 2, step_c - cols, step_c)
    zero = lambda n: np.zeros(n, np.int64)
    px = dict(coast=dict(cell=zero(cidx.size), area=row_area[cidx // cols], edge=edge.astype(float),
                         across=across.astype(float), rise=rise.astype(float), island=np.zeros(cidx.size, bool)),
              inland=dict(cell=zero(inner.size), area=row_area[ir], distance=dist.ravel()[inner].astype(float),
                          step_rows=nr - ir, step_cols=step_c, height=height.ravel()[inner],
                          water=np.full(inner.size, S.SEA), island=np.zeros(inner.size, bool)))
    e, reach = 1.1 / YEAR_S, 15.0
    r = R.belts_for(px, np.full(2048, e), transmissivity, reach)
    per_km = (r['coast_areas'].sum() + r['inland_area'].sum()) / (2 * rows * dy)
    w = float(B.width(slope, e / (2 * transmissivity), reach)) / 1000.0
    depth = dx.mean()
    assert w < 2 * depth
    assert per_km == pytest.approx(B.STAIRCASE * min(w, depth) + max(w - depth, 0.0), rel=0.01)


def test_resolution_check_sequence(product):
    """The local coasts' belts at 400, 200 and 100 m, the 16 px/deg pixels' with the same uptake."""
    for name, coast in product['belts']['resolution_check']['coasts'].items():
        for case, v in coast['cases'].items():
            assert v['belt_km2_100m'] > 0 and v['belt_km2_200m'] > 0 and v['belt_km2_400m'] > 0, (name, case)
            assert v['ratio_16ppd_over_100m'] == pytest.approx(v['belt_km2_16ppd'] / v['belt_km2_100m'], rel=1e-5)
            assert v['ratio_200m_over_100m'] == pytest.approx(v['belt_km2_200m'] / v['belt_km2_100m'], rel=1e-5)
