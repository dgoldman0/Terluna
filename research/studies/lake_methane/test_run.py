"""The methane tracer: its budgets, its transport against the lake model's own, its kinetics against an analytic
profile, its bubbles against lake_anoxia's gas physics, and lake_anoxia's numbers where the physics coincide."""
import hashlib
import json
import math
from types import SimpleNamespace

import numpy as np
import pytest

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_anoxia import run as A
from research.studies.lake_methane import biology as B
from research.studies.lake_methane import columns as C
from research.studies.lake_methane import run as M
from research.studies.lake_methane import tracer as T
from research.studies.lake_oxygen import run as R
from shared.constants import MOON_SURFACE_GRAVITY

DT = 7150.0


def synthetic(depth=50.0, layer=1.0, steps=4000, every=12, t_c=20.0, oxygen=8.0, anoxic_bottom=True,
              sediment=0.25, deep_water=0.0, mixed=None, k_o2=1.0, debt=0.0):
    """A column the tests set up by hand: a pan of uniform temperature (so the model's diffusivity sits at its cap),
    oxygen everywhere but an anoxic bottom layer, and the floor's demand only in that layer."""
    basin = L.Basin.pan(depth, layer, 1e6)
    n = basin.n
    physics = L.Physics(gravity=MOON_SURFACE_GRAVITY, dissipation=1e-10)
    scen = [L.Scenario('test', gpp=0.0, attenuation=1.0, sediment=sediment, deep_water=deep_water)]
    snaps = steps // every + 1
    ox = np.full((snaps, n, 1), oxygen, np.float32)
    if anoxic_bottom:
        ox[:, -1, 0] = 0.0
    col = C.Column(name='synthetic', basin=basin, physics=physics, scenarios=scen, names=['test'], dt=DT,
                   snap_step=np.arange(snaps) * every, t=np.full((snaps, n), t_c, np.float32), oxygen=ox,
                   debt=np.full((snaps, n, 1), debt, np.float32), age=np.zeros((snaps, n), np.float32),
                   solute=np.zeros((snaps, n), np.float32), solute_saved=True,
                   mixed_depth=np.full(steps, layer if mixed is None else mixed), surface_t=np.full(steps, t_c),
                   pressure=np.full(steps, 1.2 * 101325.0), k_o2=np.full(steps, k_o2),
                   open_water=np.ones(steps, bool), ml_o2=None)
    return col


def first_order_lane(k1, **kw):
    """A lane whose methanotrophs hold a fixed stock and oxidise at first order k1 (per s): methane far below its
    half-saturation, oxygen far above its, the stock lost at once to its seed."""
    seed_um, kc_um, y = 1.0, 1e6, 0.3
    q = k1 * (kc_um * 1e-3) / (seed_um * 1e-3)                      # mol CH4 (mol C)-1 s-1
    return T.Lane(0, mu20_per_h=q * y * 3600.0, yield_c=y, k_ch4_um=kc_um, k_o2_um=1e-6, seed_c_um=seed_um,
                  loss_per_day=1e9, activation_ev=0.0, ebullition=False, **kw)


def budget_residual(bud, inv0, inv1, j):
    """Methane in: what the water and the floors put in, and the bubbles that dissolve; out: oxidised, anaerobic,
    to the air (the air's uptake back in)."""
    s = lambda k: bud[k][:, j].sum()
    change = inv1[j] - inv0[j]
    flows = (s('made_water') + s('floor_entered') + s('bubbles_dissolved') - s('oxidised') - s('aom_water')
             - s('evaded') + s('invaded'))
    return change - flows, max(abs(change), s('made_water') + s('floor_entered'))


# --- the tracer's physics ------------------------------------------------------------------------------------------

def test_steady_diffusion_against_first_order_oxidation_matches_the_analytic_profile():
    # K C'' = k1 C above an anoxic floor layer feeding a flux F, C = 0 at the surface: C(x) = F sinh(lambda (x + d))
    # / (K lambda cosh(lambda Le)), the share reaching the air 1 / cosh(lambda Le). The surface layer empties to the
    # air after each implicit diffusion step, the model's own order, which adds the resistance of d = dz c / V above
    # its centre (c the step's exchange across its floor, V its volume), and the anoxic floor layer oxidises nothing
    # over its upper half, so the oxidising column runs Le = L + d - dz / 2.
    col = synthetic(steps=40000)
    kdiff = col.physics.max_diffusivity + L.DIFFUSIVITY_O2
    dz = col.basin.dz[0]
    length = col.basin.centre[-1] - col.basin.centre[0]
    lam = 1.5 / length
    k1 = kdiff * lam ** 2
    tr = T.Tracer(col, [first_order_lane(k1)], air_ch4=0.0, oxygen_budget=False, mixed_layer_oxygen=False)
    state = tr.initial(0)
    out = tr.run(state, start=1, stop=36000, record=False)
    late = tr.run(state, start=36000, stop=40000, record=True)
    flux = tr.floor_mol[-1, 0] / col.basin.area[0]                                # mol m-2 s-1 into the column
    escaped = late['budget']['evaded'][:, 0].sum() / (4000 * DT) / col.basin.area[0]
    oxidised = late['budget']['oxidised'][:, 0].sum() / (4000 * DT) / col.basin.area[0]
    d = dz * (kdiff * col.basin.area[1] * DT / dz) / col.basin.volume[0]
    le = length + d - 0.5 * dz
    assert escaped / flux == pytest.approx(1.0 / math.cosh(lam * le), rel=2e-3)
    assert (escaped + oxidised) / flux == pytest.approx(1.0, rel=1e-6)                 # steady
    x = col.basin.centre - col.basin.centre[0]
    analytic = flux * np.sinh(lam * (x + d)) / (kdiff * lam * math.cosh(lam * le))
    c = state['c'][:, 0]
    assert np.allclose(c[1:-1], analytic[1:-1], rtol=5e-3)              # the surface layer is emptied each step
    assert out['budget']['made_water'][:, 0].sum() < 1e-9 * out['budget']['floor_entered'][:, 0].sum()   # the floor


def test_the_methane_budget_closes_with_every_road():
    # a shallow pan whose productive floor bubbles, with anaerobic oxidation in its mud and water and the air above
    col = synthetic(depth=12.0, layer=1.0, steps=3000, sediment=2.0, deep_water=0.5, mixed=3.0, k_o2=2e-5)
    col.oxygen[:, 6:, 0] = 0.0                                                     # anoxic below 6 m
    lanes = [T.Lane(0, **B.EARTH), T.Lane(0, aom_floor_mol_m2_s=1e-8, aom_water_per_s=1e-6,
                                          aom_water_cap_mol_m2_s=1e-9, **B.EARTH),
             T.Lane(0, aerobic=False, **B.EARTH)]
    probe = T.Tracer(col, lanes)
    bub = T.Bubbles(col.basin, col.t.astype(float).mean(0), col.physics.gravity, probe.floor_production_max())
    assert bub.active.any()
    tr = T.Tracer(col, lanes, bubbles=bub)
    state = tr.initial(0)
    inv0 = tr.inventory(state['c'])
    out = tr.run(state, start=1, record=False)
    inv1 = tr.inventory(state['c'])
    for j in range(len(lanes)):
        resid, scale = budget_residual(out['budget'], inv0, inv1, j)
        assert abs(resid) < 1e-9 * scale, j
    b = out['budget']
    assert b['bubbles'][:, 0].sum() > 0 and b['aom_floor'][:, 1].sum() > 0 and b['aom_water'][:, 1].sum() > 0
    assert b['oxidised'][:, 2].sum() == 0.0
    # the floor's methane all goes somewhere: anaerobic in the mud, bubbles, into the water or oxidised at its surface
    for j in range(len(lanes)):
        parts = b['aom_floor'][:, j].sum() + b['bubbles'][:, j].sum() + b['floor_entered'][:, j].sum() \
            + b['floor_surface_oxidised'][:, j].sum()
        assert parts == pytest.approx(b['floor_total'][:, j].sum(), rel=1e-12)
        assert b['bubbles_air'][:, j].sum() + b['bubbles_dissolved'][:, j].sum() == pytest.approx(
            b['bubbles'][:, j].sum(), rel=1e-9)


def test_the_oxygen_difference_keeps_its_own_budget():
    # the model reoxidises a debt holding methane; the lanes' methanotrophs take it or leave it to the air, and the
    # difference changes only by what the model spent, what the lane used and what went back to the air
    col = synthetic(depth=30.0, layer=1.0, steps=3000, oxygen=2.0, debt=4.0, mixed=4.0, k_o2=2e-5)
    lanes = [T.Lane(0, **B.EARTH), T.Lane(0, aerobic=False, **B.EARTH)]
    tr = T.Tracer(col, lanes, mixed_layer_oxygen=False)
    state = tr.initial(0)
    d0 = state['d'].T @ tr.v
    out = tr.run(state, start=1, record=False)
    d1 = state['d'].T @ tr.v
    b = out['budget']
    for j in range(2):
        change = d1[j] - d0[j]
        flows = b['o2_model_g'][:, j].sum() - b['o2_used_g'][:, j].sum() - b['o2_returned_to_air_g'][:, j].sum()
        assert change == pytest.approx(flows, rel=1e-9, abs=1e-6 * b['o2_model_g'][:, j].sum())
    assert b['o2_used_g'][:, 1].sum() == 0.0 and b['o2_model_g'][:, 1].sum() > 0.0
    assert d1[1] + b['o2_returned_to_air_g'][:, 1].sum() > d1[0] + b['o2_returned_to_air_g'][:, 0].sum()


def test_methanotrophs_grow_on_methane_and_fall_back_to_their_seed():
    col = synthetic(depth=10.0, layer=1.0, steps=2000, sediment=1.0, mixed=1.0)
    lane = T.Lane(0, ebullition=False, **B.EARTH)
    tr = T.Tracer(col, [lane], oxygen_budget=False)
    state = tr.initial(0)
    tr.run(state, start=1, stop=1000, record=False)
    grown = state['b'][:, 0].copy()
    assert grown[-2] > 3 * tr.seed[0]                        # just above the methane-making floor
    tr.floor_mol[:] = 0.0
    tr.water_demand[:] = 0.0
    state['c'][:] = 0.0
    tr.run(state, start=1000, stop=2000, record=False)
    assert np.all(state['b'][:, 0] < grown) and np.allclose(state['b'][:, 0], tr.seed[0], rtol=0.6)


def test_mixed_layers_recover_the_models_mixing():
    basin = L.Basin.from_hypsography(352.0, 5.8e8, np.linspace(0, 1, 41), np.linspace(1, 0, 41) ** 1.5)
    rng = np.random.default_rng(1)
    k_full = rng.integers(1, basin.n, 400)
    frac = np.where(rng.random(400) < 0.3, 0.0, rng.random(400))
    depth = np.where(frac > 0, basin.z[k_full] + frac * basin.dz[k_full], basin.z[k_full])
    k, f = T.mixed_layers(basin, depth.astype(np.float32))
    assert np.all(k == k_full)
    assert np.allclose(f, frac, atol=2e-4)


def test_floor_and_water_demands_give_back_the_models_demand():
    for site in (0, 1, 2):
        lake = json.loads(R.LAKES.read_text())['chosen'][site]
        h = lake['hypsography']
        for depth in (30.0, 352.0, 1925.0):
            basin = L.Basin.from_hypsography(depth, lake['area_km2'] * 1e6, h['depth_fraction'], h['area_fraction'],
                                             max_m=8.0 if depth <= 500.0 else 25.0)
            for s in R.scenarios():
                floor = T.floor_demand(basin, s) * basin.sediment
                whole = L.demand_profile(basin, s) * basin.volume
                water = whole - floor
                assert np.all(water > -1e-9 * whole.max())
                assert floor.sum() + water.sum() == pytest.approx(whole.sum(), rel=1e-12)


def sloping(depth, layer):
    """A basin whose area falls linearly to nothing at its floor, so that every layer has floor."""
    z = np.arange(0.0, depth + 1e-9, layer)
    return L.Basin(z, 1e6 * (1.0 - z / depth) + 1.0)


def test_the_bubbles_follow_lake_anoxias_gas_physics():
    basin = sloping(60.0, 2.0)
    t_mean = np.linspace(28.0, 22.0, basin.n)
    bub = T.Bubbles(basin, t_mean, MOON_SURFACE_GRAVITY, production_max=1e-6)
    rng = np.random.default_rng(2)
    col = SimpleNamespace(basin=basin)
    for _ in range(40):
        i = int(rng.integers(0, col.basin.n))
        prod = float(10 ** rng.uniform(-8.5, -6.5))
        bottom = float(rng.uniform(0.0, 0.3))
        p = np.zeros((col.basin.n, 1))
        p[i, 0] = prod
        b = np.full((col.basin.n, 1), bottom)
        mine = bub.shares(p, b)[i, 0]
        n2 = T.N2_FRACTION * (A.SURFACE_ATM - float(L.vapour_pressure_atm(t_mean[i])))
        theirs = G.ebullition_share(prod, col.basin.centre[i], t_mean[i], A.SURFACE_ATM, MOON_SURFACE_GRAVITY, n2,
                                    bottom_ch4=bottom)
        assert mine == pytest.approx(theirs, rel=2e-3, abs=1e-6)
    # each bubble's methane is kept: what reaches the air and what dissolves on the way make the whole
    active = np.flatnonzero(bub.active)
    assert active.size
    total = bub.survival[active] + bub.weights[active].sum(1)
    assert np.allclose(total, 1.0, atol=1e-9)
    s = bub.survival[active]
    assert np.all(np.diff(s) <= 1e-12)                         # a deeper start keeps less
    assert np.all(bub.weights[active] >= -1e-12)


def test_the_porewater_balance_gives_lake_soppens_ebullition():
    # Langenegger et al. (2019): at 10 m and 5 C, production near 8 mmol m-2 d-1 sends 10-17% up as bubbles; a
    # floor twice as deep, at the same production, sends fewer; at lunar gravity the floor at 10 m sends more
    surface = 0.944 / 1.01325
    basin = sloping(24.0, 0.5)
    p = np.zeros((basin.n, 1))
    i10 = int(np.argmin(np.abs(basin.centre - 10.25)))
    i20 = int(np.argmin(np.abs(basin.centre - 20.25)))
    p[[i10, i20], 0] = 8.0 / 86400e3
    earth = T.Bubbles(basin, np.full(basin.n, 5.0), 9.81, production_max=1e-6, surface_atm=surface,
                      n2_fraction=0.78)
    share = earth.shares(p, np.zeros((basin.n, 1)))[:, 0]
    assert 0.1 < share[i10] < 0.17
    assert share[i20] < share[i10]
    moon = T.Bubbles(basin, np.full(basin.n, 5.0), MOON_SURFACE_GRAVITY, production_max=1e-6, surface_atm=surface,
                     n2_fraction=0.78)
    assert moon.shares(p, np.zeros((basin.n, 1)))[i10, 0] > share[i10]


def test_methane_crosses_the_surface_at_the_schmidt_ratio_of_oxygens_transfer():
    for t_c in (5.0, 20.0, 35.0):
        sc = float(G.schmidt('CH4', t_c)) / float(L.schmidt_o2(t_c))
        assert sc == pytest.approx(T.SCHMIDT_RATIO, rel=1e-9)
    assert T.SCHMIDT_RATIO ** -0.5 == pytest.approx(0.928, abs=0.002)


def test_the_acceptors_supply_follows_the_sediment_and_the_soils():
    s = B.acceptor_supply()
    sed = s['sediment_kg_m2_yr']
    assert 0.1 < sed[0] < sed[1] < 3.0                           # kg m-2 a-1 over the lakes
    fe_h = s['iron']['highland']
    assert fe_h['reactive_wt_percent'] == pytest.approx(5.1 * 55.845 / 71.844 * 0.43, rel=1e-9)
    assert fe_h['ch4_mol_m2_yr'][0] == pytest.approx(sed[0] * fe_h['reactive_wt_percent'] * 10.0 / 55.845 / 8.0,
                                                     rel=1e-9)
    assert s['iron']['mare']['ch4_mol_m2_yr'][1] > fe_h['ch4_mol_m2_yr'][1]
    assert s['sulphate']['ch4_mol_m2_yr'] == pytest.approx(1.4e-3 * 1.8, rel=1e-9)
    assert s['manganese']['highland']['ch4_mol_m2_yr'][1] < 0.2 * fe_h['ch4_mol_m2_yr'][1]


# --- against the lake model's saved columns and lake_anoxia -----------------------------------------------------

def anoxia_columns_match_their_product():
    """lake_anoxia's saved columns where its committed product was made from them, else None."""
    meta_path = C.ANOXIA_RUNS / 'profiles.json'
    if not meta_path.is_file():
        return None
    product = json.loads(A.OUT.read_text())
    recorded = product['producer']['inputs'].get('research/runs/lake_anoxia/profiles.json')
    if recorded != C.sha256(meta_path):
        return None
    return product


@pytest.fixture(scope='module')
def resort_column():
    product = anoxia_columns_match_their_product()
    if product is None:
        pytest.skip("lake_anoxia's saved columns are not the ones its product was made from")
    st = R.Setting()
    meta = C.anoxia_meta(C.ANOXIA_RUNS)
    name = 'deep_resort_0-15'
    return st, C.anoxia_column(st, name, meta, folder=C.ANOXIA_RUNS), product


def test_the_tracer_makes_lake_anoxias_methane_where_the_rules_coincide(resort_column):
    st, col, product = resort_column
    year = int(round(T.YEAR / col.dt))
    lanes = [T.Lane(j, ebullition=False, **B.EARTH) for j in range(4)]
    # lake_anoxia's current rule for the water methane is made in, read from its daily columns as it reads them
    tr = T.Tracer(col, lanes, mixed_layer_oxygen=False, made_below=T.lake_anoxia_made_below(), oxygen_budget=False)
    state = tr.initial(0)
    out = tr.run(state, start=1, stop=year + 1, record=False)
    snaps = int(np.searchsorted(col.snap_step, year)) + 1
    meta = C.anoxia_meta(C.ANOXIA_RUNS)
    names = list(meta['scenarios'])
    scen = C.scenarios_for(names, col.version)
    z = np.load(C.ANOXIA_RUNS / f'{col.name}.npz')
    rec = dict(t=z['t'][:snaps], oxygen=z['oxygen'][:snaps])
    for j, trophic in enumerate(A.TRAITS):
        mine = (out['budget']['made_water'][:, j].sum() + out['budget']['made_floor'][:, j].sum()) \
            / col.basin.area[0] / (year * col.dt / T.YEAR)
        theirs, _ = A.anoxic_methane(SimpleNamespace(scenarios=scen), rec, col.basin,
                                     names.index(f'{trophic}/{C.CENTRAL}'), 0.5)
        assert mine == pytest.approx(theirs, rel=0.01), trophic
    # and the overturns lake_anoxia counted
    assert len(tr.events) == product['deep_lakes'][col.name]['overturn_events']


def test_the_tracers_transport_follows_the_models_own_age_and_debt(resort_column):
    st, col, product = resort_column
    year = int(round(T.YEAR / col.dt))
    lanes = [T.Lane(0, kind='age')] + [T.Lane(j, kind='debt') for j in range(4)]
    tr = T.Tracer(col, lanes, made_below=T.PRE_FIX_MADE_BELOW, oxygen_budget=False)
    state = tr.initial(0)
    out = tr.run(state, start=1, stop=year + 1, profiles=True)
    idx = np.searchsorted(col.snap_step, out['snap_steps'])
    v = col.basin.volume
    age_model = col.age[idx].astype(float)
    age_err = out['snap_c'][:, :, 0] - age_model
    assert np.abs(age_err).mean() < 0.3                                  # days, on ages of tens to hundreds of days
    # the lake's mean age to 3% of itself (an overturn's onset, entraining old water, carries the largest error)
    lake_mean = age_model @ v / v.sum()
    assert np.all(np.abs(age_err @ v / v.sum()) <= np.maximum(0.03 * lake_mean, 1.0))
    for j in range(4):
        model = (col.debt[idx][:, :, j].astype(float) * 0.5 / 64.0) @ v / col.basin.area[0]
        clone = out['snap_c'][:, :, 1 + j] @ v / col.basin.area[0]
        assert abs(clone.mean() - model.mean()) <= max(0.02 * model.mean(), 1e-3), j      # mol m-2


# --- the product ------------------------------------------------------------------------------------------------

@pytest.fixture(scope='module')
def product():
    if not M.OUT.is_file():
        pytest.skip('the product is written once the corrected columns are traced')
    return json.loads(M.OUT.read_text())


def digest(path):
    return hashlib.sha256((M.ROOT / path).read_bytes()).hexdigest()


def traced_files(product, path):
    """The saved files the product traced from a run folder's index, by case, or None for an index it does not
    list them for."""
    if path.endswith('sea_oxygen/columns.json'):
        seas = product['seas']
        return {n: c['files'] for part in (seas, seas.get('variants', {})) for n, c in part.get('columns', {}).items()}
    if path.endswith('lake_kinds/ladder.json'):
        return product['kinds']['files']
    return None


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == M.SCHEMA
    for path, h in product['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in product['producer']['inputs'].items():
        if path.startswith('research/runs/') and not (M.ROOT / path).is_file():
            continue                                   # the saved columns live on the research drive
        if digest(path) == h:
            continue
        # an index its producer re-keyed (its files renamed by hard links) must still name what was traced
        traced = traced_files(product, path)
        assert traced is not None, path
        now = json.loads((M.ROOT / path).read_text())['cases']
        for name, files in traced.items():
            ours = [files[k] for k in sorted(files)]
            assert [now[name]['files'][k] for k in sorted(now[name]['files'])] == ours, (path, name)
