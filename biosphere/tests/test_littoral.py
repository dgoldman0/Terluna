"""The littoral community model against what can be checked: its budgets of carbon, phosphorus, nitrogen and oxygen
closing through every exchange it names, the stage-structured maturation rate against its analytic limits, the thermal
response, the bubble point at lunar and Earth gravity, and the direction of each mechanism (light and darkness,
phosphorus, idling through the night, heat, hypoxia, clipping, refuge, predation, decay without oxygen, harvest)."""
import math

import numpy as np
import pytest

from biosphere.ecology import littoral as M
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

DT = 1.0 / 12.0                      # days
LUNAR_DAY = 29.53059
FULL = [(0.0, 'phyto', 0.05), (0.0, 'peri', 0.5), (0.0, 'chara', 5.0), (0.0, 'zoo', 0.02), (0.0, 'snail', 0.5),
        (0.0, 'cray_j', 0.2), (0.0, 'cray_a', 0.5), (0.0, 'fish_j', 0.1), (0.0, 'fish_a', 0.3)]


def run(community=None, lunar_days=2.0, schedule=FULL, setting=None, forcing=None, n=None, **forcing_kw):
    community = community or M.example_community()
    setting = setting or M.Setting(depth=np.array([1.0, 2.0, 4.0]), exchange=0.3, tp_open=8.0)
    n = n or np.size(setting.depth)
    model = M.Littoral(community, setting, n)
    state = model.initial()
    steps = int(round(lunar_days * LUNAR_DAY / DT))
    forcing = forcing or M.lunar_forcing(steps, DT, **forcing_kw)
    start = model.totals({k: v.copy() for k, v in state.items()})
    out = M.drive(model, state, forcing, DT, steps, schedule, record_every=6)
    return model, state, out, start


def budgets(model, state, out, start):
    """The change of each total against the accounts that should explain it: (carbon, phosphorus, nitrogen, oxygen),
    each as (change, explained)."""
    a, seeds = out['accounts'], out['seeds']
    c1, p1, n1, o1 = model.totals(state)
    c0, p0, n0, o0 = start
    carbon = (c1 - c0, seeds['c'] + a['gpp_c'] - a['resp_c'] - a['burial_c'] - a['emerged_c'] - a['harvest_c']
              - a['ch4_c'] + a['c_exchange'])
    phosphorus = (p1 - p0, seeds['p'] - a['p_calcite'] - a['burial_p'] - a['emerged_p'] - a['harvest_p']
                  + a['p_exchange'])
    nitrogen = (n1 - n0, seeds['n'] + a['n_fixed'] - a['denitrified_n'] - a['burial_n'] - a['emerged_n']
                - a['harvest_n'] + a['n_exchange'])
    oxygen = (o1 - o0, a['o2_bio'] + a['o2_air'] - a['o2_bubbles'] + a['o2_exchange'])
    return carbon, phosphorus, nitrogen, oxygen


# ---- budgets ---------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('exchange', [0.0, 0.3])
def test_carbon_phosphorus_nitrogen_and_oxygen_budgets_close(exchange):
    setting = M.Setting(depth=np.array([1.0, 2.0, 4.0]), exchange=exchange, tp_open=8.0, harvest=0.5,
                        harvest_from_year=0.0)
    model, state, out, start = run(setting=setting, lunar_days=3.0)
    scale = dict(zip(('c', 'p', 'n', 'o'), (10.0, 100.0, 1000.0, 10.0)))
    for (change, explained), key in zip(budgets(model, state, out, start), 'cpno'):
        assert np.allclose(change, explained, rtol=0.0, atol=1e-9 * scale[key]), key
    assert np.all(out['accounts']['harvest_c'] > 0)
    for pool, series in out['series'].items():
        if pool != 'oxygen':
            assert np.all(series >= 0.0), pool


def test_budgets_close_through_anoxia_heat_and_starvation():
    """A dark, hot, loaded box: every pool drains, the oxygen goes negative as debt, and nothing is created."""
    setting = M.Setting(depth=np.array([0.5, 2.0]), exchange=0.0, tp_open=30.0)
    steps = int(round(1.5 * LUNAR_DAY / DT))
    f = M.lunar_forcing(steps, DT, noon_par=0.0, t_mean=37.0, t_swing=4.0, k=1e-7)
    model, state, out, start = run(setting=setting, lunar_days=1.5, forcing=f,
                                   schedule=FULL + [(0.0, 'det_c', 80.0)])
    carbon, phosphorus, nitrogen, oxygen = budgets(model, state, out, start)
    # the detritus seeded here carries no nutrients; its carbon is in the seeds' account
    assert np.allclose(carbon[0], carbon[1], atol=1e-8)
    assert np.allclose(phosphorus[0], phosphorus[1], atol=1e-7)
    assert np.allclose(nitrogen[0], nitrogen[1], atol=1e-6)
    assert np.allclose(oxygen[0], oxygen[1], atol=1e-8)
    assert np.all(out['o2_low'].min(0) == 0.0) and np.all(state['oxygen'] < 0.0)
    assert np.all(out['accounts']['ch4_c'] > 0.0)
    assert np.all(state['cray_a'] < 1e-3) and np.all(state['fish_a'] < 0.3)


# ---- the pieces --------------------------------------------------------------------------------------------------------
def test_thermal_response_rises_then_falls_and_kills_above_the_limit():
    th = M.Thermal(e_a=0.65, t_half_c=34.0, e_h=3.0, lethal_c=37.0, lethal_rate=1.0)
    assert float(M.performance(20.0, th)) == pytest.approx(1.0, abs=0.006)      # the reference, far below the optimum
    assert float(M.performance(34.0, th)) == pytest.approx(0.5 * float(M.arrhenius(34.0, 0.65)), rel=1e-12)
    t = np.linspace(5.0, 45.0, 401)
    perf = M.performance(t, th)
    peak = t[np.argmax(perf)]
    assert 29.0 < peak < 34.0 and perf[-1] < 0.1 * perf.max()
    # Q10 of the rising limb between 10 and 20 C for 0.65 eV
    assert float(M.arrhenius(20.0, 0.65) / M.arrhenius(10.0, 0.65)) == pytest.approx(2.48, abs=0.02)
    assert float(M.heat_mortality(36.9, th)) == 0.0 and float(M.heat_mortality(39.0, th)) == pytest.approx(2.0)


def test_maturation_rate_limits():
    z = 0.01
    assert float(M.maturation(0.05, 0.0, z)) == pytest.approx(0.05 / (1.0 - z))          # no mortality
    assert float(M.maturation(0.05, 0.05, z)) == pytest.approx(0.05 / -math.log(z))       # growth equal to mortality
    assert float(M.maturation(0.05, 0.05 * (1 + 1e-4), z)) == pytest.approx(0.05 / -math.log(z), rel=1e-3)
    assert float(M.maturation(0.0, 0.01, z)) == 0.0 and float(M.maturation(-0.02, 0.01, z)) == 0.0
    nu = np.linspace(0.001, 0.2, 50)
    g = M.maturation(nu, 0.02, z)
    assert np.all(np.diff(g) > 0) and np.all(g > 0)
    assert float(M.maturation(1e-4, 5.0, z)) < 1e-12                                      # crushed by mortality
    # the share of a cohort that survives to maturation at constant rates is z^(mu/nu); at steady state the
    # maturation flux over the births, gamma / (gamma + mu - nu), equals that share times the mass gained, z^(mu/nu - 1)
    for nu_, mu_ in ((0.05, 0.02), (0.02, 0.05), (0.1, 0.001)):
        gamma = float(M.maturation(nu_, mu_, z))
        assert gamma / (gamma + mu_ - nu_) == pytest.approx(z ** (mu_ / nu_ - 1.0), rel=1e-9)


def test_bubble_point_at_lunar_and_earth_gravity():
    moon = M.Littoral(M.Community(), M.Setting(depth=2.0, exchange=0.0, tp_open=5.0), 1)
    earth = M.Littoral(M.Community(), M.Setting(depth=2.0, exchange=0.0, tp_open=5.0, gravity=STANDARD_GRAVITY,
                                                pressure_pa=101325.0, o2_fraction=0.20946), 1)
    assert moon.bubble_ratio[0] == pytest.approx(1.0 + 998.0 * MOON_SURFACE_GRAVITY * 2.0 / (0.1746 * 121590.0))
    assert 1.14 < moon.bubble_ratio[0] < 1.16 and 1.90 < earth.bubble_ratio[0] < 1.94
    # sinking slows with gravity
    assert moon.sink_phyto / earth.sink_phyto == pytest.approx(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY)


def producers_only(**meadow):
    return M.Community(meadow=M.Meadow(**meadow))


SEED_PLANTS = [(0.0, 'phyto', 0.05), (0.0, 'peri', 0.5), (0.0, 'chara', 5.0)]


def test_light_feeds_oxygen_and_darkness_draws_it_down():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.0, tp_open=8.0)
    model, state, out, _ = run(producers_only(), lunar_days=2.0, schedule=SEED_PLANTS, setting=setting, k=3e-6)
    day = out['floor_light'][:, 0] > 50.0
    night = out['floor_light'][:, 0] < 1e-9
    sat = float(M.L.oxygen_saturation(29.5, 1.2, 0.142, 0.1746))
    assert out['o2'][day, 0].mean() > sat and out['o2'][night, 0].mean() < sat
    assert out['accounts']['o2_bubbles'][0] >= 0.0
    # in the dark nothing is fixed
    steps = int(round(LUNAR_DAY / DT))
    dark = M.lunar_forcing(steps, DT, noon_par=0.0)
    _, _, out_dark, _ = run(producers_only(), lunar_days=1.0, schedule=SEED_PLANTS, setting=setting, forcing=dark)
    assert out_dark['accounts']['gpp_c'][0] == 0.0 and out_dark['accounts']['calcite_meadow'][0] == 0.0


def test_supersaturated_water_bubbles_at_the_floor_pressure():
    """A meadow in bright light with almost no exchange to the air: the oxygen stops at the bubble point."""
    setting = M.Setting(depth=np.array([2.0]), exchange=0.0, tp_open=20.0)
    model, state, out, _ = run(producers_only(), lunar_days=0.3, schedule=[(0.0, 'chara', 10.0)], setting=setting,
                               k=1e-8, t_swing=0.0)
    sat = float(M.L.oxygen_saturation(29.5, 121590.0 / 101325.0, 0.142, 0.1746))
    assert out['accounts']['o2_bubbles'][0] > 0.0
    assert out['o2'].max() == pytest.approx(sat * model.bubble_ratio[0], rel=1e-6)


def test_without_phosphorus_in_the_water_growth_is_only_recycling():
    """A closed box with no dissolved phosphorus: the producers hold at most the phosphorus their seeds brought."""
    setting = M.Setting(depth=np.array([2.0]), exchange=0.0, tp_open=0.0)
    idling = producers_only(idle_share=0.25)          # so that the night's carbon does not decide the comparison
    model, state, out, _ = run(idling, lunar_days=3.0, schedule=SEED_PLANTS, setting=setting)
    s = out['series']
    living_p = sum(s[k][:, 0] * 1000.0 / model.ratios(k)[0] for k in ('phyto', 'peri', 'chara'))
    assert living_p.max() <= out['seeds']['p'][0] * (1 + 1e-12)
    # with phosphorus in the water the same seeds grow a larger meadow
    rich = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    _, _, grown, _ = run(idling, lunar_days=3.0, schedule=SEED_PLANTS, setting=rich)
    assert grown['series']['chara'][-1, 0] > 1.2 * s['chara'][-1, 0]


def test_idling_carries_the_meadow_through_the_night():
    """With upkeep idled in the dark the reserve lasts the night; without it the structure is burned."""
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    kept, burned = [], []
    for idle, store in ((0.25, kept), (1.0, burned)):
        com = M.Community(meadow=M.Meadow(idle_share=idle, upkeep=0.03))
        _, state, out, _ = run(com, lunar_days=4.0, schedule=[(0.0, 'chara', 5.0)], setting=setting)
        store.append(out['series']['chara'][:, 0])
        store.append(out['series']['reserve'][:, 0])
    assert kept[0][-1] > 2.0 * burned[0][-1]
    assert burned[1].min() == 0.0 and kept[1][len(kept[1]) // 2:].min() > 0.0


def test_heat_beyond_the_limit_kills_and_hypoxia_kills():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    _, mild, _, _ = run(lunar_days=0.5, setting=setting, t_mean=29.0, t_swing=1.0)
    _, hot, _, _ = run(lunar_days=0.5, setting=setting, t_mean=36.8, t_swing=0.4)
    assert hot['cray_a'][0] < 0.05 * mild['cray_a'][0]          # the crayfish's limit is 36 C
    assert hot['fish_a'][0] > 0.2 * mild['fish_a'][0]           # the fish's is 38 C: stressed, alive
    # a box closed to the air and the open water, in the dark over a loaded floor, runs out of oxygen
    closed = M.Setting(depth=np.array([2.0]), exchange=0.0, tp_open=8.0)
    steps = int(round(LUNAR_DAY / DT))
    f = M.lunar_forcing(steps, DT, noon_par=0.0, k=1e-8)
    _, anoxic, out, _ = run(lunar_days=1.0, setting=closed, forcing=f, schedule=FULL + [(0.0, 'det_c', 60.0)])
    assert out['o2_low'].min() == 0.0
    f = M.lunar_forcing(steps, DT, noon_par=0.0, k=3e-5)
    _, aired, _, _ = run(lunar_days=1.0, setting=closed, forcing=f, schedule=FULL + [(0.0, 'det_c', 60.0)])
    assert anoxic['fish_a'][0] < 0.2 * aired['fish_a'][0]


def test_crayfish_clip_the_meadow_and_stir_the_floor():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    plants = SEED_PLANTS + [(0.0, 'snail', 0.3)]
    _, without, _, _ = run(lunar_days=2.0, setting=setting, schedule=plants)
    model, with_cray, out, _ = run(lunar_days=2.0, setting=setting, schedule=plants + [(0.0, 'cray_a', 3.0)])
    assert with_cray['chara'][0] < 0.8 * without['chara'][0]
    assert with_cray['fines'][0] > 1.5 * without['fines'][0]
    a = out['accounts']
    assert a['clipped_chara_c'][0] == pytest.approx(1.5 * a['eaten_chara_c'][0], rel=1e-9)


def test_the_meadow_shelters_juveniles_and_predators_eat_the_rest():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    base = [(0.0, 'peri', 2.0), (0.0, 'det_c', 20.0), (0.0, 'cray_j', 1.0)]
    _, alone, _, _ = run(lunar_days=0.5, setting=setting, schedule=base)
    _, bare, _, _ = run(lunar_days=0.5, setting=setting, schedule=base + [(0.0, 'fish_a', 1.0)])
    _, meadow, _, _ = run(lunar_days=0.5, setting=setting, schedule=base + [(0.0, 'fish_a', 1.0), (0.0, 'chara', 60.0)])
    total = lambda s: s['cray_j'][0] + s['cray_a'][0]
    assert total(bare) < 0.8 * total(alone)
    assert total(meadow) > 1.1 * total(bare)


def test_adults_breed_only_as_far_as_the_cue_allows():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    feast = [(0.0, 'peri', 5.0), (0.0, 'det_c', 60.0), (0.0, 'snail', 3.0), (0.0, 'cray_a', 0.5)]

    def births(cue):
        com = M.example_community()
        cons = list(com.consumers)
        cons[2] = M.Consumer(**{**cons[2].__dict__, 'cue': cue})
        com.consumers = tuple(cons)
        _, _, out, _ = run(com, lunar_days=0.5, setting=setting, schedule=feast)
        return out['accounts']['births_c'][0]

    full, half, none = births(1.0), births(0.5), births(0.0)
    assert full > 0.0 and none == 0.0
    assert half == pytest.approx(0.5 * full, rel=0.05)


def test_drive_repeats_short_forcing_and_seeds_on_schedule():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    model = M.Littoral(M.example_community(), setting, 1)
    state = model.initial()
    one_day = M.lunar_forcing(int(round(LUNAR_DAY / DT)), DT)
    steps = int(round(2.2 * LUNAR_DAY / DT))
    out = M.drive(model, state, one_day, DT, steps, [(0.0, 'chara', 5.0), (0.1, 'snail', 0.4)], record_every=6)
    snail = out['series']['snail'][:, 0]
    first = np.argmax(snail > 0)
    assert out['time_day'][first] == pytest.approx(0.1 * 365.25, abs=6 * DT + 1e-9)
    assert out['seeds']['c'][0] == pytest.approx(5.0 * 1.1 + 0.4)
    assert out['time_day'][-1] == pytest.approx((steps // 6) * 6 * DT)


def test_halving_the_step_changes_a_lunar_day_little():
    setting = M.Setting(depth=np.array([2.0]), exchange=0.3, tp_open=8.0)
    ends = []
    for dt in (1.0 / 12.0, 1.0 / 24.0):
        model = M.Littoral(M.example_community(), setting, 1)
        state = model.initial()
        steps = int(round(LUNAR_DAY / dt))
        M.drive(model, state, M.lunar_forcing(steps, dt), dt, steps, FULL, record_every=steps)
        ends.append({k: float(v[0]) for k, v in state.items()})
    for pool in ('chara', 'peri', 'det_c', 'snail', 'cray_a', 'fish_a', 'oxygen'):
        assert ends[0][pool] == pytest.approx(ends[1][pool], rel=0.03, abs=2e-3), pool
