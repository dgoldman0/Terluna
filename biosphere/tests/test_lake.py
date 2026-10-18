"""The lake model against values and solutions that can be checked: oxygen at saturation, conservation of heat and
oxygen, an analytic diffusion mode, convective and wind deepening into linear stratification, reaeration and ice
growth by conduction; the bulk exchange paired as PlaSim pairs it, the density of respiration's anaerobic products, the
labile pool's decay through the night, the diffusivity's cap and the Arrhenius reference; and the deep mixing from the
wind's energy (its mass balance, its agreement with a fixed dissipation, its dimensions, the cap that no longer holds
it and its closed budgets) and from the Lake number."""
import math

import numpy as np
import pytest

from biosphere.ecology import lake as L

MOON_O2 = 21227.0 / 121590.0


def air():
    return L.Air(gas_constant=288.38, cp=1008.9)


def steady(steps, dt=3575.0, sw=0.0, lw=400.0, air_k=293.15, q=0.01, wind=3.0, p=121590.0):
    full = lambda x: np.full(steps, float(x))
    return L.Forcing(dt=dt, sw_down=full(sw), lw_down=full(lw), air_k=full(air_k), q_air=full(q), wind10=full(wind),
                     pressure=full(p), o2_fraction=MOON_O2)


def test_saturation_matches_standard_tables():
    # USGS DOTABLES (Benson & Krause): fresh water under Earth's air at one atmosphere
    assert L.oxygen_saturation(20.0) == pytest.approx(9.09, abs=0.005)
    assert L.oxygen_saturation(35.0) == pytest.approx(6.95, abs=0.005)
    assert L.oxygen_saturation(0.0) == pytest.approx(14.62, abs=0.01)
    for t in (5.0, 24.2, 34.6):
        assert L._saturation(t, 1.2, 0.14, MOON_O2) == pytest.approx(float(L.oxygen_saturation(t, 1.2, 0.14, MOON_O2)),
                                                                    rel=1e-12)


def test_lunar_air_holds_earths_oxygen_pressure_less_the_vapour():
    # 17.46% O2 at 1.2 atm is Earth's dry partial pressure; the vapour takes a smaller share of the moist air
    for t in (20.0, 35.0):
        u = float(L.vapour_pressure_atm(t))
        th = 0.000975 - 1.426e-5 * t + 6.436e-8 * t * t
        expected = MOON_O2 * (1.2 - u) * (1 - th * 1.2) / (L.O2_REFERENCE_FRACTION * (1 - u) * (1 - th))
        assert L.oxygen_saturation(t, 1.2, 0.0, MOON_O2) / L.oxygen_saturation(t) == pytest.approx(expected, rel=1e-12)
        assert 1.0 < expected < 1.02
    # salinity lowers it: at 35 C Earth's seawater holds a sixth less than fresh water, as Weiss's (1970) independent
    # salinity factor gives (USGS memorandum 2011.03, equation 3)
    tk = 35.0 + 273.15
    weiss = math.exp(35.0 * (-0.033096 + 0.014259 * tk / 100 - 0.0017 * (tk / 100) ** 2))
    assert L.oxygen_saturation(35.0, salinity=35.0) / L.oxygen_saturation(35.0) == pytest.approx(weiss, abs=0.002)


def test_basin_geometry():
    b = L.Basin.from_hypsography(352.0, 5e8, [0.0, 0.3, 1.0], [1.0, 0.5, 0.0])
    assert b.z[-1] == pytest.approx(352.0)
    assert b.sediment.sum() == pytest.approx(b.area[0])
    assert np.all(b.volume > 0)
    pan = L.Basin.pan(3.0, 0.25, 1e6)
    assert pan.n == 12 and pan.volume.sum() == pytest.approx(3e6)
    assert pan.sediment[:-1].sum() == 0.0 and pan.sediment[-1] == pytest.approx(1e6)


def test_demand_integrates_to_the_scenario():
    b = L.Basin.from_hypsography(352.0, 5e8, [0.0, 0.3, 1.0], [1.0, 0.5, 0.0])
    s = L.Scenario('x', gpp=5.3, attenuation=1.7, sediment=0.37, deep_water=0.53)
    total = float((L.demand_profile(b, s) * b.volume).sum())
    z_eu = 4.6 / 1.7
    expected = 0.8 * 5.3 * b.area[0] + 0.53 * float(np.interp(z_eu, b.z, b.area)) + 0.37 * b.area[0]
    assert total == pytest.approx(expected, rel=1e-9)
    shallow = L.Basin.pan(3.0, 0.25, 1e6)
    oligo = L.Scenario('o', gpp=0.53, attenuation=0.28, sediment=0.25, deep_water=0.07)
    # in a lake shallower than the euphotic zone the deep demand falls on the floor
    assert float((L.demand_profile(shallow, oligo) * shallow.volume).sum()) == pytest.approx(
        (0.8 * 0.53 + 0.07 + 0.25) * 1e6, rel=1e-9)


def test_heat_and_oxygen_are_conserved_through_freezing_and_thawing():
    rng = np.random.default_rng(1)
    steps = 2000
    t = np.arange(steps)
    cold = np.where(t < 900, 255.0, 300.0)                 # a cold spell that freezes the lake, then a thaw
    f = L.Forcing(dt=3575.0, sw_down=np.clip(600 * np.sin(t / 60.0), 0, None), lw_down=200 + 60 * rng.random(steps),
                  air_k=cold + 3 * rng.standard_normal(steps), q_air=np.full(steps, 0.002),
                  wind10=4 * rng.random(steps),
                  pressure=np.full(steps, 121590.0), o2_fraction=MOON_O2)
    basin = L.Basin.from_hypsography(30.0, 1e7, [0.0, 0.5, 1.0], [1.0, 0.6, 0.0], top_m=0.5, top_depth=10.0)
    sc = [L.Scenario('eu', 5.3, 1.7, 0.37, 0.53), L.Scenario('kin', 1.6, 0.57, 0.3, 0.3, activation_ev=0.62,
                                                              transfer='cole_caraco', reference_c=6.0),
          L.Scenario('pool', 5.3, 1.7, 0.37, 0.53, labile_days=4.0)]
    lake = L.Lake(basin, L.Physics(attenuation=0.57, solute_density=0.39e-3, solute_column=2,
                                   debt_density=L.anaerobic_density(0.5), iron_density_per_debt=0.009,
                                   iron_density_cap=0.017), sc, air())
    st = L.initial_state(basin, sc, 6.0, 9.0)
    h0, o0 = L.heat_content(basin, st), L.oxygen_content(basin, st)
    out = lake.run(f, st)
    assert out['ice'].max() > 0.05 and st.ice == 0.0         # ice formed and melted again
    assert L.heat_content(basin, st) - h0 == pytest.approx(out['heat_in'], rel=1e-9, abs=1e-6 * abs(h0))
    budget = L.oxygen_content(basin, st) - o0 - out['oxygen_in'] - out['produced'] + out['respired']
    assert np.all(np.abs(budget) < 1e-9 * np.abs(o0))
    assert np.all(st.oxygen >= -1e-12)


def test_diffusion_decays_an_analytic_mode():
    # a cosine mode between no-flux boundaries decays at exp(-K (pi/H)^2 t); on the layer centres it is an exact
    # eigenvector of the implicit scheme, which decays by 1 / (1 + lambda dt) a step with
    # lambda = (4 K / dz^2) sin^2(pi dz / 2H)
    basin = L.Basin.pan(10.0, 0.05, 1.0)
    k, dt, steps = 1e-4, 600.0, 200
    mode = np.cos(np.pi * basin.centre / basin.depth)
    x = mode.copy()
    c = np.full(basin.n - 1, k) * basin.area[1:-1] / np.diff(basin.centre) * dt
    for _ in range(steps):
        L._diffuse(x, c, basin.volume)
    dz = basin.dz[0]
    discrete = (1.0 + 4.0 * k / dz ** 2 * math.sin(math.pi * dz / (2.0 * basin.depth)) ** 2 * dt) ** -steps
    assert np.allclose(x, discrete * mode, rtol=1e-10, atol=1e-12)
    analytic = math.exp(-k * (math.pi / basin.depth) ** 2 * dt * steps)
    assert np.allclose(x, analytic * mode, atol=5e-3 * analytic)           # backward Euler's first-order error


def mixing_lake(c_k=0.0, c_w=0.0, c_t=0.0):
    basin = L.Basin.pan(60.0, 0.1, 1.0)
    physics = L.Physics(c_k=c_k, c_w=c_w, c_t=c_t)
    return L.Lake(basin, physics, [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)], air()), basin


def test_convection_encroaches_into_linear_stratification():
    # cooling Q at the surface with no penetration mixes to h = sqrt(2 Q t / (rho c Gamma))
    lake, basin = mixing_lake()
    gamma, q, dt = 0.05, 100.0, 600.0
    t = 25.0 - gamma * basin.centre
    ox, debt, age = np.zeros((basin.n, 1)), np.zeros((basin.n, 1)), np.zeros(basin.n)
    h = 0.0
    steps = 1500
    for _ in range(steps):
        t[0] -= q * dt / (L.RHO_CP * basin.volume[0])
        _, _, h = lake._mix(t, ox, debt, age, 0.0, dt, h, L.MOON_SURFACE_GRAVITY)
    expected = math.sqrt(2.0 * q * dt * steps / (L.RHO_CP * gamma))
    assert h == pytest.approx(expected, rel=0.02)


def test_wind_deepens_the_mixed_layer_as_the_energy_balance_sets():
    # with GLM's energy balance and no acceleration term, h^3 = 6 C_W u*^3 t / N^2 into linear stratification
    lake, basin = mixing_lake(c_w=0.23)
    gamma, u_star, dt = 0.02, 0.01, 600.0
    t = 22.0 - gamma * basin.centre
    alpha = L._expansion(20.0)
    n2 = L.MOON_SURFACE_GRAVITY * alpha * gamma
    ox, debt, age = np.zeros((basin.n, 1)), np.zeros((basin.n, 1)), np.zeros(basin.n)
    h = 0.0
    steps = 500
    for _ in range(steps):
        _, _, h = lake._mix(t, ox, debt, age, u_star, dt, h, L.MOON_SURFACE_GRAVITY)
    expected = (6.0 * 0.23 * u_star ** 3 * dt * steps / n2) ** (1.0 / 3.0)
    assert h == pytest.approx(expected, rel=0.05)


def test_lunar_gravity_lets_the_wind_mix_deeper():
    depths = []
    for g in (L.MOON_SURFACE_GRAVITY, 9.80665):
        lake, basin = mixing_lake(c_w=0.23, c_t=0.51)
        t = 22.0 - 0.02 * basin.centre
        ox, debt, age = np.zeros((basin.n, 1)), np.zeros((basin.n, 1)), np.zeros(basin.n)
        h = 0.0
        for _ in range(2000):
            _, _, h = lake._mix(t, ox, debt, age, 0.01, 600.0, h, g)
        depths.append(h)
    assert depths[0] / depths[1] == pytest.approx((9.80665 / L.MOON_SURFACE_GRAVITY) ** (1 / 3), rel=0.1)


def test_a_mixed_slab_relaxes_to_saturation():
    basin = L.Basin.pan(5.0, 1.0, 1.0)
    sc = [L.Scenario('x', 0.0, 0.5, 0.0, 0.0, transfer='cole_caraco')]
    lake = L.Lake(basin, L.Physics(full_mixing=True, heat_flow=0.0), sc, air())
    steps = 300
    # saturated air at the water's temperature and the sky's infrared at the water's own emission hold it at 20 C
    f = steady(steps, air_k=293.15, q=L.saturation_humidity(293.15, 121590.0), wind=5.0,
               lw=L.STEFAN_BOLTZMANN * 293.15 ** 4)
    st = L.initial_state(basin, sc, 20.0, 2.0)
    out = lake.run(f, st)
    t_mean = float(out['surface_t'].mean())
    k = L.transfer_velocity('cole_caraco', 5.0, 0.0, 0.0, t_mean)
    c_sat = float(out['saturation'][-1])
    expected = c_sat + (2.0 - c_sat) * math.exp(-k * 3575.0 * steps / 5.0)
    assert out['surface_o2'][-1, 0] == pytest.approx(expected, rel=0.01)


def test_ice_grows_by_conduction_from_its_surface():
    # with the lake's top water at freezing, the base grows as d(h^2)/dt = 2 k (Tf - Ts) / (rho L)
    basin = L.Basin.pan(10.0, 0.5, 1.0)
    sc = [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)]
    lake = L.Lake(basin, L.Physics(heat_flow=0.0), sc, air())
    st = L.initial_state(basin, sc, 0.0, 10.0)
    st.ice, st.ice_surface_c = 0.05, -5.0
    h0, h2 = 0.05, 0.05 ** 2
    t = st.t
    for _ in range(400):
        ice, ts, _, _ = lake._ice_step(st.ice, st.ice_surface_c, 0.0, 150.0, 240.0, 0.0005, 3.0, 121590.0, t, 3575.0)
        rate = 2.0 * L.K_ICE * (L.FREEZING_C - ts) / (L.RHO_ICE * L.LATENT_FUSION)
        h2 += rate * 3575.0 * (st.ice / max(st.ice, 1e-3))
        st.ice, st.ice_surface_c = ice, ts
    assert st.ice > 0.3
    assert st.ice == pytest.approx(math.sqrt(h2), rel=0.02)


def test_respiration_products_hold_deep_water_against_overturn():
    # bottom water 0.05 K warmer than the water above would overturn; enough dissolved CO2 keeps it in place
    basin = L.Basin.pan(10.0, 1.0, 1.0)
    held = []
    for beta in (0.0, 1.375 * 0.284e-3):
        lake = L.Lake(basin, L.Physics(c_w=0.0, c_k=0.0, solute_density=beta), [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)],
                      air())
        t = np.full(basin.n, 20.0)
        t[-1] = 20.05
        sol = np.zeros(basin.n)
        sol[-1] = 40.0                                 # g/m3 of O2 respired: 0.016 kg/m3 against 0.010 from warmth
        ox, debt, age = np.zeros((basin.n, 1)), np.zeros((basin.n, 1)), np.zeros(basin.n)
        k, _, _ = lake._mix(t, ox, debt, age, 0.0, 600.0, 0.0, L.MOON_SURFACE_GRAVITY, sol)
        held.append(k < basin.n)
    assert held == [False, True]



def test_the_exchange_pairs_like_with_like():
    # In neutral air the stress is the same whether the coefficient is taken at 10 m with the wind there or at the
    # lowest level's height with its wind (both on the log profile); the heat exchange per degree is not, and PlaSim's
    # pairing at the level gives about three quarters of the 10 m pairing's (the review's 1.33-1.39)
    a = air()
    z_low, u_low = 870.0, 4.0
    u10 = L.reduce_wind(u_low, z_low, 10.0, a.gravity)
    ts, p = 20.0, 121590.0
    neutral = (ts + 273.15) * (1.0 + 0.608 * L.saturation_humidity(ts + 273.15, p)) / (1.0 + 0.608 * 0.010)
    _, _, _, tau_level = L.bulk_fluxes(ts, neutral, 0.010, u_low, p, a, z_low)
    _, _, _, tau_10 = L.bulk_fluxes(ts, neutral, 0.010, u10, p, a)
    assert tau_level == pytest.approx(tau_10, rel=1e-3)        # to Charnock's roughness iterated at each height
    c_level = L.neutral_coefficient(u_low, z_low, a.gravity)[0] * u_low
    c_10 = L.neutral_coefficient(u10, 10.0, a.gravity)[0] * u10
    assert 1.25 < c_10 / c_level < 1.45
    # the forcing's level is what the run pairs: giving it equals setting the air's own height
    steps = 200
    f = steady(steps, wind=u_low)
    f.wind_level, f.level_height = np.full(steps, u_low), np.full(steps, z_low)
    basin = L.Basin.pan(5.0, 1.0, 1.0)
    sc = [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)]
    runs = []
    for forcing, air_ in ((f, a), (steady(steps, wind=u_low), L.Air(gas_constant=288.38, cp=1008.9,
                                                                    wind_height_m=z_low))):
        lake = L.Lake(basin, L.Physics(heat_flow=0.0), sc, air_)
        st = L.initial_state(basin, sc, 25.0, 8.0)
        runs.append(lake.run(forcing, st)['heat_in'])
    assert runs[0] == pytest.approx(runs[1], rel=1e-12)


def test_the_anaerobic_products_density():
    # Kivu's coefficients: methane lightens the water, CO2 makes it denser; a debt half methane is denser by about
    # 1.4e-4 kg/m3 for each g/m3, a debt of CO2 alone by aerobic respiration's 0.39e-3, and of methane alone lighter
    expected = 0.284e-3 * 44.01 * (0.5 / 64 + 0.5 / 32) - 1.25e-3 * 16.04 * 0.5 / 64
    assert L.anaerobic_density(0.5) == pytest.approx(expected, rel=1e-12)
    assert L.anaerobic_density(0.5) == pytest.approx(1.364e-4, rel=2e-3)
    assert L.anaerobic_density(0.0) == pytest.approx(44.01 / 32.0 * 0.284e-3, rel=1e-12)
    assert L.anaerobic_density(1.0) < 0.0
    # the debt's products hold deep water against an overturn its warmth alone would allow
    basin = L.Basin.pan(10.0, 1.0, 1.0)
    held = []
    for beta in (0.0, L.anaerobic_density(0.5)):
        lake = L.Lake(basin, L.Physics(c_w=0.0, c_k=0.0, debt_density=beta), [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)],
                      air())
        t = np.full(basin.n, 20.0)
        t[-1] = 20.02                                  # 0.004 kg/m3 lighter from warmth
        debt = np.zeros((basin.n, 1))
        debt[-1] = 50.0                                # 0.0068 kg/m3 denser from the anaerobic products
        ox, age = np.zeros((basin.n, 1)), np.zeros(basin.n)
        k, _, _ = lake._mix(t, ox, debt, age, 0.0, 600.0, 0.0, L.MOON_SURFACE_GRAVITY)
        held.append(k < basin.n)
    assert held == [False, True]
    # iron reduction's dissolved load is capped at what a ferruginous lake holds
    lake = L.Lake(basin, L.Physics(debt_density=1e-4, iron_density_per_debt=0.009, iron_density_cap=0.017),
                  [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)], air())
    d = np.array([0.0, 1.0, 10.0])
    assert lake._debt_density(d) == pytest.approx([0.0, 1e-4 + 0.009, 1e-3 + 0.017])


def test_the_labile_pool_runs_down_through_the_night():
    # in the dark the pool decays as the implicit step has it, close to exp(-t / tau)
    basin = L.Basin.pan(4.0, 1.0, 1.0)
    sc = [L.Scenario('x', 1.0, 0.5, 0.0, 0.0, labile_days=4.0)]
    lake = L.Lake(basin, L.Physics(heat_flow=0.0), sc, air())
    steps, dt = 240, 3575.0
    f = steady(steps, dt=dt)
    lake.calibrate_production(steady(10, dt=dt, sw=400.0))
    st = L.initial_state(basin, sc, 25.0, 8.0)
    st.pool[:] = 3.0
    lake.run(f, st)
    k = 1.0 / (4.0 * 86400.0)
    assert st.pool[:, 0] == pytest.approx(3.0 / (1.0 + k * dt) ** steps, rel=1e-9)
    assert st.pool[0, 0] == pytest.approx(3.0 * math.exp(-k * dt * steps), rel=0.02)
    # with oxygen to spare, what the pool loses is respired or owed, and it holds the euphotic share of production
    basin = L.Basin.pan(4.0, 1.0, 1.0)
    sc = [L.Scenario('x', 3.0, 0.5, 0.0, 0.0, labile_days=2.0)]
    lake = L.Lake(basin, L.Physics(heat_flow=0.0), sc, air())
    t = np.arange(600)
    f = L.Forcing(dt=dt, sw_down=np.clip(500.0 * np.sin(t / 30.0), 0.0, None), lw_down=np.full(600, 400.0),
                  air_k=np.full(600, 298.0), q_air=np.full(600, 0.012), wind10=np.full(600, 3.0),
                  pressure=np.full(600, 121590.0), o2_fraction=MOON_O2)
    st = L.initial_state(basin, sc, 25.0, 8.0)
    out = lake.run(f, st)
    pool_end = float((st.pool[:, 0] * basin.volume).sum())
    debt_end = float((st.debt[:, 0] * basin.volume).sum())   # the Monod limit's small debt, repaid within a day
    assert out['respired'][0] == pytest.approx(0.8 * out['produced'][0] - pool_end - debt_end, rel=1e-9)
    assert debt_end < 0.005 * out['produced'][0]


def test_the_diffusivity_is_capped_and_set_across_instability():
    ph = L.Physics(dissipation=1e-9, max_diffusivity=1e-3, unstable_diffusivity=5e-3)
    n2 = np.array([-1e-6, 0.0, 1e-8, 1e-5, 1e-3])
    k = L.diffusivity(n2, ph)
    assert k[0] == k[1] == 5e-3
    assert k[2] == 1e-3                                # 0.2 x 1e-9 / 1e-8 = 2e-2, capped
    assert k[3] == pytest.approx(0.2 * 1e-9 / 1e-5)
    assert k[4] == pytest.approx(0.2 * 1e-9 / 1e-3)
    # the defaults keep the earlier rule, the cap and the unstable value both at 1e-4 m2/s
    d = L.Physics()
    assert d.max_diffusivity == d.unstable_diffusivity == 1e-4 and d.deep_mixing == 'dissipation'
    # unstable interfaces take at least the cap, so at a high cap no stable interface mixes faster than an unstable
    # one; at a low cap they keep the unstable value
    high = L.diffusivity(n2, L.Physics(dissipation=1e-9, max_diffusivity=1e-3))
    assert high[0] == high[1] == 1e-3 and high[2] == 1e-3
    low = L.diffusivity(n2, L.Physics(dissipation=1e-9, max_diffusivity=1e-5))
    assert low[0] == low[1] == 1e-4 and low[2] == 1e-5
    assert np.all(high[2:] <= high[0]) and np.all(low[2:] <= low[0])
    # a treatment the model does not hold is refused
    with pytest.raises(NotImplementedError):
        L.diffusivity(n2, L.Physics(deep_mixing='unknown'))


def test_the_solute_is_carried_and_mixed_without_its_density():
    # with no density of its own the solute is still mixed by an overturn and still gathers respiration's products
    basin = L.Basin.pan(10.0, 1.0, 1.0)
    lake = L.Lake(basin, L.Physics(c_w=0.0, c_k=0.0, solute_density=0.0), [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)],
                  air())
    t = np.full(basin.n, 20.0)
    t[-1] = 20.05                                      # warmer below: the column overturns
    sol = np.zeros(basin.n)
    sol[-1] = 40.0
    ox, debt, age = np.full((basin.n, 1), 5.0), np.zeros((basin.n, 1)), np.zeros(basin.n)
    k, _, _ = lake._mix(t, ox, debt, age, 0.0, 600.0, 0.0, L.MOON_SURFACE_GRAVITY, sol)
    assert k == basin.n and np.allclose(sol, 4.0)
    # and a night's respiration adds to it in every layer, the surface losing it to the air
    scen = [L.Scenario('r', 0.0, 0.5, 0.5, 0.5)]
    lake = L.Lake(basin, L.Physics(solute_density=0.0), scen, air())
    state = L.initial_state(basin, scen, 20.0, 8.0)
    lake.run(steady(24), state)
    assert np.all(state.solute[1:] > 0.0)


def test_the_arrhenius_dependence_leaves_the_rates_as_given_at_its_reference():
    basin = L.Basin.pan(5.0, 1.0, 1.0)
    sc = [L.Scenario('supply', 0.0, 0.5, 0.3, 0.3), L.Scenario('warm', 0.0, 0.5, 0.3, 0.3, activation_ev=0.62,
                                                               reference_c=6.0)]
    lake = L.Lake(basin, L.Physics(heat_flow=0.0, full_mixing=True), sc, air())
    f = steady(50, air_k=279.15, q=L.saturation_humidity(279.15, 121590.0), lw=L.STEFAN_BOLTZMANN * 279.15 ** 4,
               wind=0.0)
    st = L.initial_state(basin, sc, 6.0, 8.0)
    out = lake.run(f, st)
    assert out['respired'][1] == pytest.approx(out['respired'][0], rel=2e-3)
    st = L.initial_state(basin, sc, 25.0, 8.0)
    warm = lake.run(steady(50, air_k=298.15, q=L.saturation_humidity(298.15, 121590.0),
                           lw=L.STEFAN_BOLTZMANN * 298.15 ** 4, wind=0.0), st)
    factor = math.exp(0.62 / L.BOLTZMANN_EV * (1.0 / 279.15 - 1.0 / 298.15))
    assert warm['respired'][1] / warm['respired'][0] == pytest.approx(factor, rel=0.05)


def test_budgets_close_with_the_exchange_paired_aloft():
    rng = np.random.default_rng(2)
    steps = 1500
    t = np.arange(steps)
    f = L.Forcing(dt=3575.0, sw_down=np.clip(600 * np.sin(t / 60.0), 0, None), lw_down=380 + 40 * rng.random(steps),
                  air_k=296.0 + 3 * rng.standard_normal(steps), q_air=np.full(steps, 0.012),
                  wind10=4 * rng.random(steps), pressure=np.full(steps, 121590.0), o2_fraction=MOON_O2,
                  wind_level=1.0 + 5 * rng.random(steps), level_height=np.full(steps, 870.0))
    basin = L.Basin.from_hypsography(30.0, 1e7, [0.0, 0.5, 1.0], [1.0, 0.6, 0.0], top_m=0.5, top_depth=10.0)
    sc = [L.Scenario('pool', 5.3, 1.7, 0.37, 0.53, labile_days=4.0),
          L.Scenario('warm', 1.6, 0.57, 0.3, 0.3, activation_ev=0.62, reference_c=6.0, labile_days=2.0)]
    lake = L.Lake(basin, L.Physics(attenuation=0.57, solute_density=0.39e-3, debt_density=L.anaerobic_density(0.5),
                                   max_diffusivity=1e-3), sc, air())
    st = L.initial_state(basin, sc, 24.0, 7.0)
    h0, o0 = L.heat_content(basin, st), L.oxygen_content(basin, st)
    out = lake.run(f, st)
    assert L.heat_content(basin, st) - h0 == pytest.approx(out['heat_in'], rel=1e-9, abs=1e-6 * abs(h0))
    budget = L.oxygen_content(basin, st) - o0 - out['oxygen_in'] - out['produced'] + out['respired']
    assert np.all(np.abs(budget) < 1e-9 * np.abs(o0))
    assert np.all(st.oxygen >= -1e-12) and np.all(st.pool >= 0.0)



def test_the_wind_scheme_spreads_a_share_of_the_winds_energy_over_the_water_below():
    # the drag coefficient at 10 m in Wuest & Lorke's (2003) pieces, which join where they meet
    assert L.drag_coefficient_10m(0.1) == pytest.approx(6.215e-2, rel=1e-3)
    assert float(L.drag_coefficient_10m(3.85)) == pytest.approx(-7.12e-7 * 3.85 ** 2 + 7.387e-5 * 3.85 + 6.605e-4,
                                                                 rel=2e-3)
    assert L.drag_coefficient_10m(5.0) == pytest.approx(-7.12e-7 * 25.0 + 7.387e-5 * 5.0 + 6.605e-4)
    # P10 = rho_air C10 U10^3 in W/m2, and the dissipation's mass balance: eps_h M_h = eta P10
    p10 = float(L.wind_energy_flux(5.0, 1.2))
    assert p10 == pytest.approx(1.2 * 1.01209e-3 * 125.0, rel=1e-4)
    ph = L.Physics(deep_mixing='wind', wind_share=0.003, max_diffusivity=1.0)
    mass = 1000.0 * 100.0                              # 100 m of water below the mixed layer, kg/m2
    eps = L.wind_dissipation(p10, mass, ph)
    assert eps * mass == pytest.approx(ph.wind_share * p10, rel=1e-12)
    # dimensions: (W/m2) / (kg/m2) is W/kg, m2/s3; 0.2 of it over N2 (1/s2) is a diffusivity in m2/s
    assert eps == pytest.approx(4.55e-9, rel=2e-3)
    n2 = np.array([1e-6])
    assert L.diffusivity(n2, ph, eps)[0] == pytest.approx(0.2 * eps / 1e-6)
    assert L.wind_dissipation(p10, 0.0, ph) == 0.0     # no water below the mixed layer, no dissipation there


def test_the_wind_scheme_reproduces_the_fixed_dissipation():
    n2 = np.array([-1e-7, 0.0, 1e-9, 1e-7, 1e-5, 1e-3])
    fixed = L.Physics(dissipation=3e-10, max_diffusivity=1e-2)
    wind = L.Physics(deep_mixing='wind', max_diffusivity=1e-2)
    assert np.array_equal(L.diffusivity(n2, wind, 3e-10), L.diffusivity(n2, fixed))
    # and in a run: a wind share set so that eta P10 / M_h is the fixed dissipation gives the fixed run's column
    steps = 200
    f = L.Forcing(dt=3575.0, sw_down=np.zeros(steps), lw_down=np.full(steps, 480.0), air_k=np.full(steps, 312.0),
                  q_air=np.full(steps, 0.02), wind10=np.full(steps, 5.0), pressure=np.full(steps, 121590.0),
                  o2_fraction=MOON_O2)
    basin = L.Basin.pan(40.0, 1.0, 1.0)
    sc = [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)]
    t0 = 20.0 + 6.0 * np.exp(-np.arange(basin.n) / 4.0)     # warm above, a thermocline the calm night keeps
    rho_air = 121590.0 / (air().gas_constant * 312.0)
    p10 = float(L.wind_energy_flux(5.0, rho_air))
    out = {}
    for name, ph in (('fixed', L.Physics(c_w=0.0, c_k=0.0, dissipation=1e-9, max_diffusivity=1e-2)),
                     ('wind', None)):
        if ph is None:                                 # the share that makes eta P10 / M_h equal 1e-9 W/kg
            below = L.RHO0 * 39.0                     # the warming surface and c_w = c_k = 0 keep the mixed layer
                                                       # in the top layer
            ph = L.Physics(c_w=0.0, c_k=0.0, deep_mixing='wind', wind_share=1e-9 * below / p10, max_diffusivity=1e-2)
        lake = L.Lake(basin, ph, sc, air())
        st = L.initial_state(basin, sc, t0, 8.0)
        r = lake.run(f, st)
        out[name] = (st.t.copy(), r)
    assert np.allclose(out['wind'][0], out['fixed'][0], rtol=0.0, atol=1e-6)
    assert np.allclose(out['wind'][1]['dissipation'], 1e-9, rtol=1e-3)


def test_the_raised_cap_does_not_bind_under_the_wind():
    rng = np.random.default_rng(3)
    steps = 1500
    t = np.arange(steps)
    f = L.Forcing(dt=7150.0, sw_down=np.clip(600 * np.sin(t / 30.0), 0, None), lw_down=380 + 40 * rng.random(steps),
                  air_k=296.0 + 3 * rng.standard_normal(steps), q_air=np.full(steps, 0.012),
                  wind10=5 * rng.random(steps), pressure=np.full(steps, 121590.0), o2_fraction=MOON_O2)
    basin = L.Basin.from_hypsography(100.0, 1e7, [0.0, 0.5, 1.0], [1.0, 0.6, 0.0], top_m=0.5, top_depth=10.0)
    sc = [L.Scenario('meso', 1.6, 0.57, 0.3, 0.3)]
    t0 = 20.0 + 6.0 * np.exp(-basin.centre / 10.0)
    held, column = {}, {}
    for cap in (1e-4, 1e-1, 1.0):
        lake = L.Lake(basin, L.Physics(deep_mixing='wind', max_diffusivity=cap, solute_density=0.39e-3), sc, air())
        st = L.initial_state(basin, sc, t0, 8.0)
        h0, o0 = L.heat_content(basin, st), L.oxygen_content(basin, st)
        out = lake.run(f, st)
        held[cap] = out['interfaces_held_by_cap'] / out['interfaces_below']
        column[cap] = (float(st.t[-1]), float(out['deep_o2'][-1, 0]))
        assert out['interfaces_held_by_cap'] <= out['interfaces_capped']
        # the budgets close under the wind's mixing as under the fixed dissipation's
        assert L.heat_content(basin, st) - h0 == pytest.approx(out['heat_in'], rel=1e-9, abs=1e-6 * abs(h0))
        budget = L.oxygen_content(basin, st) - o0 - out['oxygen_in'] - out['produced'] + out['respired']
        assert np.all(np.abs(budget) < 1e-9 * np.abs(o0))
    # at 0.1 m2/s the cap holds no interface: where it is reached, near-neutral water, the cap still mixes more in a
    # step than the step can, and raising it tenfold leaves the deep water as it was; the old cap held a quarter
    assert held[1e-1] == 0.0 and held[1e-4] > 0.1
    assert column[1.0][0] == pytest.approx(column[1e-1][0], abs=0.01)
    assert column[1.0][1] == pytest.approx(column[1e-1][1], abs=0.01)


def test_the_lake_number_diffusivity():
    ph = L.Physics(deep_mixing='lake_number', max_diffusivity=1.0)
    n2 = np.array([-1e-6, 1e-6, 1e-4, 4e-4])
    k1 = L.lake_number_diffusivity(n2, ph, 1.0)
    assert k1[3] == pytest.approx(200.0 * L.KAPPA_HEAT)       # at the thermocline, L_N = 1: 2.8e-5 m2/s
    assert k1[2] == pytest.approx(k1[3] / 4.0) and k1[0] == 1.0  # with N2, and unstable at the cap
    assert np.allclose(L.lake_number_diffusivity(n2, ph, 0.5)[1:], 2.0 * k1[1:])
    assert np.all(L.lake_number_diffusivity(n2, ph, math.inf)[1:] == 0.0)
    basin = L.Basin.from_hypsography(100.0, 1e7, [0.0, 0.5, 1.0], [1.0, 0.6, 0.0], top_m=0.5, top_depth=10.0)
    lake = L.Lake(basin, ph, [L.Scenario('x', 0.0, 0.5, 0.0, 0.0)], air())
    rho = L.density(20.0 + 6.0 * np.exp(-basin.centre / 10.0))
    ln = lake.lake_number(rho, 1e-5, L.MOON_SURFACE_GRAVITY)
    assert 0.0 < ln < math.inf and lake.lake_number(rho, 0.0, L.MOON_SURFACE_GRAVITY) == math.inf
    assert lake.lake_number(rho, 2e-5, L.MOON_SURFACE_GRAVITY) == pytest.approx(ln / 2.0)
