"""The lake model against values and solutions that can be checked: oxygen at saturation, conservation of heat and
oxygen, an analytic diffusion mode, convective and wind deepening into linear stratification, reaeration and ice
growth by conduction."""
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
                                                              transfer='cole_caraco')]
    lake = L.Lake(basin, L.Physics(attenuation=0.57), sc, air())
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
        sol[-1] = 40.0                                 # g/m3 of O2 respired: 0.021 kg/m3 against 0.010 from warmth
        ox, debt, age = np.zeros((basin.n, 1)), np.zeros((basin.n, 1)), np.zeros(basin.n)
        k, _, _ = lake._mix(t, ox, debt, age, 0.0, 600.0, 0.0, L.MOON_SURFACE_GRAVITY, sol)
        held.append(k < basin.n)
    assert held == [False, True]
