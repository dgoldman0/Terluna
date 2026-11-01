"""The single column's convection setups (climate/gcm/pbl_column.py) on made-up columns: the shallow-only setting keeps
the mass-flux scheme off its deep and mid-level branches, a transport's flux is recovered from the ground up past a
heating aloft that is no flux's convergence, and a step with CM1's prescribed convection closes its water and heat
budgets to round-off and carries the prescribed transport's flux."""
import numpy as np
import pytest

from climate.gcm import levels as lv
from climate.gcm import massflux as mf
from climate.gcm import pbl_column as pc
from climate.tests.test_massflux import lunar_column, run


def test_the_shallow_only_setting_keeps_the_scheme_off_its_deep_and_mid_level_branches():
    deep = lunar_column()
    assert run(deep)['type'] == 'deep'
    limited = run(deep, omega=-0.05 * np.ones(60), prm=dict(mf.PARAMETERS, **pc.SHALLOW_ONLY))
    assert limited['type'] in (None, 'shallow')
    assert limited['type'] is None or limited['shallow_trial']['depth'] <= mf.PARAMETERS['deep_dp']
    capped = lunar_column(stable_layer=True, n=200)
    kw = dict(shf=10.0, lhf=150.0)
    assert run(capped, **kw)['type'] == run(capped, prm=dict(mf.PARAMETERS, **pc.SHALLOW_ONLY), **kw)['type'] == 'shallow'
    assert run(capped, **kw)['m_base'] == run(capped, prm=dict(mf.PARAMETERS, **pc.SHALLOW_ONLY), **kw)['m_base']


def test_a_transport_s_flux_comes_from_the_ground_up_past_a_heating_aloft():
    c = pc.constants(1.6242)
    mass = np.array([10.0, 20.0, 30.0, 40.0])                       # top first
    flux = np.array([3.0, -2.0, 5.0])                                # W/m2 up through the interfaces, top first
    into = np.r_[flux, 0.0] - np.r_[0.0, flux]                       # each layer gains what comes from below less above
    dtdt = into / (c['acpd'] * mass)
    dqdt = 2.0 * into / (c['ALV'] * mass)
    heat, vapour = pc.interface_flux_from_ground(dtdt, dqdt, mass, c)
    assert heat == pytest.approx(flux) and vapour == pytest.approx(2.0 * flux)
    assert pc.interface_flux(dtdt, dqdt, mass, c)[0] == pytest.approx(flux)      # the same from the top when conservative
    dtdt[0] += 7.0 / (c['acpd'] * mass[0])                           # a heating of 7 W/m2 in the top layer alone
    assert pc.interface_flux_from_ground(dtdt, dqdt, mass, c)[0] == pytest.approx(flux)
    assert pc.interface_flux(dtdt, dqdt, mass, c)[0] == pytest.approx(flux + 7.0)


def made_up_rows(g=1.6242, ps=1.2e5):
    """Two outputs of a stable, moist column on PlaSim's fifteen levels, three hours apart after midnight, with a
    prescribed condensation and a transport that carries heat and vapour up through the lowest interfaces."""
    half = np.asarray(lv.half_levels(pc.LEVEL_SET), float)
    ph = np.r_[0.0, half * ps]
    pf = 0.5 * (ph[:-1] + ph[1:])
    sigma = pf / ps
    n = sigma.size
    c = pc.constants(g)
    mass = np.diff(ph) / g
    t = np.maximum(300.0 * sigma ** 0.25, 200.0)                     # stable throughout
    q = 0.012 * sigma ** 8                                           # unsaturated at every level
    flux = np.zeros(n - 1)
    flux[-4:] = [2.0, 6.0, 9.0, 5.0]                                 # W/m2 up through the lowest four interfaces
    into = np.r_[flux, 0.0] - np.r_[0.0, flux]
    rows = []
    for k, hour in enumerate((-150.0, -105.0)):
        rows.append(dict(t=t, q=q, p_full=pf, p_half=ph, g=g, day=60.0 + k * 0.125, hour=hour, shf=8.0, lhf=90.0,
                         ust=0.2, u=np.zeros(n), v=np.zeros(n), rad_dtdt=np.full(n, -1.0e-5), w_mean=np.zeros(n),
                         table_dtdt=np.zeros(n), table_dqdt=np.zeros(n),
                         cm1_mp_dtdt=np.where(np.arange(n) == n - 3, 2.0e-5, 0.0),
                         cm1_mp_dqdt=np.where(np.arange(n) == n - 3, -2.0e-5 * c['acpd'] / c['ALV'], 0.0),
                         cm1_transport_dtdt=into / (c['acpd'] * mass), cm1_transport_dqdt=into / (c['ALV'] * mass)))
    return rows, flux, c


def test_a_step_with_cm1_s_convection_closes_its_budgets_and_carries_the_prescribed_flux():
    rows, flux, c = made_up_rows()
    r0 = rows[0]
    pf, ph, g = r0['p_full'], r0['p_half'], r0['g']
    sigma, sigmah = pf / ph[-1], ph[1:] / ph[-1]
    mass = np.diff(ph) / g
    times = np.array([r['day'] for r in rows]) * pc.DAY
    t, q, rec = pc.step(r0['t'].copy(), r0['q'].copy(), rows, times, times[0], sigma, sigmah, pf, ph, mass, g, c, False,
                        'sea', pc.ps.LAMBDA_EARTH_M, False, 'cm1')
    b = rec['budget']
    water = b['evaporation'] + b['large_scale_water'] + b['prescribed_water'] + b['floor_water'] - b['rain'] - b['water_stored']
    heat = (b['sensible'] + b['radiation'] + b['large_scale_heat'] + b['prescribed_heat'] + c['ALV'] * b['rain']
            - b['heat_stored'])
    assert abs(water) < 1e-12 and abs(heat) < 1e-6 * abs(b['radiation'])
    assert b['prescribed_water'] == pytest.approx(float(np.sum(mass * rows[0]['cm1_mp_dqdt'])))   # the transport keeps it
    assert rec['heat_transport'] == pytest.approx(flux) and rec['vap_transport'] == pytest.approx(flux)
    assert rec['conv_type'] is None and rec['rain_conv'] == 0.0 and rec['quarter'] == 'night'
    with pytest.raises(ValueError):
        pc.step(r0['t'].copy(), r0['q'].copy(), rows, times, times[0], sigma, sigmah, pf, ph, mass, g, c, False, 'sea',
                pc.ps.LAMBDA_EARTH_M, False, 'kuo')


def test_the_shallow_night_setup_limits_the_scheme_only_while_the_sun_is_down():
    rows, _, c = made_up_rows()
    r0 = rows[0]
    pf, ph, g = r0['p_full'], r0['p_half'], r0['g']
    sigma, sigmah = pf / ph[-1], ph[1:] / ph[-1]
    mass = np.diff(ph) / g
    times = np.array([r['day'] for r in rows]) * pc.DAY
    args = (sigma, sigmah, pf, ph, mass, g, c, False, 'sea', pc.ps.LAMBDA_EARTH_M, False, 'massflux_shallow_night')
    assert pc.step(r0['t'].copy(), r0['q'].copy(), rows, times, times[0], *args)[2]['shallow_only']
    day = [dict(r, hour=h) for r, h in zip(rows, (30.0, 52.5))]
    assert not pc.step(r0['t'].copy(), r0['q'].copy(), day, times, times[0], *args)[2]['shallow_only']
