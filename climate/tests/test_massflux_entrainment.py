"""The entrainment that follows the environment (massflux.PARAMETERS entrainment_form 'environment', 27 October 2026)
on made-up columns: it is the IFS's (1.3 - RH) (qs / qs_base)^3 form above the CM1 boxes' rate aloft, shallow clouds
twice it; the plume takes the saturation at its own condensation level as the base's and the updraught carries the test
parcel's; the fixed profile stays the default and the only form the Fortran carries; the testbed takes the scheme's
parameters and keeps each step's convection when asked, and the judgement reads each criterion as it was set."""
import numpy as np
import pytest

from climate.gcm import massflux as mf
from climate.gcm import massflux_entrainment as me
from climate.gcm import massflux_fortran as mff
from climate.gcm import pbl_column as pc
from climate.tests.test_massflux import C, G, lunar_column, run
from climate.tests.test_pbl_column import made_up_rows

ENV = me.ENVIRONMENT


def test_the_environment_form_is_the_ifs_form_above_the_rate_aloft():
    eps = lambda kind, rh, qs, qsb: mf.entrainment(kind, 5000.0, rh, 3000.0, ENV, qs, qsb)
    for rh, qs, qsb in ((0.97, 0.009, 0.010), (0.7, 0.004, 0.011), (0.3, 0.011, 0.011)):
        expected = ENV['eps_bg'] + ENV['eps_env'] * (1.3 - rh) * (qs / qsb) ** 3
        assert eps('deep', rh, qs, qsb) == pytest.approx(expected, rel=1e-14)
        assert eps('shallow', rh, qs, qsb) == pytest.approx(ENV['shallow_factor'] * expected, rel=1e-14)
        assert eps('mid', rh, qs, qsb) == eps('deep', rh, qs, qsb)          # mid-level clouds take the deep rate
    assert eps('deep', 1.2, 0.01, 0.01) == eps('deep', 1.0, 0.01, 0.01)    # a supersaturated surrounding counts as 1
    assert eps('deep', 0.8, None, None) == pytest.approx(ENV['eps_bg'] + ENV['eps_env'] * 0.5)
    assert ENV['eps_env'] == 1.75e-3 and ENV['eps_bg'] == mf.PARAMETERS['eps_bg']
    # the fixed profile reads neither saturation
    fixed = mf.entrainment('deep', 5000.0, 0.8, 3000.0, mf.PARAMETERS)
    assert mf.entrainment('deep', 5000.0, 0.8, 3000.0, mf.PARAMETERS, 0.001, 0.01) == fixed
    assert fixed == pytest.approx(0.8e-4 + 2.0e-4 * np.exp(-1.0))


def test_the_updraught_takes_its_test_parcels_base_saturation():
    lc = lunar_column(rh=0.95)                       # moist enough above its base to convect deep under the form
    col = mf.Column(lc['t'], lc['q'], lc['p_full'], lc['p_half'], G, C, lc['z_full'], lc['z_half'])
    decision = mf.classify(col, 10.0, 100.0, None, ENV)
    assert decision['type'] == 'deep'
    tr = decision['trial']
    expected = tr['points'][tr['base']]['qs_env']
    # the condensation point's saturation is the surrounding air's there
    z_b, p_b = tr['points'][tr['base']]['z'], tr['points'][tr['base']]['p']
    t_env = np.interp(z_b, col.zf[::-1], col.t[::-1])
    assert expected == pytest.approx(mf.qsat(t_env, p_b, C), rel=2e-2)
    calls = []
    original = mf.entrainment

    def spy(kind, z, rh, z_base, prm=mf.PARAMETERS, qs=None, qs_base=None):
        calls.append((kind, qs, qs_base))
        return original(kind, z, rh, z_base, prm, qs, qs_base)
    mf.entrainment = spy
    try:
        up = mf.updraught(col, decision, ENV)
    finally:
        mf.entrainment = original
    assert calls and all(c[2] == expected for c in calls) and all(c[0] == 'deep' for c in calls)
    assert all(c[1] <= expected * 1.05 for c in calls)                      # the air cools above the base
    cloudy = [p['eps'] for p in up['points'] if p['eps'] > 0.0]
    assert max(cloudy) > 2.0 * max(p['eps'] for p in mf.updraught(col, decision, mf.PARAMETERS)['points'])


def test_the_fixed_profile_stays_the_default_and_the_only_form_the_fortran_carries():
    assert mf.PARAMETERS['entrainment_form'] == 'fixed'
    assert mff.PYTHON_ONLY['entrainment_form'] == 'fixed'
    with pytest.raises(ValueError):
        mff.fortran_values(dict(mf.PARAMETERS, entrainment_form='environment'))
    col = lunar_column()
    a = run(col)
    b = run(col, prm=dict(mf.PARAMETERS, eps_env=99.0))                     # the fixed profile never reads it
    assert a['type'] == b['type'] == 'deep'
    assert np.array_equal(a['dtdt'], b['dtdt']) and np.array_equal(a['dqdt'], b['dqdt'])
    moist = lunar_column(rh=0.95)
    assert run(moist)['type'] == run(moist, prm=ENV)['type'] == 'deep'
    assert not np.array_equal(run(moist)['dtdt'], run(moist, prm=ENV)['dtdt'])


def test_the_testbed_takes_the_schemes_parameters_and_keeps_its_steps():
    rows, _, _ = made_up_rows()
    rows = [dict(r, quarter=mf.quarter_name(r['hour'])) for r in rows]      # as the testbed's gridded rows carry it
    half = np.asarray(pc.grid(), float)
    plain = pc.run_column(rows, 'night', 'night', 'sea', me.SCHEME, half)
    same = pc.run_column(rows, 'night', 'night', 'sea', me.SCHEME, half, prm=mf.PARAMETERS, keep_steps=True)
    a, b = plain['quarters']['night'], same['quarters']['night']
    assert np.array_equal(a['t_k'], b['t_k']) and np.array_equal(a['q_g_kg'], b['q_g_kg'])
    assert 'steps_record' not in a and len(b['steps_record']) == b['steps']
    assert set(b['steps_record'][0]) == {'type', 'top_m', 'base_m', 'mass_flux', 'rain_conv', 'rain_ls', 'hour'}
    env = pc.run_column(rows, 'night', 'night', 'sea', me.SCHEME, half, prm=ENV, keep_steps=True)
    assert len(env['quarters']['night']['steps_record']) == b['steps']


def test_the_form_profile_is_the_form_at_each_level():
    t = np.array([300.0, 296.0, 292.0, 288.0])
    p = np.array([121000.0, 115000.0, 109000.0, 103000.0])
    qs = np.array([mf.qsat(a, b) for a, b in zip(t, p)])
    q = qs * np.array([0.7, 0.9, 0.95, 0.8])
    prof = me.form_profile(t, q, p, 1, 'deep')
    expected = ENV['eps_bg'] + ENV['eps_env'] * (1.3 - q / qs) * (qs / qs[1]) ** 3
    assert prof == pytest.approx(expected, rel=1e-12)
    assert me.form_profile(t, q, p, 1, 'shallow') == pytest.approx(2.0 * expected, rel=1e-12)


def made_up_columns(night_tops, depth, drop, lid, rains, upper_new, upper_old, above20=(0.9, 0.9)):
    """The judgement's input: the four columns with only what it reads."""
    def entry(rain, upper, tops10=0.0, tops20=None, measures=None):
        return dict(rain_ground_mm_day=rain, tops_below_10km_share=tops10, convecting_tops_above_20km_share=tops20,
                    measures=measures or {}, upper_air=[dict(minus_cm1_k=u) for u in upper])
    cols = {}
    names = ('sea afternoon', 'sea evening', 'sea night', 'land afternoon')
    for name, rain in zip(names, rains):
        cols[name] = dict(forms=dict(fixed=entry(rain, upper_old), environment=entry(rain, upper_new)))
    cols['sea night']['forms']['environment'].update(
        tops_below_10km_share=night_tops, measures=dict(mixed_depth_km=depth, vapour_drop_g_kg=drop, lid_rise_k=lid))
    for name, v in zip(('sea afternoon', 'land afternoon'), above20):
        cols[name]['forms']['environment']['convecting_tops_above_20km_share'] = v
    return cols


def test_the_judgement_reads_each_criterion_as_it_was_set():
    c = me.CRITERIA
    edge = made_up_columns(c['night_tops_below_10km_share'], c['night_mixed_depth_km'], c['night_vapour_drop_g_kg'],
                           c['night_lid_rise_k'], (c['sea_afternoon_rain_mm_day'][0], c['sea_evening_rain_mm_day'],
                                                   c['night_rain_mm_day'], c['land_afternoon_rain_mm_day'][1]),
                           upper_new=(1.1, -2.1), upper_old=(1.0, 2.0), above20=(0.5, 0.5))
    j = me.judge(edge)
    assert j['passes'] and j['night_passes'] and j['day_passes']           # every bound counts as met
    worse = made_up_columns(0.7, 6.0, 1.0, 1.0, (2.0, 2.5, 0.2, 8.0), upper_new=(1.2, 0.5), upper_old=(1.0, 0.4))
    j = me.judge(worse)
    assert not j['criteria']['night_tops_below_10km']['met'] and not j['night_passes']
    assert not j['criteria']['sea_afternoon_upper_air_41.8km']['met']      # 0.2 K further than the fixed profile
    assert j['criteria']['sea_afternoon_upper_air_55.8km']['met']
    none = made_up_columns(0.9, 6.0, 1.0, 1.0, (2.0, 2.5, 0.2, 8.0), upper_new=(0.0, 0.0), upper_old=(0.0, 0.0),
                           above20=(None, 0.9))
    j = me.judge(none)
    assert not j['criteria']['sea_afternoon_tops_above_20km']['met'] and j['night_passes'] and not j['day_passes']
    assert c['night_mixed_depth_km'] == pytest.approx(5.675) and c['night_vapour_drop_g_kg'] == pytest.approx(1.84)
    assert c['sea_evening_rain_mm_day'] == pytest.approx(2.0 / 3.0 * 2.96)
