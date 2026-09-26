"""Checks on what the sunlit exosphere loses: rates, fragments, charge exchange and the magnetosphere's partition."""
import json
import math

import numpy as np
import pytest

from shared.constants import BOLTZMANN, MOON_GM, MOON_RADIUS
from atmosphere.loss_response import absorption as ab
from atmosphere.loss_response import exosphere as ex

PARTS = ('outside', 'open', 'convected', 'recombined_escaping', 'recombined_staying', 'drained')


def exobase_profile(t=250.0, r_top=3.0, n_top=4e12):
    """An isothermal column whose top is an exobase at r_top lunar radii."""
    k = MOON_GM * ab.AIR_KG / (BOLTZMANN * t * MOON_RADIUS)
    r = np.linspace(2.0, r_top, 50)
    n = n_top * np.exp(k / r - k / r_top)
    return dict(radius_R=r, temperature_K=np.full_like(r, t), pressure_Pa=n * BOLTZMANN * t)


def test_dissociation_rates_agree_with_heays():
    quiet = ex.rates('quiet')
    for species in ('N2', 'O2'):
        assert quiet[species]['dissociation'] == pytest.approx(ex.HEAYS_1AU_S[species], rel=0.1)
    ratio = ex.rates('solar_maximum')['O2']['dissociation'] / quiet['O2']['dissociation']
    assert ratio == pytest.approx(ex.FAR_UV_ACTIVITY['solar_maximum'])


def test_escaping_share_limits():
    r_c = 3.0
    fast = ex.escaping_share(np.array([3.0, 4.0, 30.0, 1e4]), r_c, 1000.0, ex.ATOM_KG['O2'])
    assert fast[0] == pytest.approx(0.5, abs=1e-6)
    assert fast[1] == pytest.approx(0.5 * (1 + math.sqrt(1 - (3 / 4) ** 2)), rel=1e-3)   # straight paths
    assert fast[-1] == pytest.approx(1.0, abs=1e-6)
    assert np.all(np.diff(fast) > 0)
    # an atom below the escape speed stays; gravity bends a slow escaping one into the exobase more often
    assert np.all(ex.escaping_share(np.array([3.0, 4.0]), r_c, 0.05, ex.ATOM_KG['O2']) == 0.0)
    assert 0.5 < ex.escaping_share(np.array([4.0]), r_c, 0.3, ex.ATOM_KG['O2'])[0] < fast[1]


def test_diffusive_separation_lowers_the_oxygen_share():
    profile = exobase_profile()
    assert ex.exobase_shares(profile)['O2'] == pytest.approx(ab.O2_SHARE)
    assert 0.01 < ex.exobase_shares(profile, ex.HOMOPAUSE_PA)['O2'] < 0.06
    assert ex.exobase_shares(profile, homopause_pa=1e-12)['O2'] == pytest.approx(ab.O2_SHARE)


def test_ions_made_match_the_absorption_step():
    profile = exobase_profile()
    made = lambda shadow: sum(v['ions'][0] for v in ex.sunlit(profile, 'quiet', [shadow]).values())
    mean_air = lambda shadow: ab.ion_production(profile, 'quiet', [shadow], ex.OUTER_R)[0]
    for shadow in (2.9, 3.5):
        assert made(shadow) == pytest.approx(mean_air(shadow), rel=0.05)
    # far out N2, lighter than the mean air, reaches farther
    assert 1.0 < made(6.0) / mean_air(6.0) < 1.15


def test_charge_exchange_thin_limit_counts_the_swept_exosphere(monkeypatch):
    profile = exobase_profile()
    monkeypatch.setattr(ex, 'CHARGE_EXCHANGE_M2', 1e-26)
    got = ex.charge_exchange(profile)
    # every ballistic molecule the wind passes, all but those in the cylinder downstream of the exobase
    exo = ex.Exosphere(profile, 'quiet')
    swept = np.array([np.sum(exo.cell * ~((exo.mu < 0) & (ri * exo.sin_t < exo.r_c))) for ri in exo.r])
    flux = ex.lr.SOLAR_WIND['density_m3'] * ex.lr.SOLAR_WIND['speed_m_s']
    expected = flux * 1e-26 * np.trapezoid(swept * exo.molecule_mass * exo.r ** 2 * MOON_RADIUS ** 3, exo.r)
    assert got == pytest.approx(expected, rel=0.03)


def test_standoff_is_the_septembers_ten_radii_and_scales_as_the_cube_root():
    assert ex.standoff_R(ex.SEPTEMBER_MOMENT) == pytest.approx(9.95, abs=0.1)
    assert ex.standoff_R(8 * ex.SEPTEMBER_MOMENT) == pytest.approx(2 * ex.standoff_R(ex.SEPTEMBER_MOMENT))


def test_magnetosphere_accounts_for_every_ion_and_a_stronger_one_keeps_more():
    profile = exobase_profile()
    shadow = 3.0
    made = sum(v['ions'][0] for v in ex.sunlit(profile, 'quiet', [shadow]).values())
    exo = ex.Exosphere(profile, 'quiet')
    assert exo.fates(shadow, 0.0)['outside'] == pytest.approx(made, rel=0.02)
    fates = {m: exo.fates(shadow, m) for m in (1e20, ex.SEPTEMBER_MOMENT, 1e23)}
    for f in fates.values():
        assert sum(f[k] for k in PARTS) == pytest.approx(made, rel=0.02)
    lost = lambda f: f['outside'] + f['open'] + f['convected']
    assert lost(fates[1e23]) < lost(fates[ex.SEPTEMBER_MOMENT]) < lost(fates[1e20])


def test_levels_are_ordered_and_wider_shadows_never_lose_more():
    res = ex.assess(exobase_profile(), 'quiet', 2.95, budget=10.0)
    for name in ('no_magnetosphere', 'september_magnetosphere'):
        totals = [res[name][lv]['total_kg_s'] for lv in ex.LEVELS]
        assert totals[0] <= totals[1] <= totals[2]
    assert res['september_magnetosphere']['central']['total_kg_s'] < res['no_magnetosphere']['central']['total_kg_s']
    for key in ('no_magnetosphere_kg_s', 'september_magnetosphere_kg_s', 'fragments_kg_s'):
        assert np.all(np.diff(res['against_radius'][key]) <= 1e-12)
    sweep = res['moment_sweep']['loss_kg_s']
    assert sweep[-1] < sweep[0]


def test_largest_transmission_stops_at_the_budget_or_the_heating_radius():
    f = [1e-4, 1e-3, 1e-2]
    heating = [2.0, 2.5, 3.5]
    # the loss crosses the budget between the grid points (log-log)
    assert ex._largest_within(f, [0.1, 1.0, 10.0], heating, 4.0, 3.0) == pytest.approx(3e-3)
    # the heating radius passes the protected radius first
    assert ex._largest_within(f, [0.1, 1.0, 10.0], heating, 3.0, 5.0) == pytest.approx(10 ** -2.5)
    assert ex._largest_within(f, [2.0, 3.0, 4.0], heating, 4.0, 1.0) is None
    assert ex._largest_within(f, [0.1, 0.2, math.inf], [2.0, 2.5, math.inf], 4.0, 1.0) == pytest.approx(1e-3)
    assert ex._largest_within(f, [0.1, 0.2, 0.3], heating, 4.0, 1.0) == pytest.approx(1e-2)


def test_stored_results():
    data = json.loads((ex.HERE / 'results' / 'exosphere_loss.json').read_text())
    assert data['schema'] == ex.SCHEMA
    assert data['cases']
    for row in data['cases']:
        none, september = row['no_magnetosphere']['central'], row['september_magnetosphere']['central']
        assert september['total_kg_s'] <= none['total_kg_s']
        assert september['fragments_kg_s'] == pytest.approx(none['fragments_kg_s'])
