"""The lid split (climate/crm/box_lid.py) on made-up domain diagnostics: CM1's fluxes in W/m2 by part, the levels nearest
the heights asked, the means, totals and shares over a span, and the quarters of the lunar day."""
import numpy as np
import pytest

from climate.crm import box_lid as bl
from climate.crm import cm1_output as co


def made_up_output(mtime, scale=1.0):
    """One output's scalar-level and w-level files on four w levels (0, 1, 6.45 and 6.99 km)."""
    s = dict(mtime=mtime, psfc=1.2e5)
    w = dict(prsinterp=np.array([1.0e5, 1.1e5, 1.0e5, 0.95e5]), rhointerp=np.array([1.4, 1.3, 1.2, 1.1]),
             thfr=scale * np.array([0.0, 0.01, -0.001, -0.0005]), thfd=scale * np.array([0.0, 0.0, -0.0005, -0.0005]),
             thfs=scale * np.array([0.01, 0.002, -0.008, -0.009]),
             qvfr=scale * np.array([0.0, 2e-6, 3e-6, 2e-6]), qvfd=scale * np.array([0.0, 0.0, 3e-7, 4e-7]),
             qvfs=scale * np.array([4e-5, 1e-6, 5e-6, 5e-6]), khv=np.array([0.0, 20.0, 10.0, 8.0]))
    return s, w


def test_each_part_turns_into_watts_by_the_air_s_density_heat_capacity_and_exner_function():
    s, w = made_up_output(100.0)
    f = bl.output_fluxes(s, w)
    prs = np.array([1.2e5, 1.1e5, 1.0e5, 0.95e5])                        # the ground's pressure at the lowest level
    to_heat = w['rhointerp'] * co.CP * (prs / co.P00) ** (co.RD / co.CP)
    assert f['heat_subgrid'] == pytest.approx(to_heat * w['thfs'])
    assert f['heat_resolved'] == pytest.approx(to_heat * w['thfr']) and f['heat_diffusion'] == pytest.approx(to_heat * w['thfd'])
    assert f['vapour_subgrid'] == pytest.approx(w['rhointerp'] * 2.501e6 * w['qvfs'])
    assert f['mtime'] == 100.0 and f['khv'] == pytest.approx(w['khv'])


def test_the_split_averages_each_part_at_the_nearest_levels_and_shares_the_total():
    zw = [0.0, 1000.0, 6446.0, 6989.0]
    rows = [bl.output_fluxes(*made_up_output(t, scale)) for t, scale in ((1.0, 1.0), (2.0, 3.0))]
    out = bl.split(rows, zw, (6450.0, 7100.0))
    assert [lev['height_m'] for lev in out] == [6446.0, 6989.0] and out[0]['outputs'] == 2
    lev = out[0]
    heat = {p: np.mean([r[f'heat_{p}'][2] for r in rows]) for p in bl.PARTS}
    assert lev['heat']['subgrid'] == pytest.approx(heat['subgrid'])
    assert lev['heat']['total'] == pytest.approx(sum(heat.values()))
    assert lev['heat']['subgrid_share'] == pytest.approx(0.008 / 0.0095)                 # the same air at both outputs
    assert lev['heat']['subgrid_and_diffusion_share'] == pytest.approx(0.0085 / 0.0095)
    assert lev['heat']['resolved_advection'] == pytest.approx(heat['resolved'] + heat['diffusion'])
    assert lev['vapour']['total'] > 0.0 > lev['heat']['total']                       # vapour up, heat down
    assert lev['khv_m2_s'] == pytest.approx(10.0)
    totals = [r['heat_resolved'][2] + r['heat_diffusion'][2] + r['heat_subgrid'][2] for r in rows]
    assert lev['heat']['total_spread'] == pytest.approx(np.std(totals))


def test_the_quarters_run_from_sunrise_by_hour_angle_from_noon():
    assert [bl.quarter_name(h) for h in (-135.0, -45.0, 45.0, 135.0, 225.0)] == \
        ['midnight-sunrise', 'sunrise-noon', 'noon-sunset', 'sunset-midnight', 'midnight-sunrise']
