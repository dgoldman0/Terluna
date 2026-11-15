"""The spin axis's precession with the orbit's node in PlaSim (climate/gcm/precession.py): the patch builds on the design
run's chain from ExoPlaSim's untouched sources and only a run that asks for it takes it; the copy of PlaSim's orbit
holds the fixed orbit's seasons; the precession keeps the distance, turns the seasons once an eclipse year and hands
the summers at perihelion from one pole to the other in half a nodal period; and the run's setting, its first step and
its namelist entries come out as PlaSim reads them."""
from __future__ import annotations

import struct

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import levels as lv
from climate.gcm import patch_ledger as plg
from climate.gcm import precession as prc

DESIGN = dict(twilight=True, depth=True, cap=True, gust=True, similarity=True, rain_drops=True)
STEPS = 17532                                 # PlaSim's n_steps_per_year: the sidereal year in 30-minute steps
RUN = dict(eccen=0.016715, mvelp_deg=102.7, obliq_deg=1.613)          # the turned design run's planet_namelist


def test_only_a_run_that_asks_for_it_takes_the_patch():
    plain, precessing = er.patch_chain(**DESIGN), er.patch_chain(**DESIGN, precession=True)
    assert all(marker != prc.MARKER for _, marker, _ in plain)
    assert precessing == plain + prc.PATCHES
    assert ('radmod.f90', prc.MARKER) in [(n, m) for n, m, _ in er.patch_catalogue()]


@pytest.mark.skipif(not (er.UNTOUCHED_SOURCE / 'radmod.f90').exists(), reason='no untouched ExoPlaSim sources at hand')
def test_the_patch_builds_on_the_design_chain_from_the_untouched_sources():
    text = None
    for name, marker, edits in er.patch_chain(**DESIGN, precession=True):
        if name == 'radmod.f90':
            text = plg.apply(text or (er.UNTOUCHED_SOURCE / name).read_text(), edits, name)
    assert text.count(prc.MARKER) == 1
    assert 'call orb_decl(zcday, eccen, mvelpp + precang, lambm0 + precang' in text
    assert 'call gen_orb_decl(zcday, eccen, obliqr, mvelpp + precang,' in text
    assert text.count('namelist/radmod_nl/') == 2 and 'call mpbcr(precstep0)' in text
    assert all(len(line) <= 132 for _, new in prc.EDITS for line in new.splitlines())


def test_the_copy_of_the_orbit_holds_the_fixed_seasons():
    # Earth's orbit as orb_decl gives it: the equinox at day 80.5 of 365, the December solstice two days before
    # perihelion; and the run's: the south's summer at perihelion
    steps = np.arange(STEPS)
    d, f = prc.sun(steps, STEPS, 0.016715, 102.7, 23.44)
    days = steps / STEPS * 365.0
    assert abs(np.degrees(d[np.argmin(np.abs(days - 80.5))])) < 0.02
    apart = abs(days[np.argmin(d)] - days[np.argmax(f)])
    assert min(apart, 365.0 - apart) < 15.0 and f.max() == pytest.approx(1.0342, abs=5e-4)    # across the year's end
    d, f = prc.sun(steps, STEPS, **RUN)
    assert np.degrees(d[np.argmax(f)]) < -1.55 and np.degrees(d.max()) == pytest.approx(1.613, abs=1e-6)


def test_precession_keeps_the_distance_and_turns_the_seasons_once_an_eclipse_year():
    period = prc.setting(365.25636)['period_years']
    step0 = 2934936                            # the turned design run's restart at the end of its year 49
    steps = step0 + np.arange(0, 20 * STEPS, 4)
    fixed_d, fixed_f = prc.sun(steps, STEPS, **RUN)
    d, f = prc.sun(steps, STEPS, **RUN, period_years=period, step0=step0)
    assert np.allclose(f, fixed_f, rtol=0, atol=1e-14)                       # the distance keeps the orbit's year
    assert np.array_equal(d[:1], fixed_d[:1])                                # unbroken at the run's first step
    up = steps[1:][(d[:-1] < 0) & (d[1:] >= 0)]                              # the Sun crossing north over the equator
    eclipse_year = np.diff(up).mean() * 1800.0 / 86400.0                     # days of 30-minute steps
    assert eclipse_year == pytest.approx(1.0 / (1.0 / 365.25 + 1.0 / prc.node_days()), abs=0.05)
    assert eclipse_year == pytest.approx(346.62, abs=0.05)
    for year, sign in ((0.0, -1.0), (period / 2.0, 1.0), (period, -1.0)):   # whose summer the perihelion brings
        near = np.abs(steps - step0 - year * STEPS) < STEPS / 2
        at_perihelion = d[near][np.argmax(f[near])]
        assert np.sign(at_perihelion) == sign and abs(np.degrees(at_perihelion)) > 1.5


def test_the_angle_is_nil_at_the_first_step_and_turns_once_a_period():
    assert prc.angle(2934936, 2934936, 18.6, STEPS) == 0.0
    assert prc.angle(2934936 + 18.6 * STEPS, 2934936, 18.6, STEPS) == pytest.approx(2 * np.pi)
    assert not np.any(prc.angle(np.arange(5), 0, 0.0, STEPS))                # PRECYEARS 0 holds the equinox


def test_the_setting_counts_geologys_nodal_period_in_orbits():
    s = prc.setting(365.25636)
    assert s['node_days'] == 6798.38
    assert s['period_years'] == pytest.approx(6798.38 / 365.25636)
    assert 18.6 < s['period_years'] < 18.62


def test_the_first_step_and_the_namelist(tmp_path):
    path = tmp_path / 'MOST_REST.00049'
    lv.write_records(path, [('nstep', struct.pack('<i', 2934936)), ('seed', struct.pack('<4i', 1, 2, 3, 4))])
    assert prc.restart_step(path) == 2934936
    entry = dict(prc.setting(365.25636), step0=prc.restart_step(path))
    assert prc.namelist(entry) == [('PRECYEARS', repr(6798.38 / 365.25636)), ('PRECSTEP0', '2934936.0')]


def test_a_branch_from_a_precessing_run_keeps_counting_from_that_runs_first_step(tmp_path):
    # the corrected design run: counted from the target run's year-49 restart, its own year-9 restart 171,600 steps on
    import json
    model = tmp_path / 'parent' / 'model'
    model.mkdir(parents=True)
    restart = model / 'MOST_REST.00009'
    lv.write_records(restart, [('nstep', struct.pack('<i', 2934936 + 10 * 17160))])
    record = tmp_path / 'parent' / 'progress.json'
    record.write_text(json.dumps(dict(configuration=dict(precession=dict(prc.setting(365.25636), step0=2934936)))))
    assert prc.restart_step(restart) == 3106536
    assert prc.branch_step0(restart) == 2934936                  # the regression continues unbroken
    # a parent without precession, or a restart outside a run's model folder, starts the count at the restart's step
    record.write_text(json.dumps(dict(configuration=dict(model=dict(layers=15)))))
    assert prc.branch_step0(restart) == 3106536
    record.unlink()
    assert prc.branch_step0(restart) == 3106536
    loose = tmp_path / 'MOST_REST.00009'
    lv.write_records(loose, [('nstep', struct.pack('<i', 3106536))])
    assert prc.branch_step0(loose) == 3106536
    # a record that counts from after the restart belongs to another state: refused
    record.write_text(json.dumps(dict(configuration=dict(precession=dict(step0=3106537)))))
    with pytest.raises(ValueError):
        prc.branch_step0(restart)
