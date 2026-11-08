"""The second redesign's heat beyond its reference (climate/gcm/sbm_redesign2_heat.py, SBM_HEAT_PATCHES): the copy
keeps DESIGN2's decisions on every column and, where DESIGN2's limit on the rain does not bind, every output bit for
bit; where the limit binds it rains the vapour's relaxation unheld, warms the layers below their reference at the
limit's own rate and spreads the rest of the latent heat over the cloud layers by one increment, keeping the column's
moist enthalpy and water and cooling no layer. The routine with closure (c)'s, the two redesigns' and this patch,
compiled with ExoPlaSim's flags against closure (c)'s stub of PlaSim's module, gives the copy's decisions, tendencies
and rain with nheat_bm 1 on soundings and over a sequence of steps, to RTOL; with nheat_bm 0 it prints every digit the
second redesign's routine prints, for every closure. The patch applies once over the second redesign's, keeps its
lines within 132 columns and its place last in the catalogue, and the patched rainmod.f90 compiles with PlaSim's
modules; the runner records the rule and writes its namelist only when the switch is on, leaving every other run's
configuration as it was, and refuses the rule without the second redesign, before any patch. Nothing here writes to
the shared ExoPlaSim source."""
from __future__ import annotations

import shutil
import signal
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import sbm_closure as scl
from climate.gcm import sbm_redesign2 as sr2
from climate.gcm import sbm_redesign2_heat as srh
from climate.tests import test_sbm_closure as tsc
from climate.tests import test_sbm_redesign as tsr
from climate.tests import test_sbm_redesign2 as tsr2

SIGMA, DSIGMA, NLEV, AKAP, C = tsc.SIGMA, tsc.DSIGMA, tsc.NLEV, tsc.AKAP, tsc.C
STEP2 = tsr2.STEP2
RTOL = tsr2.RTOL


def cases(seed=29):
    rng = np.random.default_rng(seed)
    return [tsc.column_inputs(t, q, ps, rng) for t, q, ps in tsr2.columns()]


def entering(x):
    return x['dt'] + x['dtdt'] * STEP2, x['dq'] + x['dqdt'] * STEP2


# ------------------------------------------------------------------------------------------------------- the copy

def test_the_copy_keeps_design2_where_its_limit_does_not_bind_and_its_heat_elsewhere():
    seen = dict(limited=0, unbound=0)
    for x in cases() + [c for _, cols in tsr2.sequence_steps(nsteps=6) for c in cols]:
        te, qe = entering(x)
        a = sr2.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2)
        b = srh.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2)
        r = srh.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2, rule='reference')
        for key in ('dtdt', 'dqdt', 'rain'):
            assert np.array_equal(r[key], a[key])                   # 'reference' is DESIGN2 itself
        assert b['decision'] == a['decision']                       # the decisions are DESIGN2's
        if a['decision'] != 'deep' or not a['limited']:
            seen['unbound'] += a['decision'] == 'deep'
            for key in ('dtdt', 'dqdt', 'rain'):
                assert np.array_equal(b[key], a[key]), key          # bit for bit where the limit does not bind
            assert b['spread'] == 0.0
            continue
        seen['limited'] += 1
        m = x['dp'] * DSIGMA / C['ga']
        cp = C['acpd'] * (1.0 + C['ADV'] * qe)
        rain = b['rain'].sum()
        moist = min(a['deficit'], sr2.F_MOIST * a['excess'])
        assert np.isclose(rain, (a['excess'] - moist) / sr2.TAU_S, rtol=1e-12)     # the rain is not held
        assert rain > a['rain'].sum()
        assert abs(-(b['dqdt'] * m).sum() - rain) <= 1e-12 * rain                # water kept
        assert abs((cp * b['dtdt'] * m).sum() - C['ALV'] * rain) <= 1e-12 * C['ALV'] * rain   # moist enthalpy kept
        assert b['dtdt'].min() >= 0.0 and b['spread'] > 0.0          # no layer cools; the spread warms
        p = x['dp'] * SIGMA
        cloud = (np.arange(NLEV) >= a['top']) & (p <= a['lcl'])
        toward = np.array([max(a['t_ref'][k] - te[k], 0.0) / sr2.TAU_S * (sr2.TAU_S / STEP2) if k >= a['top'] else 0.0
                           for k in range(NLEV)])
        assert np.allclose(b['dtdt'], toward + np.where(cloud, b['spread'], 0.0), rtol=1e-12, atol=1e-20)
        assert np.all(te + STEP2 * toward <= np.maximum(a['t_ref'], te) + 1e-9)   # the toward part stops there
    assert seen['limited'] > 0 and seen['unbound'] > 0, seen         # both branches ran


def test_the_copy_refuses_an_unknown_rule():
    x = cases()[0]
    te, qe = entering(x)
    with pytest.raises(ValueError):
        srh.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2, rule='frierson')


# ------------------------------------------------------------------------------------------------------- Fortran

def sources():
    """mksbm's sources: the second redesign's (test_sbm_redesign2.sources), and with this patch over it."""
    s = tsr2.sources()
    return dict(redesign2=s['redesign2'], heat=tsc.apply(s['redesign2'], srh.SBM_HEAT_EDITS, srh.HEAT_MARKER))


def build(folder, code, nhor, switch):
    folder.mkdir(parents=True, exist_ok=True)
    declarations = scl._declarations() + sr2.declarations() + srh.declarations().replace('nheat_bm = 0',
                                                                                         f'nheat_bm = {switch}')
    files = {'stub.f90': tsc.STUB.format(nhor=nhor, nlev=NLEV, declarations=declarations),
             'gustcol.f90': tsc.gustcol_source(), 'mksbm.f90': code, 'harness.f90': tsr.HARNESS}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


@pytest.fixture(scope='module')
def executables(tmp_path_factory):
    built, texts = {}, {}

    def get(which, nhor, switch=0):
        if (which, nhor, switch) not in built:
            if not texts:
                texts.update(sources())
            built[which, nhor, switch] = build(tmp_path_factory.mktemp(f'{which}{nhor}s{switch}'),
                                               tsc.routine(texts[which], 'mksbm'), nhor, switch)
        return built[which, nhor, switch]
    return get


def python_step(x):
    te, qe = entering(x)
    return srh.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2), te, qe


@tsc.FORTRAN
def test_the_fortran_heat_rule_gives_the_copys_decisions_tendencies_and_rain(executables):
    xs = cases()
    exe = executables('heat', tsr2.NHOR, 1)
    printed, _ = tsc.run_harness(exe, tsr2.harness_text('redesign2', [(100, xs)]))
    worst, spread = {}, 0
    for j, (x, f) in enumerate(zip(xs, tsr.parsed(printed, 1, len(xs))[0])):
        py, te, qe = python_step(x)
        tsr2.assert_step_matches(f, py, x, qe, j, worst)
        spread += py['spread'] > 0.0
    assert spread > 0                                               # the spread branch ran
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


@tsc.FORTRAN
def test_the_fortran_heat_rule_follows_the_copy_step_for_step(executables):
    steps = tsr2.sequence_steps()
    ncols = len(steps[0][1])
    exe = executables('heat', ncols, 1)
    printed, _ = tsc.run_harness(exe, tsr2.harness_text('redesign2', steps))
    worst, spread = {}, 0
    for (nstep, cols), frow in zip(steps, tsr.parsed(printed, len(steps), ncols)):
        for j, (x, f) in enumerate(zip(cols, frow)):
            py, te, qe = python_step(x)
            tsr2.assert_step_matches(f, py, x, qe, (nstep, j), worst)
            spread += py['spread'] > 0.0
    assert spread > 0
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


@tsc.FORTRAN
def test_with_the_switch_off_the_patched_routine_prints_every_digit_the_second_redesigns_prints(executables):
    rng = np.random.default_rng(37)
    xs = [tsc.column_inputs(t, q, ps, rng) for t, q, ps in tsr2.columns()]
    route = tsc.route_l_steps(ncols=tsr2.NHOR, nsteps=60)
    twice = [(7, xs), (8, xs[::-1])]
    for setting, steps in (('design', twice), ('plain', twice), ('closure', twice), ('route_l', route),
                           ('redesign', twice), ('redesign2', twice)):
        text = tsr2.harness_text(setting, steps)
        mine = tsc.run_harness(executables('heat', tsr2.NHOR, 0), text)[0]
        theirs = tsc.run_harness(executables('redesign2', tsr2.NHOR, 0), text)[0]
        assert mine == theirs, setting                               # every digit


def test_the_patch_applies_once_over_the_second_redesigns_and_holds_its_place():
    marker = srh.HEAT_MARKER
    assert not any(marker in old for old, _ in srh.SBM_HEAT_EDITS)
    assert all(marker in new for _, new in srh.SBM_HEAT_EDITS[:3])
    for _, new in srh.SBM_HEAT_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines()), new[:60]
    texts = []
    if tsc.shared_rainmod() is not None:
        s = sources()
        assert tsc.undo(s['heat'], srh.SBM_HEAT_EDITS, marker) == s['redesign2']
        texts.append(s['heat'])
    if (er.UNTOUCHED_SOURCE / 'rainmod.f90').exists():
        texts.append(tsc.fresh_rainmod(dict(tsc.FRESH_SWITCHES, redesign2=True, heat=True)))
        assert marker not in tsc.fresh_rainmod(dict(tsc.FRESH_SWITCHES, redesign2=True))
    if not texts:
        pytest.skip('no ExoPlaSim source at hand')
    for text in texts:
        assert tsc.apply(text, srh.SBM_HEAT_EDITS, marker) == text
        assert sr2.REDESIGN2_MARKER in text
        assert text.count('call mpbci(nheat_bm)') == 1 and text.count('namelist/rainmod_nl/nheat_bm') == 1
    chain = er.patch_chain(**tsc.DESIGN_SWITCHES, redesign2=True, heat=True)
    assert chain[-2:] == er.SBM_REDESIGN2_PATCHES + er.SBM_HEAT_PATCHES
    assert er.SBM_HEAT_PATCHES[0] not in er.patch_chain(**tsc.DESIGN_SWITCHES, redesign2=True)
    catalogue = [(n, m) for n, m, _ in er.patch_catalogue()]
    assert catalogue.index(('rainmod.f90', marker)) > catalogue.index(('rainmod.f90', sr2.REDESIGN2_MARKER))


@tsc.FORTRAN
def test_the_heat_patched_rainmod_compiles_with_plasims_modules(tmp_path):
    path = tsc.shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(path.parent / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'rainmod.f90').write_text(sources()['heat'])
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------- runner

GOOD = dict(deep_hpa=200.0, shallow_cap=True, closure='redesign2')


def test_the_runner_records_the_rule_only_when_the_switch_is_on():
    plain = er.convection_setting(**GOOD)
    for off in (dict(GOOD, heat=None), dict(GOOD, heat='reference')):
        assert er.convection_setting(**off) == plain                # every other run's configuration as it was
    assert 'NHEAT_BM' not in dict(er.convection_namelist(plain))
    on = er.convection_setting(**dict(GOOD, heat='spread'))
    assert on['heat'] == srh.record() and on['heat']['switch'] == 'nheat_bm 1'
    assert {k: v for k, v in on.items() if k != 'heat'} == plain
    nl = dict(er.convection_namelist(on))
    assert nl['NHEAT_BM'] == '1' and nl['NCLOS_BM'] == '3'
    for bad in (dict(GOOD, heat='frierson'), dict(deep_hpa=200.0, shallow_cap=True, closure='redesign', heat='spread'),
                dict(deep_hpa=200.0, shallow_cap=True, closure='c', heat='spread'), dict(heat='spread')):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    assert er.uses_heat_patch(dict(GOOD, heat='spread')) and not er.uses_heat_patch(GOOD)
    assert not er.uses_heat_patch(None)


def test_the_configuration_refuses_the_rule_without_the_second_redesign_before_any_patch(monkeypatch):
    def reached(*a, **k):
        raise AssertionError('ensure_patched reached: the refusal must come before the patch')
    monkeypatch.setattr(er, 'ensure_patched', reached)
    with pytest.raises(ValueError, match='second redesign'):
        er.configuration('A28R', 4, convection=dict(deep_hpa=200.0, shallow_cap=True, closure='redesign',
                                                    heat='spread'), levels='l15', surface_gust='beljaars')


def test_ensure_patched_refuses_the_rule_without_the_second_redesign(tmp_path, monkeypatch):
    import sys
    import types
    package = tmp_path / 'exoplasim'
    (package / 'plasim' / 'src').mkdir(parents=True)
    monkeypatch.setitem(sys.modules, 'exoplasim', types.SimpleNamespace(__file__=str(package / '__init__.py')))
    with pytest.raises(RuntimeError, match='second redesign'):
        er.ensure_patched(heat=True)
    with pytest.raises(RuntimeError, match='second redesign'):
        er.ensure_patched(depth=True, cap=True, heat=True)
    assert not list((package / 'plasim' / 'src').iterdir())            # nothing was written


def test_the_runners_flags_take_the_rule_and_refuse_it_without_the_second_redesign(monkeypatch):
    seen = {}
    monkeypatch.setattr(er, 'run', lambda *a, **k: seen.update(args=a) or 0)
    base = ['A28R', '--convection', 'sbm', '--sbm-deep-hpa', '200', '--levels', 'l15', '--sbm-shallow-cap',
            '--surface-gust', 'beljaars_similarity', '--rain-drops', 'share']
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        assert er.main(base + ['--sbm-closure', 'redesign2', '--sbm-heat', 'spread']) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['heat'] == 'spread' and er.convection_setting(**convection)['heat']['rule'] == 'spread'
        assert er.main(base + ['--sbm-closure', 'redesign2']) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['heat'] is None and 'heat' not in er.convection_setting(**convection)
        for argv in (base + ['--sbm-closure', 'redesign', '--sbm-heat', 'spread'], base + ['--sbm-heat', 'spread'],
                     base + ['--sbm-closure', 'redesign2', '--sbm-heat', 'frierson']):
            with pytest.raises(SystemExit):
                er.main(argv)
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
