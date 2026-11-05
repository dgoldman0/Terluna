"""The second redesign's cloud layers only (climate/gcm/sbm_redesign2_subcloud.py, SBM_SUBCLOUD_PATCHES): with 'all'
the copy is sbm_redesign2.column itself; with 'cloud' it leaves every level below the condensation level of the lowest
level's air wholly alone, keeps the column's moist enthalpy and water, rains only from the cloud layers above their
reference vapour and heats only those below their reference temperature, sends a column with nothing left to rain
there to the shallow scheme with DESIGN2's top, and gives DESIGN2's outputs bit for bit wherever DESIGN2 changes no
level below that condensation level. The routine with closure (c)'s, the two redesigns', the heat rule's and this
patch, compiled with ExoPlaSim's flags against closure (c)'s stub of PlaSim's module, gives the copy's decisions,
tendencies and rain with nsubcl_bm 1 on soundings and over a sequence of steps, to RTOL; with nsubcl_bm 0 it prints
every digit the heat-patched routine prints, for every closure, with the heat rule off and on. The patch applies once
over the heat rule's, keeps its lines within 132 columns and its place after the heat rule's in the catalogue, and the
patched rainmod.f90 compiles with PlaSim's modules; the runner records the rule and writes its namelist only when the
switch is on, leaving every other run's configuration as it was, and refuses the rule without the second redesign or
beside the heat rule, before any patch. Nothing here writes to the shared ExoPlaSim source."""
from __future__ import annotations

import shutil
import signal
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import sbm
from climate.gcm import sbm_closure as scl
from climate.gcm import sbm_redesign2 as sr2
from climate.gcm import sbm_redesign2_heat as srh
from climate.gcm import sbm_redesign2_subcloud as srs
from climate.tests import test_sbm_closure as tsc
from climate.tests import test_sbm_redesign as tsr
from climate.tests import test_sbm_redesign2 as tsr2

SIGMA, DSIGMA, NLEV, AKAP, C = tsc.SIGMA, tsc.DSIGMA, tsc.NLEV, tsc.AKAP, tsc.C
STEP2 = tsr2.STEP2
RTOL = tsr2.RTOL
KEYS = ('dtdt', 'dqdt', 'rain')


def at_reference_below_the_condensation_level(t, q, ps):
    """The sounding with its lowest level's vapour at rhref_bm of its saturation and every level below that air's
    condensation level at the second redesign's references, its vapour rhref_bm of its own saturation and its
    temperature at least the plume's, so that DESIGN2 changes nothing there."""
    t, q = t.copy(), q.copy()
    p = ps * SIGMA
    q[-1] = sr2.RH_REF * sbm.saturation(t[-1], p[-1], C)
    lcl = sbm.condensation_level(t[-1], q[-1], p[-1], p[0], AKAP, C)
    t_ref = sr2.zbp_reference(t[-1], q[-1], SIGMA, ps, C, AKAP, lcl)
    for k in range(NLEV - 1):
        if p[k] > lcl:
            t[k] = max(t[k], t_ref[k])
            q[k] = sr2.RH_REF * sbm.saturation(t[k], p[k], C)
    return t, q, ps


def still(t, q, ps, rng):
    """One column's inputs with no tendencies so far, so that the state entering the convection is the sounding to
    its last digit."""
    x = tsc.column_inputs(t, q, ps, rng)
    x.update(dt=t.copy(), dq=q.copy(), dtdt=np.zeros(NLEV), dqdt=np.zeros(NLEV))
    return x


def cases(seed=29):
    """The second redesign's soundings as closure (c)'s tests give them, and each again set at DESIGN2's references
    below its condensation level."""
    rng = np.random.default_rng(seed)
    return ([tsc.column_inputs(t, q, ps, rng) for t, q, ps in tsr2.columns()]
            + [still(*at_reference_below_the_condensation_level(t, q, ps), rng) for t, q, ps in tsr2.columns()])


def entering(x):
    return x['dt'] + x['dtdt'] * STEP2, x['dq'] + x['dqdt'] * STEP2


def both(x):
    te, qe = entering(x)
    return (sr2.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2),
            srs.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2), te, qe)


# ------------------------------------------------------------------------------------------------------- the copy

def test_the_copy_leaves_the_sub_cloud_levels_alone_and_is_design2_where_design2_changes_nothing_there():
    seen = dict(changed=0, to_dry=0, untouched_deep=0)
    for x in cases() + [c for _, cols in tsr2.sequence_steps(nsteps=6) for c in cols]:
        a, b, te, qe = both(x)
        r = srs.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2, rule='all')
        assert r['decision'] == a['decision']
        for key in KEYS:
            assert np.array_equal(r[key], a[key])                       # 'all' is DESIGN2 itself
        p = x['dp'] * SIGMA
        below = p > a['lcl']
        if not (np.any(a['dtdt'][below]) or np.any(a['dqdt'][below])):  # DESIGN2 changes nothing below the lcl
            assert b['decision'] == a['decision'] and b['shallow_top'] == a['shallow_top']
            for key in KEYS:
                assert np.array_equal(b[key], a[key]), key              # bit for bit
            seen['untouched_deep'] += a['decision'] == 'deep'
            continue
        assert a['decision'] == 'deep'
        seen['changed'] += 1
        assert not (np.any(b['dtdt'][below]) or np.any(b['dqdt'][below]) or np.any(b['rain'][below]))
        if b['decision'] == 'dry':                                      # nothing left to rain in the cloud layers
            seen['to_dry'] += 1
            assert b['excess'] <= 0.0 or b['pt'] <= 0.0
            assert b['shallow_top'] == next((k for k in range(a['top'], NLEV) if p[k] >= a['base'] - sr2.DEEP_PA),
                                            NLEV - 1)
            assert not (np.any(b['dtdt']) or np.any(b['dqdt']) or np.any(b['rain']))
            continue
        assert b['decision'] == 'deep' and b['top'] == a['top']
        m = x['dp'] * DSIGMA / C['ga']
        cp = C['acpd'] * (1.0 + C['ADV'] * qe)
        rain = b['rain'].sum()
        assert rain > 0.0
        assert abs(-(b['dqdt'] * m).sum() - rain) <= 1e-12 * rain                       # water kept
        assert abs((cp * b['dtdt'] * m).sum() - C['ALV'] * rain) <= 1e-12 * C['ALV'] * rain   # moist enthalpy kept
        cloud = (np.arange(NLEV) >= b['top']) & ~below
        short = np.where(cloud, np.maximum(b['t_ref'] - te, 0.0), 0.0)
        assert np.all((b['dtdt'] > 0.0) == (short > 0.0))               # heat only into cloud layers below reference
        assert np.all((b['rain'] > 0.0) == (cloud & (qe > b['q_ref'])))  # rain only from cloud layers above it
    assert all(v > 0 for v in seen.values()), seen                      # every branch ran


def test_the_copy_refuses_an_unknown_rule_and_records_only_the_cloud_layers():
    x = cases()[0]
    te, qe = entering(x)
    with pytest.raises(ValueError):
        srs.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2, rule='export')
    with pytest.raises(ValueError):
        srs.record('all')
    assert srs.record()['constants'] == {'nsubcl_bm': 1}


# ------------------------------------------------------------------------------------------------------- Fortran

def sources():
    """mksbm's sources: the heat rule's (the second redesign's with the heat patch over it), and with this patch over
    that."""
    s = tsr2.sources()
    heat = tsc.apply(s['redesign2'], srh.SBM_HEAT_EDITS, srh.HEAT_MARKER)
    return dict(heat=heat, subcloud=tsc.apply(heat, srs.SBM_SUBCLOUD_EDITS, srs.SUBCLOUD_MARKER))


def build(folder, code, nhor, switch, heat):
    folder.mkdir(parents=True, exist_ok=True)
    declarations = (scl._declarations() + sr2.declarations()
                    + srh.declarations().replace('nheat_bm = 0', f'nheat_bm = {heat}')
                    + srs.declarations().replace('nsubcl_bm = 0', f'nsubcl_bm = {switch}'))
    files = {'stub.f90': tsc.STUB.format(nhor=nhor, nlev=NLEV, declarations=declarations),
             'gustcol.f90': tsc.gustcol_source(), 'mksbm.f90': code, 'harness.f90': tsr.HARNESS}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


@pytest.fixture(scope='module')
def executables(tmp_path_factory):
    built, texts = {}, {}

    def get(which, nhor, switch=0, heat=0):
        if (which, nhor, switch, heat) not in built:
            if not texts:
                texts.update(sources())
            built[which, nhor, switch, heat] = build(tmp_path_factory.mktemp(f'{which}{nhor}s{switch}h{heat}'),
                                                     tsc.routine(texts[which], 'mksbm'), nhor, switch, heat)
        return built[which, nhor, switch, heat]
    return get


def python_step(x):
    te, qe = entering(x)
    return srs.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2), te, qe


def differs(x):
    """Whether the cloud layers' column differs from DESIGN2's on these inputs."""
    a, b, _, _ = both(x)
    return a['decision'] != b['decision'] or not all(np.array_equal(a[k], b[k]) for k in KEYS)


@tsc.FORTRAN
def test_the_fortran_cloud_layers_give_the_copys_decisions_tendencies_and_rain(executables):
    xs = cases()
    exe = executables('subcloud', len(xs), 1)
    printed, _ = tsc.run_harness(exe, tsr2.harness_text('redesign2', [(100, xs)]))
    worst, changed = {}, 0
    for j, (x, f) in enumerate(zip(xs, tsr.parsed(printed, 1, len(xs))[0])):
        py, te, qe = python_step(x)
        tsr2.assert_step_matches(f, py, x, qe, j, worst)
        changed += differs(x)
    assert changed > 0                                              # the rule's own branch ran
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


@tsc.FORTRAN
def test_the_fortran_cloud_layers_follow_the_copy_step_for_step(executables):
    steps = tsr2.sequence_steps()
    ncols = len(steps[0][1])
    exe = executables('subcloud', ncols, 1)
    printed, _ = tsc.run_harness(exe, tsr2.harness_text('redesign2', steps))
    worst, changed = {}, 0
    for (nstep, cols), frow in zip(steps, tsr.parsed(printed, len(steps), ncols)):
        for j, (x, f) in enumerate(zip(cols, frow)):
            py, te, qe = python_step(x)
            tsr2.assert_step_matches(f, py, x, qe, (nstep, j), worst)
            changed += differs(x)
    assert changed > 0
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


@tsc.FORTRAN
def test_with_the_switch_off_the_patched_routine_prints_every_digit_the_heat_patched_routine_prints(executables):
    rng = np.random.default_rng(37)
    xs = [tsc.column_inputs(t, q, ps, rng) for t, q, ps in tsr2.columns()]
    route = tsc.route_l_steps(ncols=tsr2.NHOR, nsteps=60)
    twice = [(7, xs), (8, xs[::-1])]
    runs = [(setting, steps, 0) for setting, steps in (('design', twice), ('plain', twice), ('closure', twice),
                                                         ('route_l', route), ('redesign', twice), ('redesign2', twice))]
    for setting, steps, heat in runs + [('redesign2', twice, 1)]:     # and with the heat rule on, where it acts
        text = tsr2.harness_text(setting, steps)
        mine = tsc.run_harness(executables('subcloud', tsr2.NHOR, 0, heat), text)[0]
        theirs = tsc.run_harness(executables('heat', tsr2.NHOR, 0, heat), text)[0]
        assert mine == theirs, (setting, heat)                         # every digit


def test_the_patch_applies_once_over_the_heat_rules_and_holds_its_place():
    marker = srs.SUBCLOUD_MARKER
    assert not any(marker in old for old, _ in srs.SBM_SUBCLOUD_EDITS)
    assert all(marker in new for _, new in srs.SBM_SUBCLOUD_EDITS)
    for _, new in srs.SBM_SUBCLOUD_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines()), new[:60]
    texts = []
    if tsc.shared_rainmod() is not None:
        s = sources()
        assert tsc.undo(s['subcloud'], srs.SBM_SUBCLOUD_EDITS, marker) == s['heat']
        texts.append(s['subcloud'])
    if (er.UNTOUCHED_SOURCE / 'rainmod.f90').exists():
        texts.append(tsc.fresh_rainmod(dict(tsc.FRESH_SWITCHES, redesign2=True, subcloud=True)))
        assert marker not in tsc.fresh_rainmod(dict(tsc.FRESH_SWITCHES, redesign2=True, heat=True))
    if not texts:
        pytest.skip('no ExoPlaSim source at hand')
    for text in texts:
        assert tsc.apply(text, srs.SBM_SUBCLOUD_EDITS, marker) == text
        assert srh.HEAT_MARKER in text and sr2.REDESIGN2_MARKER in text
        assert text.count('call mpbci(nsubcl_bm)') == 1 and text.count('namelist/rainmod_nl/nsubcl_bm') == 1
        assert text.count('if(nsubcl_bm == 1') == 1
    chain = er.patch_chain(**tsc.DESIGN_SWITCHES, redesign2=True, subcloud=True)
    assert chain[-3:] == er.SBM_REDESIGN2_PATCHES + er.SBM_HEAT_PATCHES + er.SBM_SUBCLOUD_PATCHES
    assert er.SBM_SUBCLOUD_PATCHES[0] not in er.patch_chain(**tsc.DESIGN_SWITCHES, redesign2=True, heat=True)
    catalogue = [(n, m) for n, m, _ in er.patch_catalogue()]
    assert catalogue.index(('rainmod.f90', marker)) > catalogue.index(('rainmod.f90', srh.HEAT_MARKER))


@tsc.FORTRAN
def test_the_subcloud_patched_rainmod_compiles_with_plasims_modules(tmp_path):
    path = tsc.shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(path.parent / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'rainmod.f90').write_text(sources()['subcloud'])
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------- runner

GOOD = dict(deep_hpa=200.0, shallow_cap=True, closure='redesign2')


def test_the_runner_records_the_rule_only_when_the_switch_is_on():
    plain = er.convection_setting(**GOOD)
    for off in (dict(GOOD, subcloud=None), dict(GOOD, subcloud='all')):
        assert er.convection_setting(**off) == plain                # every other run's configuration as it was
    assert 'NSUBCL_BM' not in dict(er.convection_namelist(plain))
    on = er.convection_setting(**dict(GOOD, subcloud='cloud'))
    assert on['subcloud'] == srs.record() and on['subcloud']['switch'] == 'nsubcl_bm 1'
    assert {k: v for k, v in on.items() if k != 'subcloud'} == plain
    nl = dict(er.convection_namelist(on))
    assert nl['NSUBCL_BM'] == '1' and nl['NCLOS_BM'] == '3' and 'NHEAT_BM' not in nl
    for bad in (dict(GOOD, subcloud='export'), dict(GOOD, subcloud='cloud', heat='spread'),
                dict(deep_hpa=200.0, shallow_cap=True, closure='redesign', subcloud='cloud'),
                dict(deep_hpa=200.0, shallow_cap=True, closure='c', subcloud='cloud'), dict(subcloud='cloud')):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    assert er.uses_subcloud_patch(dict(GOOD, subcloud='cloud')) and not er.uses_subcloud_patch(GOOD)
    assert not er.uses_subcloud_patch(dict(GOOD, subcloud='all')) and not er.uses_subcloud_patch(None)


def test_the_configuration_refuses_the_rule_without_the_second_redesign_before_any_patch(monkeypatch):
    def reached(*a, **k):
        raise AssertionError('ensure_patched reached: the refusal must come before the patch')
    monkeypatch.setattr(er, 'ensure_patched', reached)
    for convection in (dict(deep_hpa=200.0, shallow_cap=True, closure='redesign', subcloud='cloud'),
                       dict(GOOD, subcloud='cloud', heat='spread')):
        with pytest.raises(ValueError, match='second redesign'):
            er.configuration('A28R', 4, convection=convection, levels='l15', surface_gust='beljaars')


def test_ensure_patched_refuses_the_rule_without_the_second_redesign(tmp_path, monkeypatch):
    import sys
    import types
    package = tmp_path / 'exoplasim'
    (package / 'plasim' / 'src').mkdir(parents=True)
    monkeypatch.setitem(sys.modules, 'exoplasim', types.SimpleNamespace(__file__=str(package / '__init__.py')))
    with pytest.raises(RuntimeError, match='second redesign'):
        er.ensure_patched(subcloud=True)
    with pytest.raises(RuntimeError, match='second redesign'):
        er.ensure_patched(depth=True, cap=True, subcloud=True)
    assert not list((package / 'plasim' / 'src').iterdir())            # nothing was written


def test_the_runners_flags_take_the_rule_and_refuse_it_without_the_second_redesign(monkeypatch):
    seen = {}
    monkeypatch.setattr(er, 'run', lambda *a, **k: seen.update(args=a) or 0)
    base = ['A28R', '--convection', 'sbm', '--sbm-deep-hpa', '200', '--levels', 'l15', '--sbm-shallow-cap',
            '--surface-gust', 'beljaars_similarity', '--rain-drops', 'share']
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        assert er.main(base + ['--sbm-closure', 'redesign2', '--sbm-subcloud', 'cloud']) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['subcloud'] == 'cloud'
        assert er.convection_setting(**convection)['subcloud']['rule'] == 'cloud'
        assert er.main(base + ['--sbm-closure', 'redesign2']) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['subcloud'] is None and 'subcloud' not in er.convection_setting(**convection)
        for argv in (base + ['--sbm-closure', 'redesign', '--sbm-subcloud', 'cloud'],
                     base + ['--sbm-subcloud', 'cloud'],
                     base + ['--sbm-closure', 'redesign2', '--sbm-subcloud', 'export'],
                     base + ['--sbm-closure', 'redesign2', '--sbm-subcloud', 'cloud', '--sbm-heat', 'spread']):
            with pytest.raises(SystemExit):
                er.main(argv)
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
