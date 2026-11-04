"""The redesigned closure of the simplified Betts-Miller scheme (climate/gcm/sbm_redesign.py, SBM_REDESIGN_PATCHES): the
copy keeps the column's moist enthalpy and water for every decision, forms no rain where it moistens, warms a layer only
toward a reference above it, takes the design scheme's trigger, depth test and shallow top, and relaxes toward Romps's
humidity only in the buoyant cloud. The routine with closure (c)'s patch and this one, compiled with ExoPlaSim's flags
against closure (c)'s stub of PlaSim's module and the gust patch's own gustcol, gives the Python copy's decisions,
tendencies, rain, shallow tops and convective cloud layers on soundings and over a sequence of steps that moves each
column through the decisions, to 1e-10, without array temporaries under -fcheck=all, and keeps the moist enthalpy and
water there too; with nclos_bm 0 or 1 it prints every digit the routine without this patch prints, and with nclos_bm 0
every digit the routine without either patch prints. The patch applies once to the shared rainmod.f90 and to the
untouched source with a redesign run's chain, after closure (c)'s, keeps its lines within 132 columns, and the patched
rainmod.f90 compiles with PlaSim's modules; the runner records the redesign, writes its namelist and refuses what it
cannot run with, before any patch."""
from __future__ import annotations

import signal
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import sbm
from climate.gcm import sbm_closure as scl
from climate.gcm import sbm_redesign as srd
from climate.tests import test_sbm_closure as tsc      # closure (c)'s soundings, stub, harness and sources

SIGMA, DSIGMA, HALF, NLEV, AKAP, C = tsc.SIGMA, tsc.DSIGMA, tsc.HALF, tsc.NLEV, tsc.AKAP, tsc.C
STEP2 = 2.0 * scl.STEP_S
# The Fortran against the copy: the same operations in the same order, apart from the rain's units (dsigma, then
# dp / ga in the harness, where the copy takes the layer's mass) and the compiler's evaluation of the same
# expressions, so the two agree to a few units in the last place; 1e-10 leaves room for that and catches any change
# of formula, which shows at 1e-6 or more.
RTOL = 1e-10
# (tau_bm, rh_bm, dpdeep_bm, nshcap_bm, nclos_bm, nroutel_bm), closure (c)'s settings with the redesign's beside them
SETTINGS = dict(tsc.SETTINGS, redesign=(srd.TAU_S, 0.7, srd.DEEP_PA, 1, 2, 0))
VARIANTS = ((1.0, 0.0), (0.9, 0.0), (1.05, 1.0), (1.1, 2.0), (0.4, 0.0), (0.7, 0.0), (1.0, -3.0), (0.6, -2.0),
            (1.2, 3.0))                                     # (vapour factor, warming of the lowest levels, K)
# closure (c)'s harness printing each column's count of convective cloud layers after its decision
HARNESS = tsc.HARNESS.replace("write(*,'(4i6,es25.16e3)') kshal(jh), ktopo(jh), nsbmroute(jh), nsbmsrc(jh), zlast",
                              "write(*,'(4i6,es25.16e3,i6)') kshal(jh), ktopo(jh), nsbmroute(jh), nsbmsrc(jh), zlast, "
                              "icctot(jh)")
assert 'icctot(jh)' in HARNESS
FORTRAN = tsc.FORTRAN


def columns():
    """The soundings of closure (c)'s tests, each moistened or dried and its lowest levels warmed or cooled."""
    return [tsc.sounding(name, moist, warm) for name in tsc.SOUNDINGS for moist, warm in VARIANTS]


NHOR = len(VARIANTS) * len(tsc.SOUNDINGS)                  # every harness here runs this many columns


# ------------------------------------------------------------------------------------------------------- the copy

def test_the_copy_keeps_moist_enthalpy_and_water_and_moistens_without_rain():
    seen = dict.fromkeys(srd.DECISIONS, 0)
    clipped = warmed = 0
    for t, q, ps in columns():
        o = srd.column(t, q, SIGMA, DSIGMA, ps, C, AKAP)
        seen[o['decision']] += 1
        b = srd.budgets(o, q, SIGMA, DSIGMA, ps, C)
        flux = srd.enthalpy_flux(o, q, SIGMA, DSIGMA, ps, C)
        scale = max(1.0, C['ALV'] * b['rain'], abs(o['pt']))
        assert abs(b['enthalpy']) <= 1e-9 * scale and abs(flux[-1]) <= 1e-9 * scale
        if o['decision'] in ('none', 'thin'):
            assert not np.any(o['dtdt']) and not np.any(o['dqdt']) and not np.any(o['rain'])
            continue
        relaxed = np.arange(NLEV) >= o['top']
        assert not np.any(o['dtdt'][~relaxed]) and not np.any(o['dqdt'][~relaxed])
        warm = o['dtdt'][relaxed] - o['shift']                         # each layer's heating toward its reference
        assert np.all(warm >= -1e-12 * np.max(np.abs(o['dtdt'])))
        clipped += int(np.sum(o['t_ref'][relaxed] < t[relaxed]))
        warmed += int(np.sum(o['t_ref'][relaxed] > t[relaxed]))
        m = ps * DSIGMA / C['ga']
        held = float(np.sum(q * m)) / srd.TAU_S                     # the column's vapour over the relaxation time
        if o['decision'] == 'deep':
            assert b['rain'] > 0.0 and abs(b['vapour'] + b['rain']) <= 1e-12 * b['rain']
            assert abs(b['heating'] - C['ALV'] * b['rain']) <= 1e-9 * C['ALV'] * b['rain']
            assert np.all((o['rain'] > 0.0) == (o['dqdt'] < 0.0))
        else:
            assert not np.any(o['rain']) and 0.0 < o['factor'] <= 1.0
            assert abs(b['vapour']) <= 1e-12 * held and abs(b['heating']) <= 1e-9 * max(1.0, abs(o['pt']))
    assert all(seen.values()), seen                                  # every decision, on these soundings
    assert clipped and warmed                                        # the reference both below and above the air


def test_the_copy_takes_the_design_schemes_trigger_and_romps_humidity_only_in_the_buoyant_cloud():
    for t, q, ps in columns():
        o = srd.column(t, q, SIGMA, DSIGMA, ps, C, AKAP)
        d = sbm.column(t, q, SIGMA, DSIGMA, ps, tau=srd.TAU_S, rh=0.7, step=srd.STEP2_S, akap=AKAP, c=C,
                       deep_pa=srd.DEEP_PA, shallow_cap=True)
        assert np.array_equal(o['trigger_t'], d['reference_t'])      # the design scheme's parcel, bit for bit
        top = d.get('parcel_top', d['top'])
        assert o['top'] == top
        if top is None:
            continue
        p = ps * SIGMA
        thin = d['base'] - p[top] <= srd.DEEP_PA
        assert (o['decision'] == 'thin') == thin and (not thin or o['shallow_top'] == top)
        if thin:
            continue
        assert np.array_equal(o['t_ref'], sbm.diluted_parcel(t, q, SIGMA, ps, srd.DILUTION, AKAP, C,
                                                             max_mix=srd.MAX_MIX))
        k = np.arange(NLEV)
        cloud = (k >= o['top']) & (k < NLEV - 1) & (p <= o['base']) & (k <= o['lfc'])
        assert np.all(np.isfinite(o['rh_ref'][cloud])) and not np.any(np.isfinite(o['rh_ref'][~cloud]))
        below = (k >= o['top']) & ~cloud
        assert np.array_equal(o['q_ref'][below], q[below])            # the air's own vapour beneath the cloud
        r = scl.romps_reference(t, q, SIGMA, ps, o['t_ref'], o['top'], C)
        assert np.array_equal(o['rh_ref'][cloud], r['rh'][cloud])


# ------------------------------------------------------------------------------------------------------- Fortran

def sources():
    """mksbm's source, the shared rainmod.f90 without closure (c)'s patch and this one, with closure (c)'s, and with
    both."""
    path = tsc.shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    base = tsc.undo(tsc.without_redesign(path.read_text()), scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER)
    closure = tsc.apply(base, scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER)
    both = tsc.apply(closure, srd.SBM_REDESIGN_EDITS, srd.REDESIGN_MARKER)
    return base, closure, both


def build(folder, code, nhor):
    folder.mkdir(parents=True, exist_ok=True)
    files = {'stub.f90': tsc.STUB.format(nhor=nhor, nlev=NLEV, declarations=scl._declarations()),
             'gustcol.f90': tsc.gustcol_source(), 'mksbm.f90': code, 'harness.f90': HARNESS}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def harness_text(setting, steps):
    """The harness's input as closure (c)'s tests write it, with the redesign's setting beside theirs."""
    tau, rh, deep, cap, nclos, nroutel = SETTINGS[setting]
    lines = [tsc.num(np.r_[SIGMA, DSIGMA, HALF]),
             tsc.num([AKAP, C['ga'], C['acpd'], C['ADV'], C['ALV'], C['ALS'], C['TMELT'], C['rdbrv'], C['ra1'],
                      C['ra2'], C['ra4'], C['gascon'], scl.STEP_S]),
             f'{tau!r} {rh!r} {deep!r} {cap} {nclos} {nroutel} {len(steps)}']
    for nstep, cols in steps:
        lines.append(str(nstep))
        for x in cols:
            lines.append(tsc.num([x['dp'], x['dshfl'], x['devap'], x['dust3']]))
            lines += [tsc.num(x[k]) for k in ('dt', 'dq', 'dtdt', 'dqdt', 'du', 'dv', 'dudt', 'dvdt', 'swr', 'lwr',
                                              'dqt', 'dtt')]
    return '\n'.join(lines) + '\n'


def parsed(printed, nsteps, ncols):
    out = tsc.parsed(printed, nsteps, ncols)
    rows = printed.strip('\n').split('\n')[4::5]
    for s, cols in enumerate(out):
        for j, col in enumerate(cols):
            col['layers'] = int(rows[s * ncols + j].split()[5])
    return out


@pytest.fixture(scope='module')
def executables(tmp_path_factory):
    """The harness built once on each source, as it is first needed: 'base' without either patch, 'closure' with
    closure (c)'s, 'both' with this one over it, for nhor columns."""
    built = {}

    def get(which, nhor):
        if (which, nhor) not in built:
            code = tsc.routine(dict(zip(('base', 'closure', 'both'), sources()))[which], 'mksbm')
            built[which, nhor] = build(tmp_path_factory.mktemp(f'{which}{nhor}'), code, nhor)
        return built[which, nhor]
    return get


def python_step(x):
    """The copy on one column's inputs as the routine takes them: the state entering the convection."""
    te, qe = x['dt'] + x['dtdt'] * STEP2, x['dq'] + x['dqdt'] * STEP2
    return srd.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2), te, qe


def assert_step_matches(f, py, x, te, qe, label, worst):
    """The routine's outputs against the copy's, with the moist enthalpy and water of its tendencies; worst gathers
    the largest relative difference of the convection's own tendencies and rain."""
    assert f['route'] == srd.DECISION_CODES[py['decision']], label
    assert np.allclose(f['dtdt'], x['dtdt'] + py['dtdt'], rtol=RTOL, atol=1e-18), label
    assert np.allclose(f['dqdt'], x['dqdt'] + py['dqdt'], rtol=RTOL, atol=1e-20), label
    assert np.allclose(f['rain'], py['rain'], rtol=RTOL, atol=1e-18), label
    assert f['shallow'] == (py['decision'] == 'thin'), label
    if py['decision'] == 'thin':
        assert f['top'] == py['top'] + 1, label                      # Fortran's levels count from 1
    assert f['layers'] == int(np.sum(py['rain'] > 0.0)), label        # cloud layers only where it rains
    for key, got in (('dtdt', f['dtdt'] - x['dtdt']), ('dqdt', f['dqdt'] - x['dqdt']), ('rain', f['rain'])):
        scale = float(np.max(np.abs(py[key])))
        if scale > 0.0:
            worst[key] = max(worst.get(key, 0.0), float(np.max(np.abs(got - py[key]))) / scale)
    if py['decision'] in ('deep', 'moisten'):                        # the moist enthalpy and water, in the Fortran
        m = x['dp'] * DSIGMA / C['ga']
        cp = C['acpd'] * (1.0 + C['ADV'] * qe)
        rain = f['rain'].sum()
        vapour = ((f['dqdt'] - x['dqdt']) * m).sum()
        heat = (cp * (f['dtdt'] - x['dtdt']) * m).sum()
        if py['decision'] == 'deep':
            assert abs(-vapour - rain) <= 1e-9 * rain and abs(heat - C['ALV'] * rain) <= 1e-8 * C['ALV'] * rain, label
        else:
            held = (qe * m).sum() / srd.TAU_S
            assert rain == 0.0 and abs(vapour) <= 1e-10 * held, label
            assert abs(heat) <= 1e-8 * max(1.0, py['pt']), label


@FORTRAN
def test_the_fortran_redesign_gives_the_python_copys_decisions_tendencies_and_rain(executables):
    rng = np.random.default_rng(13)
    cases = [tsc.column_inputs(t, q, ps, rng) for t, q, ps in columns()]
    exe = executables('both', NHOR)
    printed, _ = tsc.run_harness(exe, harness_text('redesign', [(100, cases)]))
    seen, worst = set(), {}
    for j, (x, f) in enumerate(zip(cases, parsed(printed, 1, len(cases))[0])):
        py, te, qe = python_step(x)
        assert_step_matches(f, py, x, te, qe, j, worst)
        seen.add(py['decision'])
    assert seen == set(srd.DECISIONS), seen                           # the cases do what they were chosen for
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


def sequence_steps(ncols=NHOR, nsteps=24, seed=17):
    """A sequence of steps on the soundings in which each column moves through the decisions: its sounding, vapour
    factor and low-level warming change every few steps, with small changes of the state from step to step."""
    rng = np.random.default_rng(seed)
    names = list(tsc.SOUNDINGS)
    steps = []
    for s in range(nsteps):
        cols = []
        for j in range(ncols):
            moist, warm = VARIANTS[(j + s // 3) % len(VARIANTS)]
            t, q, ps = tsc.sounding(names[(j + s // 9) % len(names)], moist, warm)
            wobble = 1.0 + 0.02 * np.sin(0.3 * s + j + np.arange(NLEV))
            cols.append(tsc.column_inputs(t + 0.2 * np.sin(0.5 * s + j), q * wobble, ps, rng))
        steps.append((7000 + s, cols))
    return steps


@FORTRAN
def test_the_fortran_redesign_follows_the_copy_step_for_step_through_every_decision(executables):
    steps = sequence_steps()
    ncols = len(steps[0][1])
    exe = executables('both', ncols)
    printed, _ = tsc.run_harness(exe, harness_text('redesign', steps))
    worst, decisions = {}, np.zeros((len(steps), ncols), int)
    for s, ((nstep, cols), frow) in enumerate(zip(steps, parsed(printed, len(steps), ncols))):
        for j, (x, f) in enumerate(zip(cols, frow)):
            py, te, qe = python_step(x)
            assert_step_matches(f, py, x, te, qe, (nstep, j), worst)
            decisions[s, j] = srd.DECISION_CODES[py['decision']]
    assert set(decisions.ravel()) == set(srd.DECISION_CODES.values())
    changes = {(int(a), int(b)) for col in decisions.T for a, b in zip(col[:-1], col[1:]) if a != b}
    assert any(a != 0 and b == 0 for a, b in changes)               # a decision left behind is cleared
    assert {(1, 4), (4, 1), (3, 4), (4, 3)} & changes
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


@FORTRAN
def test_with_the_switch_off_the_patched_routine_prints_every_digit_the_routine_without_it_prints(executables):
    rng = np.random.default_rng(19)
    cases = [tsc.column_inputs(t, q, ps, rng) for t, q, ps in columns()]
    route = tsc.route_l_steps(ncols=NHOR, nsteps=60)
    for setting, steps in (('design', [(7, cases), (8, cases[::-1])]), ('plain', [(7, cases), (8, cases[::-1])]),
                           ('closure', [(7, cases), (8, cases[::-1])]), ('route_l', route)):
        text = harness_text(setting, steps)
        printed = {which: tsc.run_harness(executables(which, NHOR), text)[0]
                   for which in (('base', 'closure', 'both') if setting in ('design', 'plain') else ('closure', 'both'))}
        assert printed['both'] == printed['closure'], setting       # every printed digit, so every bit
        if 'base' in printed:
            assert printed['both'] == printed['base'], setting
        rows = printed['both'].strip('\n').split('\n')
        assert any(float(v) > 0.0 for r in rows[2::5] for v in r.split()), setting   # the deep branch ran


def test_the_patch_applies_once_to_the_shared_source_and_to_a_redesign_runs_chain():
    marker = srd.REDESIGN_MARKER
    assert not any(marker in old for old, _ in srd.SBM_REDESIGN_EDITS)
    assert all(marker in new for _, new in srd.SBM_REDESIGN_EDITS)
    assert scl.CLOSURE_MARKER not in marker and marker not in scl.CLOSURE_MARKER
    for _, new in srd.SBM_REDESIGN_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines()), new[:60]
    patched = []
    if tsc.shared_rainmod() is not None:
        _, closure, both = sources()
        assert tsc.undo(both, srd.SBM_REDESIGN_EDITS, marker) == closure
        patched.append(both)
    if (er.UNTOUCHED_SOURCE / 'rainmod.f90').exists():
        patched.append(tsc.fresh_rainmod(dict(tsc.FRESH_SWITCHES, redesign=True)))
        for switches in (tsc.FRESH_SWITCHES, dict(tsc.FRESH_SWITCHES, closure=True)):
            assert marker not in tsc.fresh_rainmod(switches)        # a run without the redesign keeps its source
    if not patched:
        pytest.skip('no ExoPlaSim source at hand')
    for text in patched:
        assert tsc.apply(text, srd.SBM_REDESIGN_EDITS, marker) == text
        assert scl.CLOSURE_MARKER in text and text.count('      subroutine sbmred(') == 1
        for name in ('sbmclosc', 'sbmapply', 'sbmdilp', 'sbmromps', 'sbmrlhist'):
            assert text.count(f'      subroutine {name}(') + text.count(f'      subroutine {name}\n') == 1
        assert text.count('       if(nclos_bm == 2) then\n        call sbmred(jhor)\n        cycle\n') == 1
        assert text.count('elseif(lbuoy .and. nclos_bm == 1) then') == 1 and 'elseif(lbuoy) then' not in text
    chain = er.patch_chain(**tsc.DESIGN_SWITCHES, redesign=True)
    assert chain[-len(er.SBM_CLOSURE_PATCHES) - 1:] == er.SBM_CLOSURE_PATCHES + er.SBM_REDESIGN_PATCHES
    assert er.SBM_REDESIGN_PATCHES[0] not in er.patch_chain(**tsc.DESIGN_SWITCHES, closure=True)
    assert er.patch_chain(**tsc.DESIGN_SWITCHES, closure=True, redesign=True) == chain


@FORTRAN
def test_the_patched_rainmod_compiles_with_plasims_modules(tmp_path):
    path = tsc.shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(path.parent / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'rainmod.f90').write_text(sources()[2])
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------- runner

GOOD = dict(deep_hpa=200.0, shallow_cap=True, closure='redesign')


def test_the_runner_records_the_redesign_and_writes_its_namelist():
    r = er.convection_setting(**GOOD)
    assert r['tau_s'] == srd.TAU_S and 'rh' not in r and 'route_l' not in r
    record = r['closure']
    assert record['form'] == 'redesign' and record['switch'] == 'nclos_bm 2'
    assert record['constants'] == {name: value for name, value, _ in scl.CLOSURE_VALUES}
    assert record['departures'] == srd.DEPARTURES and record['drive_module_sha256'] == srd.DRIVE_SHA256
    nl = dict(er.convection_namelist(r))
    assert nl['NCLOS_BM'] == '2' and nl['TAU_BM'] == repr(srd.TAU_S) and 'RH_BM' not in nl and 'NROUTEL_BM' not in nl
    assert nl['DPDEEP_BM'] == '20000.0' and nl['NSHCAP_BM'] == '1' and nl['DIL_BM'] == repr(srd.DILUTION)
    assert {name.upper() for name, _, _ in scl.CLOSURE_VALUES} <= set(nl)
    assert not {name.upper() for name, _, _ in scl.ROUTE_L_VALUES} & set(nl)
    c = er.convection_setting(deep_hpa=200.0, shallow_cap=True, closure='c')
    assert dict(er.convection_namelist(c))['NCLOS_BM'] == '1' and c['closure']['form'] == 'c'
    for bad in (dict(closure='redesign'), dict(deep_hpa=200.0, closure='redesign'), dict(GOOD, deep_hpa=150.0),
                dict(GOOD, rh=0.7), dict(GOOD, tau_s=7200.0), dict(GOOD, eps=0.5), dict(GOOD, cin_j_kg=10.0),
                dict(GOOD, cape_j_kg=10.0), dict(GOOD, eps=0.5, mixing_steps=2), dict(GOOD, diagnostic_night=True),
                dict(GOOD, route_l=True)):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    assert er.convection_setting(**dict(GOOD, tau_s=srd.TAU_S)) == r
    assert er.uses_redesign_patch(GOOD) and not er.uses_closure_patch(GOOD)
    assert not er.uses_redesign_patch(dict(closure='c')) and not er.uses_redesign_patch(None)


def test_the_configuration_refuses_the_redesign_without_the_gust_or_with_route_l_before_any_patch(monkeypatch):
    def reached(*a, **k):
        raise AssertionError('ensure_patched reached: the refusal must come before the patch')
    monkeypatch.setattr(er, 'ensure_patched', reached)
    with pytest.raises(ValueError, match='gust'):
        er.configuration('A28R', 4, convection=GOOD, levels='l15')
    with pytest.raises(ValueError, match='route L'):
        er.configuration('A28R', 4, convection=dict(GOOD, route_l=True), levels='l15', surface_gust='beljaars')
    with pytest.raises(ValueError, match='rh'):
        er.configuration('A28R', 4, convection=dict(GOOD, rh=0.7), levels='l15', surface_gust='beljaars')


def test_ensure_patched_refuses_the_redesign_without_the_depth_test_the_cap_and_the_gust(tmp_path, monkeypatch):
    import sys
    import types
    package = tmp_path / 'exoplasim'
    (package / 'plasim' / 'src').mkdir(parents=True)
    monkeypatch.setitem(sys.modules, 'exoplasim', types.SimpleNamespace(__file__=str(package / '__init__.py')))
    for switches in (dict(), dict(depth=True), dict(depth=True, cap=True)):
        with pytest.raises(RuntimeError, match='redesigned closure'):
            er.ensure_patched(redesign=True, **switches)
    assert not list((package / 'plasim' / 'src').iterdir())                # nothing was written


def test_the_runners_flags_take_the_redesign_and_refuse_what_it_cannot_run(monkeypatch):
    seen = {}
    monkeypatch.setattr(er, 'run', lambda *a, **k: seen.update(args=a) or 0)
    base = ['A28R', '--convection', 'sbm', '--sbm-deep-hpa', '200', '--levels', 'l15', '--sbm-shallow-cap',
            '--surface-gust', 'beljaars_similarity', '--rain-drops', 'share']
    redesign = base + ['--sbm-closure', 'redesign']
    bad = (redesign + ['--sbm-route-l'], redesign + ['--sbm-rh', '0.7'], redesign + ['--sbm-tau', '7200'],
           redesign + ['--sbm-eps', '0.5'], redesign + ['--sbm-cape', '10'], redesign + ['--sbm-diagnostic-night'],
           redesign + ['--surface-gust', 'none'], [a for a in redesign if a != '--sbm-shallow-cap'],
           ['A28R', '--sbm-closure', 'redesign'], ['A28R', '--convection', 'massflux', '--sbm-closure', 'redesign'])
    # main sets the runner's signal handlers before it checks, so the test puts the process's own back
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        assert er.main(redesign) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['closure'] == 'redesign' and convection['route_l'] is None and convection['rh'] is None
        for argv in bad:
            with pytest.raises(SystemExit):
                er.main(argv)
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
