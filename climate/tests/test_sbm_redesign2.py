"""The second redesign of the simplified Betts-Miller scheme's closure (climate/gcm/sbm_redesign2.py,
SBM_REDESIGN2_PATCHES): the copy keeps the column's moist enthalpy and water for every decision, warms a layer only
toward a reference above it and never past it within deltsec2, rains only from the layers above their reference vapour
and moistens the layers below theirs with at most its share, takes the design scheme's parcel with its level of free
convection at or above the condensation level, and anchors its reference to the lowest level's air, the plume at no
entrainment being the finely integrated undiluted parcel. The routine with closure (c)'s, the first redesign's and this
patch, compiled with ExoPlaSim's flags against closure (c)'s stub of PlaSim's module and the gust patch's own gustcol,
gives the Python copy's decisions, tendencies, rain, shallow tops and convective cloud layers on soundings and over a
sequence of steps that moves each column through the decisions, to RTOL, without array temporaries under -fcheck=all,
and keeps the moist enthalpy and water there too; with nclos_bm 0, 1 or 2 it prints every digit the routine without this
patch prints, and with nclos_bm 0 every digit the routine without any of the three. The patch applies once to the shared
rainmod.f90 and to the untouched source with a second-redesign run's chain, after the first redesign's, keeps its lines
within 132 columns, holds its place in the patch catalogue, and the patched rainmod.f90 compiles with PlaSim's modules;
the runner records the form, writes its namelist and refuses what it cannot run with, before any patch."""
from __future__ import annotations

import re
import shutil
import signal
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import sbm
from climate.gcm import sbm_closure as scl
from climate.gcm import sbm_redesign as srd
from climate.gcm import sbm_redesign2 as sr2
from climate.tests import test_sbm_closure as tsc      # closure (c)'s soundings, stub, harness and sources
from climate.tests import test_sbm_redesign as tsr     # the first redesign's harness, its cloud layers printed

SIGMA, DSIGMA, HALF, NLEV, AKAP, C = tsc.SIGMA, tsc.DSIGMA, tsc.HALF, tsc.NLEV, tsc.AKAP, tsc.C
STEP2 = 2.0 * scl.STEP_S
# The Fortran against the copy: the same operations in the same order, with the same library's exp, log and pow and
# the sums accumulated level by level in both, apart from the rain's units (dsigma times the rain's share, then dp / ga
# in the harness, where the copy takes the layer's mass). A unit in the last place of the plume's temperature, near
# 250 K, becomes about 1e-14 of a tendency through the shortfall T_ref - T, two temperatures a kelvin or so apart, and
# the harness shows differences of that size (5e-14 at most); 1e-10 leaves room for that through the plume's 224 steps
# a column and catches any change of formula, which shows at 1e-6 or more.
RTOL = 1e-10
# (tau_bm, rh_bm, dpdeep_bm, nshcap_bm, nclos_bm, nroutel_bm): the earlier closures' settings with this one beside them
SETTINGS = dict(tsr.SETTINGS, redesign2=(sr2.TAU_S, 0.7, sr2.DEEP_PA, 1, 3, 0))
VARIANTS = tsr.VARIANTS + ((1.0, -6.0), (0.55, 0.0), (1.3, 0.0), (0.8, 2.0))   # (vapour factor, low-level warming, K)
COOL_NEXT_TO_SOURCE = 1.5          # K: the level above the source cooled, so that the parcel is buoyant there


def columns():
    """The soundings of closure (c)'s tests, each moistened or dried and its lowest levels warmed or cooled, and each
    with the level above the source cooled so that its parcel is buoyant next to the source."""
    out = [tsc.sounding(name, moist, warm) for name in tsc.SOUNDINGS for moist, warm in VARIANTS]
    for name in tsc.SOUNDINGS:
        t, q, ps = tsc.sounding(name)
        t = t.copy()
        t[-2] -= COOL_NEXT_TO_SOURCE
        out.append((t, q, ps))
    return out


NHOR = len(VARIANTS) * len(tsc.SOUNDINGS) + len(tsc.SOUNDINGS)       # every harness here runs this many columns


# ------------------------------------------------------------------------------------------------------- the copy

def test_the_copy_keeps_moist_enthalpy_and_water_and_heats_only_toward_its_reference():
    seen = dict.fromkeys(sr2.DECISIONS, 0)
    limited, at_share = set(), set()
    for t, q, ps in columns():
        o = sr2.column(t, q, SIGMA, DSIGMA, ps, C, AKAP)
        seen[o['decision']] += 1
        b = sr2.budgets(o, q, SIGMA, DSIGMA, ps, C)
        flux = sr2.enthalpy_flux(o, q, SIGMA, DSIGMA, ps, C)
        scale = max(1.0, C['ALV'] * b['rain'], abs(o['pt']))
        assert abs(b['enthalpy']) <= 1e-9 * scale and abs(flux[-1]) <= 1e-9 * scale
        if o['decision'] != 'deep':
            assert not np.any(o['dtdt']) and not np.any(o['dqdt']) and not np.any(o['rain'])
            continue
        relaxed = np.arange(NLEV) >= o['top']
        assert not np.any(o['dtdt'][~relaxed]) and not np.any(o['dqdt'][~relaxed])
        short = np.where(relaxed, np.maximum(o['t_ref'] - t, 0.0), 0.0)   # each layer's shortfall from its reference
        assert np.all((o['dtdt'] > 0.0) == (short > 0.0))                 # heat only into the layers below it
        reach = o['dtdt'][short > 0.0] * sr2.STEP2_S / short[short > 0.0]  # the share of it closed within deltsec2
        assert np.all(reach <= 1.0 + 1e-12) and np.ptp(reach) <= 1e-12   # never past it, in proportion to it
        assert (abs(reach[0] - 1.0) <= 1e-12) == o['limited']            # the rain held where the heat reaches it
        assert b['rain'] > 0.0 and abs(b['vapour'] + b['rain']) <= 1e-12 * b['rain']
        assert abs(b['heating'] - C['ALV'] * b['rain']) <= 1e-9 * C['ALV'] * b['rain']
        above = relaxed & (q > o['q_ref'])
        assert np.all((o['rain'] > 0.0) == above) and np.all(o['dqdt'][above] < 0.0)
        assert np.all(o['dqdt'][relaxed & ~above] >= 0.0)
        m = ps * DSIGMA / C['ga']
        moistening = float(np.sum(np.maximum(o['dqdt'], 0.0) * m)) * sr2.TAU_S
        assert abs(moistening - o['moist']) <= 1e-12 * o['excess']
        assert o['moist'] <= sr2.F_MOIST * o['excess'] * (1.0 + 1e-12) and o['moist'] <= o['deficit'] * (1.0 + 1e-12)
        assert b['rain'] * sr2.TAU_S <= (o['excess'] - o['moist']) * (1.0 + 1e-12)
        limited.add(o['limited'])
        at_share.add(o['moist'] < o['deficit'])
    assert all(seen.values()), seen                                  # every decision, on these soundings
    assert limited == {True, False} and at_share == {True, False}    # the heat's limit and the share, both ways


def test_the_copy_takes_the_variants_trigger_and_a_reference_anchored_to_the_lowest_level():
    moved = 0
    for t, q, ps in columns():
        o = sr2.column(t, q, SIGMA, DSIGMA, ps, C, AKAP)
        p = ps * SIGMA
        tp = sbm.parcel(t, q, SIGMA, p, 0.0, AKAP, C)
        assert np.array_equal(o['trigger_t'], tp)                        # the design scheme's parcel, bit for bit
        lcl = sbm.condensation_level(t[-1], q[-1], p[-1], p[0], AKAP, C)
        assert o['lcl'] == lcl
        buoyant = [k for k in range(NLEV - 2, -1, -1) if tp[k] > t[k] and p[k] <= lcl]
        old_lfc, old_top = srd.cloud_layer(t, tp)                        # the design scheme's first buoyant layer
        if old_lfc is not None and p[old_lfc] > lcl and buoyant and buoyant[0] != old_lfc:
            moved += 1
        if not buoyant:
            assert o['decision'] == 'none' and o['top'] is None
            continue
        top = o['top']
        assert o['lfc'] == buoyant[0] and all(tp[k] > t[k] for k in range(top, o['lfc'] + 1))
        assert top == 0 or tp[top - 1] <= t[top - 1]
        assert o['base'] == sbm.condensation_level(t[-1], q[-1], p[-1], p[top], AKAP, C)
        thin = o['base'] - p[top] <= sr2.DEEP_PA
        assert (o['decision'] == 'thin') == thin and (not thin or o['shallow_top'] == top)
        if thin:
            continue
        fine = sbm.diluted_parcel(t, q, SIGMA, ps, 0.0, AKAP, C, base=lcl, min_steps=sr2.NSUB)
        assert np.array_equal(sr2.zbp_reference(t[-1], q[-1], SIGMA, ps, C, AKAP, lcl, eps=0.0), fine)
        assert np.all(o['t_ref'][:-1] < fine[:-1] + 1e-9)                 # entraining, it is never warmer
        k = np.arange(NLEV)
        assert np.array_equal(o['q_ref'][k >= top], sr2.RH_REF * np.array(
            [sbm.saturation(t[j], p[j], C) for j in range(top, NLEV)]))
        assert not np.any(np.isfinite(o['q_ref'][k < top]))
        if o['decision'] == 'dry':
            assert o['shallow_top'] == next((j for j in range(top, NLEV) if p[j] >= o['base'] - sr2.DEEP_PA), NLEV - 1)
        for shift in (-2.0, 2.0):                                         # the air above the lowest level moved
            s = sr2.column(np.where(k < NLEV - 1, t + shift, t), q, SIGMA, DSIGMA, ps, C, AKAP)
            if s['t_ref'] is not None:
                assert np.array_equal(s['t_ref'], o['t_ref'])            # the reference does not follow it
    assert moved                                                         # the variant departs from the old trigger


# ------------------------------------------------------------------------------------------------------- Fortran

def sources():
    """mksbm's sources: the shared rainmod.f90 without closure (c)'s, the first redesign's and this patch ('base'),
    with closure (c)'s, with the first redesign's over it, and with this one over both."""
    path = tsc.shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    base = tsc.undo(tsc.without_redesign(path.read_text()), scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER)
    closure = tsc.apply(base, scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER)
    redesign = tsc.apply(closure, srd.SBM_REDESIGN_EDITS, srd.REDESIGN_MARKER)
    return dict(base=base, closure=closure, redesign=redesign,
                redesign2=tsc.apply(redesign, sr2.SBM_REDESIGN2_EDITS, sr2.REDESIGN2_MARKER))


def build(folder, code, nhor):
    folder.mkdir(parents=True, exist_ok=True)
    files = {'stub.f90': tsc.STUB.format(nhor=nhor, nlev=NLEV, declarations=scl._declarations() + sr2.declarations()),
             'gustcol.f90': tsc.gustcol_source(), 'mksbm.f90': code, 'harness.f90': tsr.HARNESS}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def harness_text(setting, steps):
    """The harness's input as closure (c)'s tests write it, with this form's setting beside theirs; the form's own
    constants are its namelist variables' defaults in the stub, the copy's values."""
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


@pytest.fixture(scope='module')
def executables(tmp_path_factory):
    """The harness built once on each source, as it is first needed, for nhor columns."""
    built, texts = {}, {}

    def get(which, nhor):
        if (which, nhor) not in built:
            if not texts:
                texts.update(sources())
            built[which, nhor] = build(tmp_path_factory.mktemp(f'{which}{nhor}'), tsc.routine(texts[which], 'mksbm'),
                                       nhor)
        return built[which, nhor]
    return get


def python_step(x):
    """The copy on one column's inputs as the routine takes them: the state entering the convection."""
    te, qe = x['dt'] + x['dtdt'] * STEP2, x['dq'] + x['dqdt'] * STEP2
    return sr2.column(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, STEP2), te, qe


def assert_step_matches(f, py, x, qe, label, worst):
    """The routine's outputs against the copy's, with the moist enthalpy and water of its tendencies; worst gathers
    the largest relative difference of the convection's own tendencies and rain."""
    assert f['route'] == sr2.DECISION_CODES[py['decision']], label
    assert np.allclose(f['dtdt'], x['dtdt'] + py['dtdt'], rtol=RTOL, atol=1e-18), label
    assert np.allclose(f['dqdt'], x['dqdt'] + py['dqdt'], rtol=RTOL, atol=1e-20), label
    assert np.allclose(f['rain'], py['rain'], rtol=RTOL, atol=1e-18), label
    assert f['shallow'] == (py['decision'] in ('thin', 'dry')), label
    if py['decision'] in ('thin', 'dry'):
        assert f['top'] == py['shallow_top'] + 1, label              # Fortran's levels count from 1
    assert f['layers'] == int(np.sum(py['rain'] > 0.0)), label        # cloud layers only where it rains
    for key, got in (('dtdt', f['dtdt'] - x['dtdt']), ('dqdt', f['dqdt'] - x['dqdt']), ('rain', f['rain'])):
        scale = float(np.max(np.abs(py[key])))
        if scale > 0.0:
            worst[key] = max(worst.get(key, 0.0), float(np.max(np.abs(got - py[key]))) / scale)
    if py['decision'] == 'deep':                                     # the moist enthalpy and water, in the Fortran
        m = x['dp'] * DSIGMA / C['ga']
        cp = C['acpd'] * (1.0 + C['ADV'] * qe)
        rain = f['rain'].sum()
        vapour = ((f['dqdt'] - x['dqdt']) * m).sum()
        heat = (cp * (f['dtdt'] - x['dtdt']) * m).sum()
        assert abs(-vapour - rain) <= 1e-9 * rain and abs(heat - C['ALV'] * rain) <= 1e-8 * C['ALV'] * rain, label


@tsc.FORTRAN
def test_the_fortran_second_redesign_gives_the_python_copys_decisions_tendencies_and_rain(executables):
    rng = np.random.default_rng(29)
    cases = [tsc.column_inputs(t, q, ps, rng) for t, q, ps in columns()]
    exe = executables('redesign2', NHOR)
    printed, _ = tsc.run_harness(exe, harness_text('redesign2', [(100, cases)]))
    seen, worst = set(), {}
    for j, (x, f) in enumerate(zip(cases, tsr.parsed(printed, 1, len(cases))[0])):
        py, te, qe = python_step(x)
        assert_step_matches(f, py, x, qe, j, worst)
        seen.add(py['decision'])
    assert seen == set(sr2.DECISIONS), seen                          # the cases do what they were chosen for
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


def sequence_steps(ncols=NHOR, nsteps=24, seed=31):
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


@tsc.FORTRAN
def test_the_fortran_second_redesign_follows_the_copy_step_for_step_through_every_decision(executables):
    steps = sequence_steps()
    ncols = len(steps[0][1])
    exe = executables('redesign2', ncols)
    printed, _ = tsc.run_harness(exe, harness_text('redesign2', steps))
    worst, decisions = {}, np.zeros((len(steps), ncols), int)
    for s, ((nstep, cols), frow) in enumerate(zip(steps, tsr.parsed(printed, len(steps), ncols))):
        for j, (x, f) in enumerate(zip(cols, frow)):
            py, te, qe = python_step(x)
            assert_step_matches(f, py, x, qe, (nstep, j), worst)
            decisions[s, j] = sr2.DECISION_CODES[py['decision']]
    assert set(decisions.ravel()) == set(sr2.DECISION_CODES.values())
    changes = {(int(a), int(b)) for col in decisions.T for a, b in zip(col[:-1], col[1:]) if a != b}
    assert any(a != 0 and b == 0 for a, b in changes)               # a decision left behind is cleared
    assert {(1, 5), (5, 1)} <= changes and ({(1, 3), (3, 1)} & changes)
    print('largest relative difference of the convection\'s tendencies and rain:', worst)
    assert all(v <= RTOL for v in worst.values()), worst


@tsc.FORTRAN
def test_with_the_switch_off_the_patched_routine_prints_every_digit_the_routine_without_it_prints(executables):
    rng = np.random.default_rng(37)
    cases = [tsc.column_inputs(t, q, ps, rng) for t, q, ps in columns()]
    route = tsc.route_l_steps(ncols=NHOR, nsteps=60)
    twice = [(7, cases), (8, cases[::-1])]
    for setting, steps, others in (('design', twice, ('base', 'closure', 'redesign')),
                                   ('plain', twice, ('base', 'closure', 'redesign')),
                                   ('closure', twice, ('closure', 'redesign')), ('route_l', route, ('closure', 'redesign')),
                                   ('redesign', twice, ('redesign',))):
        text = harness_text(setting, steps)
        mine = tsc.run_harness(executables('redesign2', NHOR), text)[0]
        for which in others:
            assert mine == tsc.run_harness(executables(which, NHOR), text)[0], (setting, which)   # every digit
        rows = mine.strip('\n').split('\n')
        assert any(float(v) > 0.0 for r in rows[2::5] for v in r.split()), setting   # the deep branch ran


def test_the_patch_applies_once_to_the_shared_source_and_to_a_second_redesign_runs_chain():
    marker = sr2.REDESIGN2_MARKER
    assert not any(marker in old for old, _ in sr2.SBM_REDESIGN2_EDITS)
    assert all(marker in new for _, new in sr2.SBM_REDESIGN2_EDITS)
    for other in (scl.CLOSURE_MARKER, srd.REDESIGN_MARKER):
        assert other not in marker and marker not in other
    for _, new in sr2.SBM_REDESIGN2_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines()), new[:60]
    patched = []
    if tsc.shared_rainmod() is not None:
        s = sources()
        assert tsc.undo(s['redesign2'], sr2.SBM_REDESIGN2_EDITS, marker) == s['redesign']
        patched.append(s['redesign2'])
    if (er.UNTOUCHED_SOURCE / 'rainmod.f90').exists():
        patched.append(tsc.fresh_rainmod(dict(tsc.FRESH_SWITCHES, redesign2=True)))
        for switches in (tsc.FRESH_SWITCHES, dict(tsc.FRESH_SWITCHES, closure=True),
                         dict(tsc.FRESH_SWITCHES, redesign=True)):
            assert marker not in tsc.fresh_rainmod(switches)        # a run without the form keeps its source
    if not patched:
        pytest.skip('no ExoPlaSim source at hand')
    for text in patched:
        assert tsc.apply(text, sr2.SBM_REDESIGN2_EDITS, marker) == text
        assert scl.CLOSURE_MARKER in text and srd.REDESIGN_MARKER in text
        for name in ('sbmred2', 'sbmzbp', 'sbmred', 'sbmapply', 'sbmclosc'):
            assert text.count(f'      subroutine {name}(') == 1, name
        assert text.count('       if(nclos_bm == 3) then\n        call sbmred2(jhor)\n        cycle\n') == 1
        assert text.index('call sbmred(jhor)') < text.index('call sbmred2(jhor)')
        for name, _, _ in sr2.VALUES:
            assert text.count(f"call mpbc{'i' if name in sr2.INTEGER_NAMES else 'r'}({name})") == 1, name
            assert len(re.findall(rf'namelist/rainmod_nl/.*\b{name}\b', text)) == 1, name
    chain = er.patch_chain(**tsc.DESIGN_SWITCHES, redesign2=True)
    assert chain[-len(er.SBM_CLOSURE_PATCHES) - 2:] == (er.SBM_CLOSURE_PATCHES + er.SBM_REDESIGN_PATCHES
                                                        + er.SBM_REDESIGN2_PATCHES)
    assert er.SBM_REDESIGN2_PATCHES[0] not in er.patch_chain(**tsc.DESIGN_SWITCHES, redesign=True)
    assert er.patch_chain(**tsc.DESIGN_SWITCHES, closure=True, redesign=True, redesign2=True) == chain


def test_the_catalogue_holds_the_patch_after_the_first_redesigns_in_every_runs_order():
    catalogue = [(n, m) for n, m, _ in er.patch_catalogue()]
    assert len(set(catalogue)) == len(catalogue)
    places = [catalogue.index(('rainmod.f90', m)) for m in (scl.CLOSURE_MARKER, srd.REDESIGN_MARKER,
                                                             sr2.REDESIGN2_MARKER)]
    assert places == sorted(places)
    for switches in (dict(tsc.DESIGN_SWITCHES, redesign2=True), dict(tsc.FRESH_SWITCHES, redesign2=True),
                     dict(tsc.DESIGN_SWITCHES, closure=True, redesign2=True)):
        chain = [(n, m) for n, m, _ in er.patch_chain(**switches)]
        assert set(chain) <= set(catalogue)
        for name in {n for n, _ in chain}:
            order = [catalogue.index(e) for e in dict.fromkeys(e for e in chain if e[0] == name)]
            assert order == sorted(order), name


@tsc.FORTRAN
def test_the_patched_rainmod_compiles_with_plasims_modules(tmp_path):
    path = tsc.shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(path.parent / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'rainmod.f90').write_text(sources()['redesign2'])
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *tsc.GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------- runner

GOOD = dict(deep_hpa=200.0, shallow_cap=True, closure='redesign2')


def test_the_runner_records_the_second_redesign_and_writes_its_namelist():
    r = er.convection_setting(**GOOD)
    assert r['tau_s'] == sr2.TAU_S and 'rh' not in r and 'route_l' not in r
    record = r['closure']
    assert record['form'] == 'redesign2' and record['switch'] == 'nclos_bm 3'
    assert record['constants'] == {name: value for name, value, _ in sr2.VALUES}
    assert record['departures'] == sr2.DEPARTURES and record['drive_module_sha256'] == sr2.DRIVE_SHA256
    assert record['drive_trigger_sha256'] == sr2.DRIVE_TRIGGER_SHA256
    nl = dict(er.convection_namelist(r))
    assert nl['NCLOS_BM'] == '3' and nl['TAU_BM'] == repr(sr2.TAU_S) and 'RH_BM' not in nl and 'NROUTEL_BM' not in nl
    assert nl['DPDEEP_BM'] == '20000.0' and nl['NSHCAP_BM'] == '1'
    assert {k: v for k, v in nl.items() if k in sr2.namelist_values()} == sr2.namelist_values()
    assert nl['NZBP_BM'] == '16' and nl['EPSZ_BM'] == '0.5' and nl['RHREF_BM'] == '0.7' and nl['FMOIST_BM'] == '0.1'
    assert not {name.upper() for name, _, _ in scl.CLOSURE_VALUES + scl.ROUTE_L_VALUES} & set(nl)
    first = er.convection_setting(deep_hpa=200.0, shallow_cap=True, closure='redesign')
    assert dict(er.convection_namelist(first))['NCLOS_BM'] == '2' and first['tau_s'] == scl.TAU_S
    for bad in (dict(closure='redesign2'), dict(deep_hpa=200.0, closure='redesign2'), dict(GOOD, deep_hpa=150.0),
                dict(GOOD, rh=0.7), dict(GOOD, tau_s=7200.0), dict(GOOD, tau_s=scl.TAU_S), dict(GOOD, eps=0.5),
                dict(GOOD, cin_j_kg=10.0), dict(GOOD, cape_j_kg=10.0), dict(GOOD, eps=0.5, mixing_steps=2),
                dict(GOOD, diagnostic_night=True), dict(GOOD, route_l=True), dict(GOOD, closure='redesign3')):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    assert er.convection_setting(**dict(GOOD, tau_s=sr2.TAU_S)) == r
    assert er.uses_redesign2_patch(GOOD) and not er.uses_redesign_patch(GOOD) and not er.uses_closure_patch(GOOD)
    assert not er.uses_redesign2_patch(dict(closure='redesign')) and not er.uses_redesign2_patch(None)


def test_the_configuration_refuses_the_second_redesign_without_the_gust_or_with_route_l_before_any_patch(monkeypatch):
    def reached(*a, **k):
        raise AssertionError('ensure_patched reached: the refusal must come before the patch')
    monkeypatch.setattr(er, 'ensure_patched', reached)
    with pytest.raises(ValueError, match='gust'):
        er.configuration('A28R', 4, convection=GOOD, levels='l15')
    with pytest.raises(ValueError, match='route L'):
        er.configuration('A28R', 4, convection=dict(GOOD, route_l=True), levels='l15', surface_gust='beljaars')
    with pytest.raises(ValueError, match='rh'):
        er.configuration('A28R', 4, convection=dict(GOOD, rh=0.7), levels='l15', surface_gust='beljaars')
    with pytest.raises(ValueError, match='tau_s'):
        er.configuration('A28R', 4, convection=dict(GOOD, tau_s=7200.0), levels='l15', surface_gust='beljaars')


def test_ensure_patched_refuses_the_second_redesign_without_the_depth_test_the_cap_and_the_gust(tmp_path, monkeypatch):
    import sys
    import types
    package = tmp_path / 'exoplasim'
    (package / 'plasim' / 'src').mkdir(parents=True)
    monkeypatch.setitem(sys.modules, 'exoplasim', types.SimpleNamespace(__file__=str(package / '__init__.py')))
    for switches in (dict(), dict(depth=True), dict(depth=True, cap=True)):
        with pytest.raises(RuntimeError, match='second redesign'):
            er.ensure_patched(redesign2=True, **switches)
    assert not list((package / 'plasim' / 'src').iterdir())                # nothing was written


def test_the_runners_flags_take_the_second_redesign_and_refuse_what_it_cannot_run(monkeypatch):
    seen = {}
    monkeypatch.setattr(er, 'run', lambda *a, **k: seen.update(args=a) or 0)
    base = ['A28R', '--convection', 'sbm', '--sbm-deep-hpa', '200', '--levels', 'l15', '--sbm-shallow-cap',
            '--surface-gust', 'beljaars_similarity', '--rain-drops', 'share']
    redesign2 = base + ['--sbm-closure', 'redesign2']
    bad = (redesign2 + ['--sbm-route-l'], redesign2 + ['--sbm-rh', '0.7'], redesign2 + ['--sbm-tau', '7200'],
           redesign2 + ['--sbm-tau', '64800'], redesign2 + ['--sbm-eps', '0.5'], redesign2 + ['--sbm-cape', '10'],
           redesign2 + ['--sbm-diagnostic-night'], redesign2 + ['--surface-gust', 'none'],
           [a for a in redesign2 if a != '--sbm-shallow-cap'], ['A28R', '--sbm-closure', 'redesign2'],
           ['A28R', '--convection', 'massflux', '--sbm-closure', 'redesign2'])
    # main sets the runner's signal handlers before it checks, so the test puts the process's own back
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        assert er.main(redesign2) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['closure'] == 'redesign2' and convection['route_l'] is None and convection['rh'] is None
        assert er.convection_setting(**convection)['closure']['form'] == 'redesign2'
        for argv in bad:
            with pytest.raises(SystemExit):
                er.main(argv)
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
