"""Closure (c) and route L in the simplified Betts-Miller scheme (climate/gcm/sbm_closure.py, SBM_CLOSURE_PATCHES): the
closure keeps the column's moist enthalpy and water, relaxes every level from the parcel's source to its top toward
Romps's humidity and hands a thin or too dry cloud to the shallow scheme with the cap's top; the copy's gradient,
interpolation and history follow numpy and the single column; route L starts only where a candidate has been moistened
and runs on for 12 hours after its last start. The routine with the patch, compiled with ExoPlaSim's flags against a
stub of PlaSim's module and with the gust patch's own gustcol, gives the Python copy's tendencies, rain, routes, source
levels, shallow tops, route L's history and its starts, column by column and step by step over a sequence that starts,
runs on and stops the route, to 1e-10, without array temporaries under -fcheck=all; with both switches off it prints
every digit the unpatched routine prints. The patch applies once to the shared rainmod.f90 and to the untouched source
with a closure run's chain, keeps its lines within 132 columns, and the patched rainmod.f90 compiles with PlaSim's
modules; the runner records closure (c) and route L, writes their namelist and refuses what they cannot run with."""
from __future__ import annotations

import shutil
import signal
import subprocess
import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import levels as lv
from climate.gcm import pbl_fortran as pbf
from climate.gcm import sbm
from climate.gcm import sbm_closure as scl
from climate.gcm import surface_gust as sg

GCM_FLAGS = ('-O3', '-cpp', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffpe-summary=none',
             '-finit-real=zero', '-fdefault-real-8')
FORTRAN = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')
HALF = np.asarray(lv.half_levels('l15'), float)
SIGMA = np.asarray(lv.full_levels(HALF), float)
DSIGMA = np.diff(np.r_[0.0, HALF])
NLEV = SIGMA.size
AKAP = 0.286
# PlaSim's constants as the design run takes them (its planet namelist; sbm_offline.run_constants' rule)
GASCON, GA = 288.376720919906, 1.6242188765608128
C = dict(sbm.PLASIM, gascon=GASCON, ga=GA, acpd=GASCON / AKAP, rdbrv=GASCON / 461.51,
         ADV=1870.0 / (GASCON / AKAP) - 1.0, akap=AKAP)
num = lambda xs: ' '.join(repr(float(x)) for x in np.atleast_1d(xs))

# Soundings from the single column's own runs under CM1's forcing (trigger_column.py with CLOSURE 'romps' and L_START
# 'closure', 3 November 2026; session-tools/agents/gap_boxes/fortran_cl/soundings.json on the research drive): the
# state entering the convection, top first, and the surface pressure, on PlaSim's fifteen levels (l15).
SOUNDINGS = {
    'L': dict(  # box_near_side_sea_l15_day4, the column's step 444
        ps=121292.890625,
        t=[167.4649824383741, 191.7245115407954, 214.6872966323238, 232.2926813085267, 250.73139795761884,
           265.7599729025477, 274.28408163781035, 280.78504261163374, 285.1063795564669, 288.30566649715047,
           291.5578101905033, 294.76989983959703, 297.4797706273561, 298.7966733732267, 299.00337362295505],
        q=[7.738917186224517e-08, 2.8010418607473018e-06, 1.001049069336602e-09, 0.0001473305391259981,
           0.0007039514299307748, 0.0023547078929638914, 0.004215596481079151, 0.005885842652915804,
           0.007839931939022424, 0.007719388172662102, 0.007719055410032019, 0.007710653657248875,
           0.007693117551157466, 0.007927326800446218, 0.01265600621488694],
    ),
    'L_unbuoyant_low': dict(  # box_near_side_sea_l15_day4, the column's step 592
        ps=121292.890625,
        t=[167.5409049203864, 193.5112193954941, 214.37110790347256, 232.14782978818718, 250.0516813906892,
           264.81272946418756, 275.70168539339005, 281.71988524627346, 284.9447315246402, 286.70692230497696,
           289.24531392841044, 292.32545573965655, 295.0108322140844, 297.07665529802557, 298.51771523459],
        q=[7.886348313415272e-08, 4.269813521291837e-06, 2.9344384853987755e-05, 0.00021335956698296388,
           0.0008215367353304496, 0.0025626773086275746, 0.004603165801803678, 0.005974636502017326,
           0.008712917568356503, 0.008569192294899877, 0.009504531524362782, 0.009358735879600589,
           0.00934884099036115, 0.009332158569836544, 0.009395540809410368],
    ),
    'low': dict(  # box_0e_design, the column's step 333
        ps=119192.28125,
        t=[169.52989147920294, 191.68567166151672, 216.6150277873958, 239.17243118158416, 257.71976414090614,
           269.57772431784423, 279.22225924668254, 285.414865626781, 288.8885446339109, 290.57954517636597,
           292.26257017105445, 293.937248900987, 295.964816626989, 298.0465826558138, 297.1776324701156],
        q=[1.0000000000000003e-09, 1.0000000000025467e-09, 1.0052813554942353e-09, 0.0005245214709830003,
           0.0020838590875818896, 0.004034923355007238, 0.0061314309289395305, 0.007797794922035407,
           0.010109902719470244, 0.010146214469419098, 0.011400285127540044, 0.01146199015827175,
           0.012063063005363708, 0.012048731211412776, 0.013031057970424679],
    ),
    'shallow': dict(  # box_near_side_sea_l15_day4, the column's step 37
        ps=121292.890625,
        t=[167.14384913175948, 192.09166544512567, 215.86787825340414, 233.33786140545263, 250.91839216519213,
           265.1094889920172, 274.0066651194345, 280.77067043890474, 285.2203907276592, 288.26517796985223,
           289.75671791568055, 291.34297288208796, 293.75754670833345, 295.8100946965954, 297.26208749881283],
        q=[7.082830387694053e-08, 1.0000000009080155e-09, 1.0003785128149426e-09, 9.707662055296598e-05,
           0.0007914648865405444, 0.0024512635151851428, 0.003950182680305508, 0.005057216660042478,
           0.005990025865352665, 0.006506764756074004, 0.008999432137843846, 0.011673272350416838,
           0.011684705643620646, 0.011763517771833867, 0.011828021087327154],
    ),
}


def sounding(name, moist=1.0, warm_ground=0.0):
    s = SOUNDINGS[name]
    t = np.array(s['t'], float) + warm_ground * (SIGMA > 0.95)
    return t, np.array(s['q'], float) * moist, float(s['ps'])


# ------------------------------------------------------------------------------------------------------- the copy

def deep_columns():
    out = []
    for name in SOUNDINGS:
        for moist, warm in ((1.0, 0.0), (1.05, 1.0), (1.1, 2.0)):
            t, q, ps = sounding(name, moist, warm)
            o = scl.closure(t, q, SIGMA, DSIGMA, ps, C, AKAP)
            if o['deep']:
                out.append((t, q, ps, o))
    return out


def test_closure_c_keeps_the_columns_moist_enthalpy_and_water_and_relaxes_from_its_source_to_its_top():
    cases = deep_columns()
    assert len(cases) >= 3
    for t, q, ps, o in cases:
        m = ps * DSIGMA / C['ga']
        cp = C['acpd'] * (1.0 + C['ADV'] * q)
        rain = o['rain'].sum()
        assert rain > 0.0 and abs(-(o['dqdt'] * m).sum() - rain) <= 1e-12 * rain
        assert abs((cp * o['dtdt'] * m).sum() - C['ALV'] * rain) <= 1e-9 * C['ALV'] * rain
        top = o['top']
        relaxed = np.arange(NLEV) >= top
        assert np.all(o['dtdt'][~relaxed] == 0.0) and np.all(o['dqdt'][~relaxed] == 0.0)
        assert np.all(o['dqdt'][relaxed] != 0.0)                    # the whole column from the source to the top
        ref = o['rh_ref'][relaxed]
        assert np.all((ref > 0.0) & (ref <= 1.0))
        r = scl.romps_reference(t, q, SIGMA, ps, o['reference_t'], top, C)
        assert o['rh_ref'][-1] == min(1.0, scl.DELTA_BODY / (scl.DELTA_BODY + r['gamma'][-1]))
        k = slice(top, NLEV - 1)
        assert np.all(r['delta'][k] >= scl.DELTA_BODY)
        expect = np.where(r['gamma'][k] > 0.0, np.minimum(1.0, r['delta'][k] / (r['delta'][k] + r['gamma'][k])), 1.0)
        assert np.array_equal(o['rh_ref'][k], expect)


def test_a_thin_cloud_keeps_the_parcels_top_and_one_too_dry_to_rain_takes_the_shallow_caps():
    found = {'thin': 0, 'dry': 0}
    for name in SOUNDINGS:
        for moist in (0.6, 0.8, 0.9, 1.0):
            t, q, ps = sounding(name, moist)
            o = scl.closure(t, q, SIGMA, DSIGMA, ps, C, AKAP)
            if o['top'] is None or o['deep']:
                continue
            assert np.all(o['dtdt'] == 0.0) and np.all(o['rain'] == 0.0)
            p = ps * SIGMA
            if o['base'] - p[o['top']] <= scl.DEEP_PA:
                found['thin'] += 1
                assert o['shallow_top'] == o['top']
            else:
                found['dry'] += 1
                k = o['shallow_top']
                assert k >= o['top'] and p[k] >= o['base'] - scl.DEEP_PA
                assert k == o['top'] or p[k - 1] < o['base'] - scl.DEEP_PA
    assert found['thin'] >= 1


def test_the_copys_gradient_and_interpolation_are_numpys():
    rng = np.random.default_rng(3)
    for _ in range(50):
        z = np.cumsum(rng.uniform(200.0, 9000.0, NLEV))[::-1]
        f = rng.normal(250.0, 30.0, NLEV)
        assert np.array_equal(scl.gradient(f, z), np.gradient(f, z))
        xp = np.cumsum(rng.uniform(10.0, 5000.0, 12))
        fp = rng.normal(11.0, 0.5, 12)
        for x in np.r_[rng.uniform(xp[0] - 100.0, xp[-1] + 100.0, 20), xp[3], xp[-1]]:
            assert scl.interp(float(x), xp, fp) == float(np.interp(x, xp, fp))


def test_the_history_entry_is_the_linear_change_of_relative_humidity():
    t, q, ps = sounding(next(iter(SOUNDINGS)))
    p = ps * SIGMA
    rad, dyn_t, dyn_q = np.full(NLEV, -1.2e-5), np.full(NLEV, 4.0e-6), q * 1.0e-6
    entry = scl.history_entry(t, q, p, rad, dyn_t, dyn_q, C)
    step = scl.STEP_S
    rh = lambda tt, qq: qq / np.array([sbm.saturation(a, b, C) for a, b in zip(tt, p)])
    exact = rh(t + step * (rad + dyn_t), q + step * dyn_q) - rh(t, q)
    # the column's linearisation, with its own latent heat and gas constant in place of the slope of PlaSim's
    # saturation: within 1% of the exact change through route L's layers (12% off at the top level's 167 K)
    lift = (SIGMA >= scl.LIFT_SIGMA[0]) & (SIGMA <= scl.LIFT_SIGMA[1])
    assert np.allclose(entry[lift], exact[lift], rtol=0.01, atol=0.0)


def test_route_l_starts_only_when_a_candidate_is_moistened_and_runs_on_twelve_hours():
    t, q, ps = sounding('L')
    state = scl.RouteL(NLEV)
    base = scl.convect(t, q, SIGMA, DSIGMA, ps, C, AKAP, route=dict(nstep=0, step=scl.STEP_S, drh24=np.zeros(NLEV),
                                                                     p_bl=ps, state=state))
    assert base['candidates'] and not base['start'] and not base['on'] and base['route'] != 'L'
    routes = []
    for n in range(1, 40):
        drh = np.full(NLEV, 0.01 if n == 5 else -0.01)
        out = scl.convect(t, q, SIGMA, DSIGMA, ps, C, AKAP, route=dict(nstep=n, step=scl.STEP_S, drh24=drh, p_bl=ps,
                                                                        state=state))
        routes.append(out['route'])
        if out['route'] == 'L':
            best = max(out['candidates'], key=lambda c: c[1])
            assert out['source'] == best[0]
    on = [n for n, r in zip(range(1, 40), routes) if r == 'L']
    assert on == list(range(5, 5 + 25))                             # the start and 24 steps of 30 minutes after it


# ------------------------------------------------------------------------------------------------------- Fortran

def shared_rainmod():
    path = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'rainmod.f90'
    return path if path.exists() else None


def apply(text, edits, marker):
    if marker in text:
        return text
    for old, new in edits:
        assert text.count(old) == 1, old[:60]
        text = text.replace(old, new)
    return text


def undo(text, edits, marker):
    if marker not in text:
        return text
    for old, new in reversed(edits):
        assert text.count(new) == 1, new[:60]
        text = text.replace(new, old)
    return text


def routine(text, name):
    a = text.index(f'      subroutine {name}\n')
    b = text.index(f'      end subroutine {name}\n') + len(f'      end subroutine {name}\n')
    return text[a:b]


STUB = """
      module rainmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEM={nlev}-1, NLEP={nlev}+1
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dqt(NHOR,NLEP), dp(NHOR)
      real :: du(NHOR,NLEP), dv(NHOR,NLEP), dudt(NHOR,NLEP), dvdt(NHOR,NLEP)
      real :: dshfl(NHOR), devap(NHOR), dust3(NHOR), dtdtswr(NHOR,NLEV), dtdtlwr(NHOR,NLEV), dttdyn(NHOR,NLEV)
      real :: sigma(NLEV), dsigma(NLEV), sigmah(NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV), dgp3d(NHOR,NLEV,20)
      integer :: icclev(NHOR,NLEV), icctot(NHOR)
      real :: akap, ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon, deltsec, deltsec2, pdeep
      real :: tau_bm, rh_bm, eps_bm=0., cin_bm=0., cape_bm=0., epst_bm=0., dpdeep_bm
      integer :: nsub_bm = 1, nshcap_bm = 0, ndiagnight_bm = 0
      integer :: nshallow = 1, ndiaggp = 0, nstep = 0, mypid = 0
      integer :: kshal(NHOR), ktopo(NHOR)
{declarations}      end module rainmod

      subroutine mkshallow(pt,pq,kshallow,klift,ktop,pdtdt,pdqdt)
      use rainmod
      integer kshallow(NHOR), klift(NHOR), ktop(NHOR)
      real pt(NHOR,NLEV), pq(NHOR,NLEV), pdtdt(NHOR,NLEV), pdqdt(NHOR,NLEV)
      kshal(:)=kshallow(:)
      ktopo(:)=ktop(:)
      pdtdt(:,:)=0.
      pdqdt(:,:)=0.
      end subroutine mkshallow
      subroutine sbmnight(pmuz)
      use rainmod
      real pmuz(NHOR)
      pmuz(:)=1.
      end subroutine sbmnight
      subroutine sbmdtt(pdtt)
      use rainmod
      real pdtt(NHOR,NLEV)
      pdtt(:,:)=dttdyn(:,:)
      end subroutine sbmdtt
      subroutine mpabort(ym)
      character(len=*) :: ym
      write(0,*) ym
      stop 1
      end subroutine mpabort
"""
HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: k, jh, js, jw, nsteps
      real :: zlast, zsum(NLEV), t1, t2, tsum
      read(*,*) sigma, dsigma, sigmah
      read(*,*) akap, ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon, deltsec
      read(*,*) tau_bm, rh_bm, dpdeep_bm, nshcap_bm, nclos_bm, nroutel_bm, nsteps
      deltsec2=deltsec+deltsec
      pdeep=999999.
      nshallow=1
      if(nclos_bm == 1 .and. nroutel_bm == 1) then
       allocate(drlhist(NHOR,NLEV,NRLWIN))
       allocate(drllast(NHOR))
       drlhist(:,:,:)=0.
       drllast(:)=-1.e30
      endif
      tsum=0.
      do js=1,nsteps
       read(*,*) nstep
       do jh=1,NHOR
        read(*,*) dp(jh), dshfl(jh), devap(jh), dust3(jh)
        read(*,*) (dt(jh,k), k=1,NLEV)
        read(*,*) (dq(jh,k), k=1,NLEV)
        read(*,*) (dtdt(jh,k), k=1,NLEV)
        read(*,*) (dqdt(jh,k), k=1,NLEV)
        read(*,*) (du(jh,k), k=1,NLEV)
        read(*,*) (dv(jh,k), k=1,NLEV)
        read(*,*) (dudt(jh,k), k=1,NLEV)
        read(*,*) (dvdt(jh,k), k=1,NLEV)
        read(*,*) (dtdtswr(jh,k), k=1,NLEV)
        read(*,*) (dtdtlwr(jh,k), k=1,NLEV)
        read(*,*) (dqt(jh,k), k=1,NLEV)
        read(*,*) (dttdyn(jh,k), k=1,NLEV)
       enddo
       kshal(:)=0
       ktopo(:)=0
       call cpu_time(t1)
       call mksbm
       call cpu_time(t2)
       tsum=tsum+(t2-t1)
       do jh=1,NHOR
        zlast=-1.e30
        zsum(:)=0.
        if(allocated(drllast)) then
         zlast=drllast(jh)
         do jw=1,NRLWIN
          zsum(:)=zsum(:)+drlhist(jh,:,jw)
         enddo
        endif
        write(*,'(*(es25.16e3))') (dtdt(jh,k), k=1,NLEV)
        write(*,'(*(es25.16e3))') (dqdt(jh,k), k=1,NLEV)
        write(*,'(*(es25.16e3))') (dprcl(jh,k)*dp(jh)/ga, k=1,NLEV)
        write(*,'(*(es25.16e3))') zsum
        write(*,'(4i6,es25.16e3)') kshal(jh), ktopo(jh), nsbmroute(jh), nsbmsrc(jh), zlast
       enddo
      enddo
      write(0,'(a,es14.6)') 'mksbm seconds ', tsum
      end program harness
"""


def gustcol_source() -> str:
    """The gust patch's gustcol with a stub of fluxmod holding the nonlocal scheme's constants at their defaults."""
    a = sg.GUST_CODE.index('      subroutine gustcol(')
    b = sg.GUST_CODE.index('      end subroutine gustcol\n') + len('      end subroutine gustcol\n')
    return ('      module fluxmod\n      real, parameter :: vonkarman=0.4\n' + pbf.declarations()
            + '      end module fluxmod\n\n' + sg.GUST_CODE[a:b])


def build(folder, code, nhor, flags=GCM_FLAGS):
    folder.mkdir(parents=True, exist_ok=True)
    files = {'stub.f90': STUB.format(nhor=nhor, nlev=NLEV, declarations=scl._declarations()),
             'gustcol.f90': gustcol_source(), 'mksbm.f90': code, 'harness.f90': HARNESS}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *flags, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def source_routines():
    """mksbm as the shared rainmod.f90 holds it without this patch, and with it."""
    path = shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = undo(path.read_text(), scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER)
    return routine(text, 'mksbm'), routine(apply(text, scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER), 'mksbm')


# (tau_bm, rh_bm, dpdeep_bm, nshcap_bm, nclos_bm, nroutel_bm): the design run's scheme, the scheme without the depth
# test and the cap, closure (c) alone and closure (c) with route L
SETTINGS = dict(design=(7200.0, 0.7, 20000.0, 1, 0, 0), plain=(7200.0, 0.7, 0.0, 0, 0, 0),
                closure=(scl.TAU_S, 0.7, scl.DEEP_PA, 1, 1, 0), route_l=(scl.TAU_S, 0.7, scl.DEEP_PA, 1, 1, 1))


def harness_text(setting, steps):
    """The harness's input: the levels, the constants, the setting, then each step's number and each column's
    surface pressure and fluxes, state entering the physics, tendencies before the convection, winds and their
    tendencies, radiation, and the dynamics' vapour and heating."""
    tau, rh, deep, cap, nclos, nroutel = SETTINGS[setting]
    lines = [num(np.r_[SIGMA, DSIGMA, HALF]),
             num([AKAP, C['ga'], C['acpd'], C['ADV'], C['ALV'], C['ALS'], C['TMELT'], C['rdbrv'], C['ra1'], C['ra2'],
                  C['ra4'], C['gascon'], scl.STEP_S]),
             f'{tau!r} {rh!r} {deep!r} {cap} {nclos} {nroutel} {len(steps)}']
    for nstep, cols in steps:
        lines.append(str(nstep))
        for x in cols:
            lines.append(num([x['dp'], x['dshfl'], x['devap'], x['dust3']]))
            lines += [num(x[k]) for k in ('dt', 'dq', 'dtdt', 'dqdt', 'du', 'dv', 'dudt', 'dvdt', 'swr', 'lwr', 'dqt',
                                          'dtt')]
    return '\n'.join(lines) + '\n'


def run_harness(exe, text):
    run = subprocess.run([str(exe)], cwd=exe.parent, input=text, capture_output=True, text=True, check=True)
    assert 'array temporary' not in run.stderr, 'the patched routine makes array temporaries under -fcheck=all'
    return run.stdout, run.stderr


def parsed(printed, nsteps, ncols):
    rows = printed.strip('\n').split('\n')
    assert len(rows) == 5 * nsteps * ncols
    out = []
    for s in range(nsteps):
        cols = []
        for j in range(ncols):
            r = rows[5 * (s * ncols + j):5 * (s * ncols + j + 1)]
            kshal, ktop, route, src = (int(v) for v in r[4].split()[:4])
            cols.append(dict(dtdt=np.array(r[0].split(), float), dqdt=np.array(r[1].split(), float),
                             rain=np.array(r[2].split(), float), drh24=np.array(r[3].split(), float),
                             shallow=bool(kshal), top=ktop, route=route, source=src, last=float(r[4].split()[4])))
        out.append(cols)
    return out


ROUTES = {None: 0, 'low': 1, 'L': 2, 'shallow': 3}


def column_inputs(t, q, ps, rng, unstable=True, moisten=0.0):
    """One column's inputs as PlaSim's physics finds them at the convection: the state entering the physics and the
    tendencies so far such that the state entering the convection is the sounding, a boundary layer heated (or
    cooled) from below, winds, the radiation's heating and the dynamics' vapour and heating; moisten scales the
    dynamics' moistening aloft (per second, as a share of the vapour)."""
    dtdt = np.full(NLEV, -1.0e-5) + 2.0e-5 * (SIGMA > 0.95)
    dqdt = np.where(SIGMA > 0.97, 4.0e-8, 0.0)
    swr = np.full(NLEV, 6.0e-6)
    lwr = dtdt - swr - 1.0e-6 * (SIGMA > 0.95)
    step2 = 2.0 * scl.STEP_S
    u = 4.0 + 6.0 * (1.0 - SIGMA) + rng.normal(0.0, 0.3, NLEV)
    return dict(dp=ps, dshfl=-60.0 if unstable else 15.0, devap=-4.0e-8 if unstable else -1.0e-9, dust3=0.2 ** 3,
                dt=t - step2 * dtdt, dq=q - step2 * dqdt, dtdt=dtdt, dqdt=dqdt, du=u, dv=0.5 * u[::-1],
                dudt=np.full(NLEV, -1.0e-5), dvdt=np.full(NLEV, 2.0e-5), swr=swr, lwr=lwr,
                dqt=q * moisten * ((SIGMA > 0.6) & (SIGMA < 0.95)), dtt=np.full(NLEV, -2.0e-6))


def python_step(x, nstep, setting, state):
    """The copy on one column's inputs as the routine takes them: the state entering the convection, route L's
    history and boundary-layer top."""
    step2 = 2.0 * scl.STEP_S
    te, qe = x['dt'] + x['dtdt'] * step2, x['dq'] + x['dqdt'] * step2
    tau = SETTINGS[setting][0]
    route = None
    drh24 = np.zeros(NLEV)
    if setting == 'route_l':
        p = x['dp'] * SIGMA
        entry = scl.history_entry(x['dt'], x['dq'], p, x['swr'] + x['lwr'], x['dtt'], x['dqt'], C)
        drh24 = state.record(nstep, entry)
        wth0, wq0 = scl.surface_fluxes(x['dshfl'], x['devap'], te[-1], qe[-1], p[-1], C)
        p_bl = scl.boundary_layer_pressure(te, qe, x['du'] + x['dudt'] * step2, x['dv'] + x['dvdt'] * step2, SIGMA,
                                           HALF, x['dp'], wth0, wq0, x['dust3'] ** (1.0 / 3.0), C)
        route = dict(nstep=nstep, step=scl.STEP_S, drh24=drh24, p_bl=p_bl, state=state)
    assert tau == scl.TAU_S
    out = scl.convect(te, qe, SIGMA, DSIGMA, x['dp'], C, AKAP, step2, route=route)
    out['drh24'] = drh24
    return out


def assert_step_matches(fortran, python, x, label):
    assert fortran['route'] == ROUTES[python['route']], label
    assert np.allclose(fortran['dtdt'], x['dtdt'] + python['dtdt'], rtol=1e-10, atol=1e-18), label
    assert np.allclose(fortran['dqdt'], x['dqdt'] + python['dqdt'], rtol=1e-10, atol=1e-20), label
    assert np.allclose(fortran['rain'], python['rain'], rtol=1e-10, atol=1e-18), label
    if python['route'] in ('low', 'L'):
        assert fortran['source'] == python['source'] + 1, label
    assert fortran['shallow'] == (python['route'] == 'shallow'), label
    if python['route'] == 'shallow':
        assert fortran['top'] == python['shallow_top'] + 1, label            # Fortran's levels count from 1
    assert np.allclose(fortran['drh24'], python['drh24'], rtol=1e-10, atol=1e-16), label


def closure_cases():
    rng = np.random.default_rng(11)
    cases = []
    for name in SOUNDINGS:
        for moist, warm in ((1.0, 0.0), (0.9, 0.0), (1.05, 1.0), (1.1, 2.0), (0.4, 0.0)):
            t, q, ps = sounding(name, moist, warm)
            cases.append(column_inputs(t, q, ps, rng))
    return cases


@FORTRAN
def test_the_fortran_closure_gives_the_python_copys_tendencies_rain_and_routes(tmp_path):
    cases = closure_cases()
    _, code = source_routines()
    exe = build(tmp_path / 'closure', code, nhor=len(cases))
    printed, _ = run_harness(exe, harness_text('closure', [(100, cases)]))
    fortran = parsed(printed, 1, len(cases))[0]
    seen = set()
    for j, (x, f) in enumerate(zip(cases, fortran)):
        py = python_step(x, 100, 'closure', None)
        assert_step_matches(f, py, x, j)
        seen.add(py['route'])
        if py['route'] == 'low':                                       # the moist enthalpy and water, in the Fortran
            m = x['dp'] * DSIGMA / C['ga']
            te, qe = x['dt'] + x['dtdt'] * 2.0 * scl.STEP_S, x['dq'] + x['dqdt'] * 2.0 * scl.STEP_S
            cp = C['acpd'] * (1.0 + C['ADV'] * qe)
            rain = f['rain'].sum()
            assert abs(-((f['dqdt'] - x['dqdt']) * m).sum() - rain) <= 1e-9 * rain
            assert abs((cp * (f['dtdt'] - x['dtdt']) * m).sum() - C['ALV'] * rain) <= 1e-8 * C['ALV'] * rain
    assert seen >= {'low', 'shallow', None}                            # the cases do what they were chosen for


def route_l_steps(ncols=4, nsteps=110, seed=5):
    """A sequence of steps on the soundings that starts route L, lets it run on after the moistening turns to drying
    and stops it, then starts it again: the columns' dynamics moisten their candidate layers for 12 steps, dry them
    from step 12, and moisten them again from step 80, with small changes of the state from step to step; one column
    in four takes a stable surface layer, where every layer in the range is a candidate."""
    rng = np.random.default_rng(seed)
    names = list(SOUNDINGS)
    steps = []
    for s in range(nsteps):
        cols = []
        for j in range(ncols):
            t, q, ps = sounding(names[j % len(names)])
            wobble = 1.0 + 0.02 * np.sin(0.3 * s + j + np.arange(NLEV))
            moisten = 2.0e-6 if (s < 12 or s >= 80) else -3.0e-6
            cols.append(column_inputs(t, q * wobble, ps, rng, unstable=(j % 4 != 3), moisten=moisten))
        steps.append((5000 + s, cols))
    return steps


@FORTRAN
def test_the_fortran_route_l_gives_the_python_copys_decisions_and_history_step_for_step(tmp_path):
    steps = route_l_steps()
    ncols = len(steps[0][1])
    _, code = source_routines()
    exe = build(tmp_path / 'route_l', code, nhor=ncols)
    printed, _ = run_harness(exe, harness_text('route_l', steps))
    fortran = parsed(printed, len(steps), ncols)
    states = [scl.RouteL(NLEV) for _ in range(ncols)]
    routes = []
    for (nstep, cols), frow in zip(steps, fortran):
        for j, (x, f) in enumerate(zip(cols, frow)):
            py = python_step(x, nstep, 'route_l', states[j])
            assert_step_matches(f, py, x, (nstep, j))
            if py['route'] in ('low', 'L'):                             # the moist enthalpy and water, in the Fortran
                k = py['source']
                te, qe = x['dt'] + x['dtdt'] * 2.0 * scl.STEP_S, x['dq'] + x['dqdt'] * 2.0 * scl.STEP_S
                qk = qe.copy()
                qk[k] = sbm.saturation(te[k], x['dp'] * SIGMA[k], C)      # route L's source, brought to saturation
                m = x['dp'] * DSIGMA / C['ga']
                cp = C['acpd'] * (1.0 + C['ADV'] * (qk if py['route'] == 'L' else qe))
                rain = f['rain'].sum()
                assert abs(-((f['dqdt'] - x['dqdt']) * m).sum() - rain) <= 1e-9 * rain, (nstep, j)
                assert abs((cp * (f['dtdt'] - x['dtdt']) * m).sum() - C['ALV'] * rain) <= 1e-8 * C['ALV'] * rain
            assert f['last'] == (float(states[j].last) if states[j].last is not None else -1.0e30), (nstep, j)
            routes.append((nstep, j, py['route'], py['start']))
    # the sequence starts the route, runs it on after its last start and stops it, then starts it again
    on = [(n, j) for n, j, r, _ in routes if r == 'L']
    starts = [(n, j) for n, j, _, s in routes if s]
    assert on and starts
    late = {j for n, j, r, s in routes if r == 'L' and not s}
    assert late, 'route L never ran on after a start'
    stopped = {j for n, j, r, _ in routes if r != 'L' and any(m < n and k == j for m, k in on)}
    assert stopped, 'route L never stopped'


@FORTRAN
def test_with_the_switches_off_the_patched_routine_prints_every_digit_the_unpatched_one_prints(tmp_path):
    cases = closure_cases()
    old, new = source_routines()
    exe_old = build(tmp_path / 'old', old, nhor=len(cases))
    exe_new = build(tmp_path / 'new', new, nhor=len(cases))
    rained = shallow = 0
    for setting in ('design', 'plain'):
        text = harness_text(setting, [(7, cases), (8, cases[::-1])])
        a, _ = run_harness(exe_old, text)
        b, _ = run_harness(exe_new, text)
        assert a == b                                               # every printed digit, so every bit
        rows = a.strip('\n').split('\n')
        rained += sum(any(float(v) > 0.0 for v in r.split()) for r in rows[2::5])
        shallow += sum(int(r.split()[0]) for r in rows[4::5])
    assert rained and shallow                                       # the deep and the shallow branches both ran


def fresh_rainmod(switches) -> str:
    """The untouched rainmod.f90 with the rainmod patches of a run's chain, in its order."""
    path = er.UNTOUCHED_SOURCE / 'rainmod.f90'
    if not path.exists():
        pytest.skip('no untouched ExoPlaSim source at hand')
    text = path.read_text()
    for name, marker, edits in er.patch_chain(**switches):
        if name == 'rainmod.f90':
            text = apply(text, edits, marker)
    return text


# The design run's switches; a fresh source takes the mass-flux patch too, whose text the convection's vapour
# diagnostic (pbl_fortran.DIAGNOSTIC_PATCHES, in a gust run's chain) edits, as the shared copy carries it.
DESIGN_SWITCHES = dict(twilight=True, depth=True, cap=True, gust=True, similarity=True, rain_drops=True)
FRESH_SWITCHES = dict(DESIGN_SWITCHES, massflux=True)


def test_the_patch_applies_once_to_the_shared_source_and_to_a_closure_runs_chain():
    marker = scl.CLOSURE_MARKER
    assert not any(marker in old for old, _ in scl.SBM_CLOSURE_EDITS)
    assert all(marker in new or 'Terluna' in new for _, new in scl.SBM_CLOSURE_EDITS)
    for _, new in scl.SBM_CLOSURE_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines()), new[:60]
    patched = []
    path = shared_rainmod()
    if path is not None:
        patched.append(apply(undo(path.read_text(), scl.SBM_CLOSURE_EDITS, marker), scl.SBM_CLOSURE_EDITS, marker))
    if (er.UNTOUCHED_SOURCE / 'rainmod.f90').exists():
        patched.append(fresh_rainmod(dict(FRESH_SWITCHES, closure=True)))
        assert marker not in fresh_rainmod(FRESH_SWITCHES)                 # a run without closure (c) keeps its source
    if not patched:
        pytest.skip('no ExoPlaSim source at hand')
    for text in patched:
        assert apply(text, scl.SBM_CLOSURE_EDITS, marker) == text
        for name in ('sbmclosc', 'sbmapply', 'sbmrlhist', 'sbmdtt'):
            assert text.count(f'      subroutine {name}\n') + text.count(f'      subroutine {name}(') == 1
        for name in scl.NAMELIST_NAMES:
            assert text.count(f'call mpbc{"i" if name in scl.INTEGER_NAMES else "r"}({name})') == 1
        assert text.count('if(nclos_bm == 1) then\n        call sbmclosc(jhor)') == 1
        assert text.count("call mpputgp('dsbmrlhist'") == 1 and text.count("call mpgetgp('dsbmrlhist'") == 1
        assert text.count("restart_has_array('dsbmrlhist'") == 1
    chain = er.patch_chain(**DESIGN_SWITCHES, closure=True)
    assert chain[-len(er.SBM_CLOSURE_PATCHES):] == er.SBM_CLOSURE_PATCHES             # last, over every other patch
    assert ('restartmod.f90', er.SWITCH_MARKER, er.SWITCH_RESTART_EDITS) in er.SBM_CLOSURE_PATCHES
    assert not any(p in er.patch_chain(**DESIGN_SWITCHES) for p in er.SBM_CLOSURE_PATCHES[:1])


@FORTRAN
def test_the_patched_rainmod_compiles_with_plasims_modules(tmp_path):
    path = shared_rainmod()
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = path.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(src / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    text = apply(undo(path.read_text(), scl.SBM_CLOSURE_EDITS, scl.CLOSURE_MARKER), scl.SBM_CLOSURE_EDITS,
                 scl.CLOSURE_MARKER)
    (tmp_path / 'rainmod.f90').write_text(text)
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------- runner

def test_the_runner_records_closure_c_and_route_l_and_writes_their_namelist():
    design = er.convection_setting(rh=0.7, deep_hpa=200.0, shallow_cap=True)
    assert 'closure' not in design and 'route_l' not in design
    assert er.convection_setting(rh=0.7, deep_hpa=200.0, shallow_cap=True, closure='design') == design
    c = er.convection_setting(deep_hpa=200.0, shallow_cap=True, closure='c')
    assert c['tau_s'] == scl.TAU_S and 'rh' not in c and c['closure']['form'] == 'c' and 'route_l' not in c
    assert c['closure']['constants'] == {name: value for name, value, _ in scl.CLOSURE_VALUES}
    cl = er.convection_setting(deep_hpa=200.0, shallow_cap=True, closure='c', route_l=True)
    assert cl['route_l']['constants'] == {name: value for name, value, _ in scl.ROUTE_L_VALUES}
    assert cl['route_l']['history_steps'] == scl.WINDOW and cl['route_l']['departures']
    nl = dict(er.convection_namelist(cl))
    assert nl['NCLOS_BM'] == '1' and nl['NROUTEL_BM'] == '1' and nl['TAU_BM'] == repr(scl.TAU_S)
    assert 'RH_BM' not in nl and nl['DIL_BM'] == repr(scl.DILUTION) and nl['XLV_BM'] == repr(scl.XLV)
    assert set(scl.namelist_values()) <= set(nl) and nl['DPDEEP_BM'] == '20000.0' and nl['NSHCAP_BM'] == '1'
    assert 'NROUTEL_BM' not in dict(er.convection_namelist(c))
    assert 'NCLOS_BM' not in dict(er.convection_namelist(design))
    good = dict(deep_hpa=200.0, shallow_cap=True, closure='c')
    for bad in (dict(closure='c'), dict(deep_hpa=200.0, closure='c'), dict(good, deep_hpa=150.0), dict(good, rh=0.7),
                dict(good, tau_s=7200.0), dict(good, eps=0.5), dict(good, cin_j_kg=10.0), dict(good, cape_j_kg=10.0),
                dict(good, diagnostic_night=True), dict(good, closure='design', route_l=True), dict(good, closure='b')):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    assert er.convection_setting(**dict(good, tau_s=scl.TAU_S)) == c
    assert er.uses_closure_patch(dict(closure='c')) and not er.uses_closure_patch(dict(closure='design'))
    assert not er.uses_closure_patch(None)


def test_the_configuration_refuses_closure_c_without_the_gust_or_on_another_step_before_any_patch(monkeypatch):
    def reached(*a, **k):
        raise AssertionError('ensure_patched reached: the refusal must come before the patch')
    monkeypatch.setattr(er, 'ensure_patched', reached)
    c = dict(deep_hpa=200.0, shallow_cap=True, closure='c', route_l=True)
    with pytest.raises(ValueError, match='gust'):
        er.configuration('A28R', 4, convection=c, levels='l15')
    with pytest.raises(ValueError, match='30-minute'):
        er.configuration('A28R', 4, overrides=dict(timestep_min=15.0), convection=c, levels='l15',
                         surface_gust='beljaars')
    with pytest.raises(ValueError, match='rh'):
        er.configuration('A28R', 4, convection=dict(c, rh=0.7), levels='l15', surface_gust='beljaars')


def test_ensure_patched_refuses_closure_c_without_the_depth_test_the_cap_and_the_gust(tmp_path, monkeypatch):
    import sys
    import types
    package = tmp_path / 'exoplasim'
    (package / 'plasim' / 'src').mkdir(parents=True)
    monkeypatch.setitem(sys.modules, 'exoplasim', types.SimpleNamespace(__file__=str(package / '__init__.py')))
    for switches in (dict(), dict(depth=True), dict(depth=True, cap=True)):
        with pytest.raises(RuntimeError, match='closure'):
            er.ensure_patched(closure=True, **switches)
    assert not list((package / 'plasim' / 'src').iterdir())                # nothing was written


def test_the_runners_flags_take_closure_c_and_route_l_and_refuse_what_they_cannot_run(monkeypatch):
    seen = {}
    monkeypatch.setattr(er, 'run', lambda *a, **k: seen.update(args=a) or 0)
    base = ['A28R', '--convection', 'sbm', '--sbm-deep-hpa', '200', '--levels', 'l15', '--sbm-shallow-cap',
            '--surface-gust', 'beljaars_similarity', '--rain-drops', 'share']
    bad = (base + ['--sbm-route-l'], base + ['--sbm-closure', 'c', '--sbm-rh', '0.7'],
           base + ['--sbm-closure', 'c', '--sbm-tau', '7200'],
           base + ['--sbm-closure', 'c', '--sbm-route-l', '--surface-gust', 'none'],
           ['A28R', '--sbm-closure', 'c'], ['A28R', '--convection', 'massflux', '--sbm-route-l'])
    # main sets the runner's signal handlers before it checks, so the test puts the process's own back
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        assert er.main(base + ['--sbm-closure', 'c', '--sbm-route-l']) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['closure'] == 'c' and convection['route_l'] is True
        assert er.main(base + ['--sbm-closure', 'c']) == 0
        convection = next(a for a in seen['args'] if isinstance(a, dict) and 'closure' in a)
        assert convection['closure'] == 'c' and convection['route_l'] is None
        for argv in bad:
            with pytest.raises(SystemExit):
                er.main(argv)
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
