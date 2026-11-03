"""The hollows' store in PlaSim's land (climate/gcm/land_store.py): the pools' rate is Penman's exact form of the
no-storage balance with PlaSim's surface transfer over water; the store's curve is exact for paraboloids and nil at nil
water; the runoff share is bilinear within its grids; each land step conserves water exactly, keeps the store's water
from going negative and spills only above capacity; PlaSim's soil step with the patch, compiled with ExoPlaSim's flags
against PlaSim's own module, gives the Python copy's water to 1e-10 with the switch on and prints every digit the
unpatched soil step prints with it off; every edit only inserts lines and the inserted statements run only with
nlstore 1; the patch applies once to the shared landmod.f90 and compiles with PlaSim's modules; and the runner records
the store, writes its tables and namelist, and refuses a run without the 10 m wind."""
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import land_store as ls
from climate.tests.test_pbl_fortran import FORTRAN, G_MOON, GCM_FLAGS, LUNAR, apply, num

PLASIM = dict(LUNAR, ga=G_MOON, akap=0.286, **ls.PLASIM_CONSTANTS)
PC = ls.pool_constants(PLASIM, G_MOON, 121590.0)
DT = 1800.0


def twilight_source(name):
    path = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / name
    return path if path.exists() else None


def tables(n, seed=3):
    """A store's tables for n cells: paraboloid bowls (area as the root of the fill) and a deeper depression whose
    area grows more slowly, a seal-like seepage, routing shares and a runoff table rising with depth and wetness."""
    rng = np.random.default_rng(seed)
    fill = np.array([0.0, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1.0])
    bowls = rng.uniform(0.05, 0.2, n)[:, None] * np.sqrt(fill)[None, :]
    deps = rng.uniform(0.0, 0.03, n)[:, None] * fill[None, :] ** 0.35
    area = bowls + deps
    seep = area * rng.uniform(2.5e-8, 4e-8, n)[:, None] * (1.0 + 0.3 * fill[None, :])
    depths = np.array([0.0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100, 150, 200, 300.0])
    wets = np.array([0.0, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.5, 0.75, 1.0])
    table = np.clip((depths[:, None] - 1.0) / (depths[:, None] + 15.0), 0.0, 1.0) * (0.6 + 0.4 * wets[None, :])
    shares = np.stack([rng.uniform(0.4, 0.8, n), rng.uniform(0.0, 0.01, n), rng.uniform(0.0, 0.02, n)], axis=1)
    cap = rng.uniform(0.5, 20.0, n)
    return dict(fill=fill, area=area, seep=seep, cap=cap, shares=shares, depths=depths, wets=wets, table=table,
                wmax=np.full(n, 0.5))


def forcing(n, seed=5):
    rng = np.random.default_rng(seed)
    return dict(prc=rng.choice([0.0, 2e-7, 3e-6, 2e-5], n), prl=rng.choice([0.0, 1e-7, 8e-7], n),
                devap_bucket=-rng.uniform(0.0, 4e-8, n), pool_rate=rng.uniform(-2e-9, 1.2e-7, n))


# ------------------------------------------------------------------------------------------------------- physics

def test_the_pools_rate_is_penmans_exact_form_with_plasims_transfer_over_water():
    sw = np.array([650.0, 300.0, 0.0, 0.0, 900.0])
    lw = np.array([440.0, 430.0, 455.0, 400.0, 420.0])
    ta = np.array([304.0, 303.0, 299.0, 296.0, 310.0])
    qa = np.array([0.009, 0.010, 0.012, 0.006, 0.004])
    p = np.full(5, 116000.0)
    wind = np.array([2.0, 1.2, 0.4, 0.05, 6.0])
    e, tw, v = ls.pool_evaporation(sw, lw, ta, qa, p, wind, PC)
    epsv = 1.0 / PC[2] - 1.0
    rho = p / (PC[1] * ta * (1.0 + epsv * qa))
    s = np.maximum(wind, PC[14])
    # the balance closes at the water's temperature with the transfer there
    assert np.abs(ls.pool_residual(tw, ls.pool_transfer(tw, ta, qa, p, s, PC), sw, lw, ta, qa, p, rho, PC)).max() < 1e-6
    # an independent root of the same balance by bisection
    for k in range(5):
        f = lambda t: float(ls.pool_residual(np.array(t), ls.pool_transfer(np.array(t), ta[k], qa[k], p[k], s[k], PC),
                                             sw[k], lw[k], ta[k], qa[k], p[k], rho[k], PC))
        lo, hi = ta[k] - 39.0, ta[k] + 39.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0.0 else (lo, mid)
        assert tw[k] == pytest.approx(0.5 * (lo + hi), abs=1e-6)
    assert e == pytest.approx(rho * v * (ls.saturation_humidity(tw, p, PC) - qa) / 1000.0, rel=1e-12)
    # sunlit pools evaporate and run warm; the night's pool runs cooler than its air
    assert (e[[0, 1, 4]] > 0.0).all() and (tw[[0, 4]] > ta[[0, 4]]).all() and tw[3] < ta[3]


def test_the_transfer_is_louis_form_at_the_winds_height_with_the_free_convection_floor():
    cn = (0.4 / math.log(PC[13] / PC[12] + 1.0)) ** 2
    ta, qa, p, s = 300.0, 0.01, 116000.0, 3.0
    epsv = 1.0 / PC[2] - 1.0
    qs = lambda t: float(ls.saturation_humidity(np.array(t), p, PC))
    # neutral: the water's virtual temperature equals the air's
    tw = 300.0
    for _ in range(50):
        tw = ta * (1.0 + epsv * qa) / (1.0 + epsv * qs(tw))
    assert float(ls.pool_transfer(np.array(tw), ta, qa, p, s, PC)) == pytest.approx(s * cn, rel=1e-9)
    # unstable and stable as surflx forms Louis's factors (b = c = d = 5)
    for t in (303.0, 297.0):
        thva, thvw = ta * (1.0 + epsv * qa), t * (1.0 + epsv * qs(t))
        ri = PC[0] * PC[13] * (thva - thvw) / (s * s * thva)
        fh = (1.0 - 15.0 * ri / (1.0 + 75.0 * cn * math.sqrt(-ri * (PC[13] / PC[12] + 1.0))) if ri <= 0.0
              else 1.0 / (1.0 + 15.0 * ri * math.sqrt(1.0 + 5.0 * ri)))
        assert float(ls.pool_transfer(np.array(t), ta, qa, p, s, PC)) == pytest.approx(s * cn * fh, rel=1e-12)
    # in calm air a warm pool's transfer is the tank law's at lunar gravity and 1.2 atm
    v = float(ls.pool_transfer(np.array(305.0), ta, qa, p, PC[14], PC))
    dthv = 305.0 * (1.0 + epsv * qs(305.0)) - ta * (1.0 + epsv * qa)
    assert v == pytest.approx(PC[15] * dthv ** (1.0 / 3.0), rel=1e-12)
    assert PC[15] == pytest.approx(0.0016 * (G_MOON / 9.80665 * 101325.0 / 121590.0) ** (1.0 / 3.0), rel=1e-12)


def test_the_curve_is_exact_for_paraboloids_and_nil_at_nil_water():
    t = tables(4)
    fill = t['fill']
    area = 0.2 * np.sqrt(fill)[None, :] * np.ones((4, 1))
    seep = 3e-8 * area
    cap = np.array([1.0, 2.0, 5.0, 0.0])
    for f in (0.0, 2e-7, 4e-5, 0.002, 0.07, 0.5, 1.0, 3.0):
        fw, q = ls.curve(fill, area, seep, cap, f * cap)
        expect = np.where(cap > 0.0, 0.2 * math.sqrt(min(f, 1.0)), 0.0)
        assert fw == pytest.approx(expect, abs=1e-15) and q == pytest.approx(3e-8 * expect, abs=1e-22)
    fw, q = ls.curve(fill, t['area'], t['seep'], t['cap'], np.zeros(4))
    assert (fw == 0.0).all() and (q == 0.0).all()
    # at the table's own fills the curve returns the table
    for k, f in enumerate(fill[1:], 1):
        fw, q = ls.curve(fill, t['area'], t['seep'], t['cap'], f * t['cap'])
        assert fw == pytest.approx(t['area'][:, k], rel=1e-12) and q == pytest.approx(t['seep'][:, k], rel=1e-12)


def test_the_runoff_share_is_bilinear_within_its_grids():
    t = tables(1)
    d, w, tab = t['depths'], t['wets'], t['table']
    for i in range(d.size):
        for j in range(w.size):
            assert float(ls.runoff_share(d, w, tab, d[i], w[j])) == pytest.approx(tab[i, j], abs=1e-15)
    mid = float(ls.runoff_share(d, w, tab, 0.5 * (d[6] + d[7]), 0.25 * w[3] + 0.75 * w[4]))
    expect = 0.5 * (0.25 * tab[6, 3] + 0.75 * tab[6, 4]) + 0.5 * (0.25 * tab[7, 3] + 0.75 * tab[7, 4])
    assert mid == pytest.approx(expect, rel=1e-12)
    assert float(ls.runoff_share(d, w, tab, 1e4, 2.0)) == tab[-1, -1]
    assert float(ls.runoff_share(d, w, tab, -1.0, -1.0)) == tab[0, 0] == 0.0


def test_each_land_step_conserves_water_and_keeps_the_stores_limits():
    n = 400
    t = tables(n)
    rng = np.random.default_rng(11)
    state = dict(w=rng.uniform(0.0, 0.03, n), v=np.where(rng.random(n) < 0.3, 0.0, rng.uniform(0.0, 0.2, n)), fw=None)
    state['v'][:20] = t['cap'][:20] * 0.999                  # some nearly full stores, which spill
    state['fw'], _ = ls.curve(t['fill'], t['area'], t['seep'], t['cap'], state['v'])
    spilled = np.zeros(n)
    for step in range(300):
        f = forcing(n, seed=step)
        if step < 3:
            f['prc'][:20] = 5e-4                            # a downpour onto the full stores
        v0 = state['v'].copy()
        state, flux, b = ls.land_step(DT, state, f, t, 0.3)
        lhs = b['rain']
        rhs = (b['bucket_change'] + b['bucket_evaporation'] + b['to_rivers'] + b['store_change']
               + b['store_evaporation'] + b['slow_store'])
        # exact to rounding on the largest of the budget's terms and the store's own water (up to 20 m here)
        scale = np.max([np.abs(b[k]) for k in b] + [v0, state['v']], axis=0)
        assert (np.abs(rhs - lhs) <= 1e-14 * np.maximum(scale, 1e-3)).all()
        assert (b['clamp'] == 0.0).all() and (state['v'] >= 0.0).all() and (state['v'] <= t['cap'] + 1e-15).all()
        assert (state['fw'][state['v'] == 0.0] == 0.0).all() and (flux['es'] * DT <= v0 + 1e-18).all()
        # the store spills only what lies above its capacity, and only when full
        assert (flux['spill'][state['v'] < t['cap'] * (1.0 - 1e-12)] == 0.0).all()
        spilled += flux['spill']
    assert (spilled[:20] > 0.0).all()


def test_the_stores_evaporation_is_at_most_its_water_and_none_under_snow():
    es = ls.store_evaporation(np.array([0.1, 0.1, 0.1, 0.0]), np.array([1e-6, 1.0, 1.0, 1.0]),
                              np.array([1e-6, 1e-7, -1e-8, 1e-7]), DT, snow=np.array([0.0, 0.0, 0.0, 0.0]))
    assert es[0] == pytest.approx(1e-6 / DT) and es[1] == pytest.approx(1e-8) and es[2] == pytest.approx(-1e-9)
    assert es[3] == 0.0
    assert ls.store_evaporation(np.array([0.1]), np.array([1.0]), np.array([1e-7]), DT, snow=np.array([0.01]))[0] == 0.0


# ------------------------------------------------------------------------------------------------------- Fortran

def patched(text):
    return apply(text, ls.LANDMOD_EDITS, ls.LAND_STORE_MARKER)


def soil_routines(text):
    start = text.index('      module landmod\n')
    end = text.index('      end module landmod\n') + len('      end module landmod\n')
    a = text.index('      subroutine tands\n')
    b = text.index('      end subroutine soilstep\n') + len('      end subroutine soilstep\n')
    return text[start:end] + text[a:b]


STUBS = """
      subroutine mpgagp(pf,pp,klev)
      real pf(*),pp(*)
      integer klev
      end subroutine mpgagp
      subroutine mpscgp(pf,pp,klev)
      real pf(*),pp(*)
      integer klev
      end subroutine mpscgp
      subroutine mpbcrn(p,n)
      real p(*)
      integer n
      end subroutine mpbcrn
"""
HARNESS = """
      program harness
      use landmod
      implicit none
      integer :: mode, n, nsteps, j, js, k
      real :: zsig, zevap(NHOR)
      read(*,*) mode, n, nsteps
      read(*,*) ga, gascon, akap, rdbrv, acpd, ra1, ra2, ra4, deltsec, zsig
      sigma(NLEV) = zsig
      dsigma(NLEV) = 2.*(1.-zsig)
      dls(:) = 0.
      do j=1,n
       read(*,*) dls(j), dprc(j), dprl(j), dprs(j), zevap(j), dwatc(j), dwmax(j), dsnowz(j), dglac(j)
       read(*,*) dts(j), dtsm(j), dshfl(j), dlhfl(j), dflux(j,NLEP), (dsoilt(j,k),k=1,NLSOIL)
      enddo
{store_input}
      do js=1,nsteps
       devap(:) = zevap(:)
{store_step}
       call soilstep
       do j=1,n
        write(*,'(*(es25.16e3))') dwatc(j), drunoff(j), dwater(j), dts(j), dsnowz(j), dsoilt(j,1), dsoilt(j,NLSOIL) &
     &   {store_out}
       enddo
      enddo
      end program harness
"""
STORE_INPUT = """
      if (mode == 1) then
       nlstore = 1
       ndiaggp2d = NLSOUT
       allocate(dgp2d(NHOR,NLSOUT))
       nafter = 0
       read(*,*) dlswet, dlssbk, dlsalbw, dlsemw, dlsz0w, dlszw, dlsumin, dlscfree, dlslb, dlslc, dlsld, alv, tmelt
       read(*,*) dlsfill
       read(*,*) dlsdgrid
       read(*,*) dlswgrid
       read(*,*) dlsrtab
       do j=1,n
        read(*,*) dlscap(j), dlsv(j), dlsfw(j), (dlsarea(j,k),k=1,NLSF), (dlsseep(j,k),k=1,NLSF), (dlsshare(j,k),k=1,3)
        read(*,*) dswfl(j,NLEP), dlwfl(j,NLEP), dalb(j), dt(j,NLEP), dt(j,NLEV), dq(j,NLEV), dp(j), du10(j), dv10(j)
       enddo
      endif
"""
STORE_STEP = "       if (nlstore == 1) call lstoreevap\n"
STORE_OUT = ", dlsv(j), dlsfw(j), dlses(j), devap(j), (dlsflx(j,k),k=1,NLSOUT)"


def compile_soil(folder, text, with_store: bool):
    src = twilight_source('plasimmod.f90')
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n      parameter(NLEV_ATM = 15)\n'
                                       '      parameter(NPRO_ATM = 1)\n      end module resmod\n')
    shutil.copy(src, folder / 'plasimmod.f90')
    (folder / 'landsoil.f90').write_text(soil_routines(text))
    (folder / 'stubs.f90').write_text(STUBS)
    harness = HARNESS.replace('{store_input}', STORE_INPUT if with_store else '').replace(
        '{store_step}', STORE_STEP if with_store else '').replace('{store_out}', STORE_OUT if with_store else '')
    (folder / 'harness.f90').write_text(harness)
    files = ['resmod.f90', 'plasimmod.f90', 'landsoil.f90', 'stubs.f90', 'harness.f90']
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def cells(n, seed=7):
    """Land cells as soilstep finds them: dry and wet buckets, some under snow and one a glacier, with rain."""
    rng = np.random.default_rng(seed)
    rows = []
    for j in range(n):
        snow = 0.02 if j % 7 == 3 else 0.0
        rows.append(dict(dls=1.0, dprc=float(rng.choice([0.0, 1e-6, 2e-5])), dprl=float(rng.choice([0.0, 5e-7])),
                         dprs=1e-7 if snow else 0.0, devap=-float(rng.uniform(0.0, 3e-8)),
                         dwatc=float(rng.choice([0.0, 0.002, 0.015, 0.3, 0.4999])), dwmax=0.5, dsnowz=snow,
                         dglac=1.0 if j == 5 else 0.0, dts=float(rng.uniform(295.0, 320.0)) if not snow else 272.0,
                         dtsm=300.0 if not snow else 271.5, dshfl=-float(rng.uniform(0, 120)),
                         dlhfl=-float(rng.uniform(0, 60)), dflux=float(rng.uniform(-50, 500)),
                         dsoilt=list(rng.uniform(290.0, 305.0, 5))))
    return rows


def harness_input(mode, rows, nsteps, store=None):
    c = PLASIM
    lines = [f'{mode} {len(rows)} {nsteps}',
             num([c['ga'], c['gascon'], c['akap'], c['rdbrv'], c['acpd'], c['ra1'], c['ra2'], c['ra4'], DT, 0.99415])]
    for r in rows:
        lines.append(num([r[k] for k in ('dls', 'dprc', 'dprl', 'dprs', 'devap', 'dwatc', 'dwmax', 'dsnowz', 'dglac')]))
        lines.append(num([r['dts'], r['dtsm'], r['dshfl'], r['dlhfl'], r['dflux'], *r['dsoilt']]))
    if store is not None:
        t, st, air = store
        lines.append(num([0.3, PC[9], PC[10], PC[11], PC[12], PC[13], PC[14], PC[15], PC[16], PC[17], PC[18],
                          PLASIM['alv'], PLASIM['tmelt']]))
        lines += [num(t['fill']), num(t['depths']), num(t['wets']), num(t['table'].T.ravel())]
        for j in range(len(rows)):
            lines.append(num([t['cap'][j], st['v'][j], st['fw'][j], *t['area'][j], *t['seep'][j], *t['shares'][j]]))
            lines.append(num([air[k][j] for k in ('sw_net', 'lw_net', 'alb', 'ts', 'tl', 'ql', 'p', 'u10', 'v10')]))
    return '\n'.join(lines) + '\n'


def run(exe, text):
    out = subprocess.run([str(exe)], input=text, capture_output=True, text=True, check=True).stdout
    return [line for line in out.strip('\n').split('\n')]


@FORTRAN
def test_with_the_switch_off_the_patched_soil_step_prints_every_digit_the_old_one_prints(tmp_path):
    source = twilight_source('landmod.f90')
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    old = text if ls.LAND_STORE_MARKER not in text else None
    if old is None:
        pytest.skip('the shared landmod.f90 carries the patch already; compare against a fresh copy')
    rows = cells(14)
    a = run(compile_soil(tmp_path / 'old', old, False), harness_input(0, rows, 6))
    b = run(compile_soil(tmp_path / 'new', patched(old), False), harness_input(0, rows, 6))
    c = run(compile_soil(tmp_path / 'new_store', patched(old), True), harness_input(0, rows, 6))
    assert a == b and [r.split()[:7] for r in c] == [r.split() for r in a]


@FORTRAN
def test_the_fortran_store_gives_the_python_copys_water(tmp_path):
    source = twilight_source('landmod.f90')
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    text = patched(text) if ls.LAND_STORE_MARKER not in text else text
    rows = [r for r in cells(14) if r['dsnowz'] == 0.0 and r['dglac'] == 0.0][:8]
    for r in rows:
        r['dwatc'] = min(r['dwatc'], 0.05)
    n = len(rows)
    t = tables(n)
    rng = np.random.default_rng(13)
    v = np.where(np.arange(n) % 3 == 0, 0.0, rng.uniform(0.0, 0.1, n))
    v[1] = t['cap'][1] * 0.9999
    rows[1]['dprc'] = 4e-4                                                  # a downpour that spills it
    fw, _ = ls.curve(t['fill'], t['area'], t['seep'], t['cap'], v)
    air = dict(sw_net=rng.uniform(0, 600, n), lw_net=rng.uniform(-90, -5, n), alb=np.full(n, 0.2),
               ts=rng.uniform(295, 325, n), tl=rng.uniform(280, 300, n), ql=rng.uniform(0.004, 0.014, n),
               p=rng.uniform(110000, 121000, n), u10=rng.uniform(-3, 3, n), v10=rng.uniform(-3, 3, n))
    nsteps = 40
    out = run(compile_soil(tmp_path / 'store', text, True), harness_input(1, rows, nsteps, (t, dict(v=v, fw=fw), air)))
    sw, lw = ls.radiation_from_plasim(air['sw_net'], air['alb'], air['lw_net'], air['ts'])
    ta = air['tl'] * 0.99415 ** (-PLASIM['akap'])
    rate, _, _ = ls.pool_evaporation(sw, lw, ta, air['ql'], air['p'], np.hypot(air['u10'], air['v10']), PC)
    state = dict(w=np.array([r['dwatc'] for r in rows]), v=v, fw=fw)
    f = dict(prc=np.array([r['dprc'] for r in rows]), prl=np.array([r['dprl'] for r in rows]),
             devap_bucket=np.array([r['devap'] for r in rows]), pool_rate=rate)
    for step in range(nsteps):
        state, flux, _ = ls.land_step(DT, state, f, t, 0.3)
        for j in range(n):
            vals = [float(x) for x in out[step * n + j].split()]
            assert vals[0] == pytest.approx(state['w'][j], rel=1e-10, abs=1e-18)
            assert vals[1] == pytest.approx(flux['drunoff'][j], rel=1e-10, abs=1e-22)
            assert vals[7] == pytest.approx(state['v'][j], rel=1e-10, abs=1e-18)
            assert vals[8] == pytest.approx(state['fw'][j], rel=1e-10, abs=1e-18)
            assert vals[9] == pytest.approx(flux['es'][j], rel=1e-10, abs=1e-22)
            assert vals[10] == pytest.approx(flux['devap'][j], rel=1e-10, abs=1e-22)
            expect = [flux['runoff'][j], state['v'][j], state['fw'][j], -flux['es'][j], flux['seepage'][j],
                      flux['recharge'][j], flux['rivers'][j] + flux['spill'][j]]
            assert vals[11:18] == pytest.approx(expect, rel=1e-10, abs=1e-22)


def inserted(old, new):
    """The lines an edit adds, if it only inserts whole lines (old's lines appear in new in order), else None."""
    a, b = old.splitlines(), new.splitlines()
    out, i = [], 0
    for line in b:
        if i < len(a) and line == a[i]:
            i += 1
        else:
            out.append(line)
    return out if i == len(a) else None


def test_every_edit_only_inserts_lines_and_its_statements_run_only_with_the_switch():
    for old, new in ls.LANDMOD_EDITS:
        assert ls.LAND_STORE_MARKER not in old
        assert all(len(line) <= 132 for line in new.splitlines())
        if old.startswith('     &') and 'dgroundalb' in old:        # the namelist group gains a continuation
            assert new.startswith(old[:-1] + '  &\n') and all(n in new for n in ls.NAMELIST_NAMES)
            continue
        added = inserted(old, new)
        assert added is not None, old[:60]
        if 'subroutine lstoreini' in new or 'integer :: nlstore' in new:
            continue                                                 # new routines and declarations
        statements = [line for line in added if line.strip() and not line.lstrip().startswith('!')]
        if old == ls._BROADCAST:                                     # the namelist values' broadcasts
            assert all(re.fullmatch(r'      call mpbc[ir]\(\w+\)( +!.*)?', s) for s in statements)
            continue
        # every other inserted statement is a switch's test or lies within one
        depth = 0
        for s in statements:
            body = s.split('!')[0].strip()
            if body.startswith('if (nlstore == 1)') or body.startswith('if (nlstore /= 1)'):
                depth += body.endswith('then')
                continue
            assert depth > 0, s
            if body.startswith('if') and body.endswith('then'):
                depth += 1
            elif body == 'endif':
                depth -= 1


def test_the_patch_applies_once_to_the_shared_landmod():
    source = twilight_source('landmod.f90')
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    if ls.LAND_STORE_MARKER not in text:
        for old, _ in ls.LANDMOD_EDITS:
            assert text.count(old) == 1, old[:60]
    text = patched(text)
    assert patched(text) == text
    for routine in ('lstoreini', 'lstoreevap', 'lstorewat', 'lstoreroff', 'lspool', 'lscurve', 'lsstep'):
        assert text.count(f'      subroutine {routine}(') + text.count(f'      subroutine {routine}\n') == 1
    for name in ls.NAMELIST_NAMES:
        assert len(re.findall(rf'call mpbc[ri]\({name}\)', text)) == 1
    assert text.count('call lstorewat') == 1 and text.count('call lstoreroff') == 1 and text.count('call lstoreevap') == 1
    assert [p[0] for p in ls.LAND_STORE_PATCHES] == ['landmod.f90']   # outmod.f90, which every key hashes, is untouched
    # the Fortran's grids and slots are the module's, and its σ is the radiation's own
    radmod = twilight_source('radmod.f90')
    assert re.search(r'parameter\(SBK = 5\.67E-8\)', radmod.read_text())
    assert f'NLSF = {ls.NLSF}, NLSD = {ls.NLSD}, NLSW = {ls.NLSW}' in text and f'NLSOUT = {ls.NLSOUT}' in text


@FORTRAN
def test_the_patched_landmod_compiles_with_plasims_modules(tmp_path):
    source = twilight_source('landmod.f90')
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = source.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    for name in ('plasimmod.f90', 'specblock.f90', 'radmod.f90'):
        shutil.copy(src / name, tmp_path / name)
    (tmp_path / 'landmod.f90').write_text(patched(source.read_text()))
    for name in ('resmod.f90', 'plasimmod.f90', 'specblock.f90', 'radmod.f90', 'landmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


# ------------------------------------------------------------------------------------------------------- runner

def product(tmp_path, cells_=((5, 10), (6, 11))):
    """A synthetic store-tables product on T21: two dryland cells with stores."""
    t = tables(len(cells_))
    rows = []
    for k, (j, i) in enumerate(cells_):
        store = {h: dict(capacity_m=float(t['cap'][k]), area_share=[float(x) for x in t['area'][k]],
                         seepage_m_s={f: [float(x) * (1.0 + q) for x in t['seep'][k]] for q, f in enumerate(ls.FLOORS)})
                 for h in ls.HOLLOWS}
        shares = {b: {s: {h: dict(hollows=0.6, rivers=0.001, recharge=0.002, soil=0.397) for h in ls.HOLLOWS}
                      for s in ls.SURFACES} for b in ls.LIBRARIES}
        rows.append(dict(j=j, i=i, store=store, shares=shares))
    p = dict(schema=ls.PRODUCT_SCHEMA, grid=dict(nlon=64, nlat=32), floors=list(ls.FLOORS), surfaces=list(ls.SURFACES),
             libraries=list(ls.LIBRARIES), hollows=list(ls.HOLLOWS),
             grids=dict(fill=list(t['fill']), rain_depth_mm=list(t['depths']), wetness=list(t['wets'])),
             runoff_tables={s: (t['table'] * (1 + q)).tolist() for q, s in enumerate(ls.SURFACES)},
             producer=dict(files={'research/studies/storm_water/store_tables.py': 'abc'}), basis=dict(run='test'),
             cells=rows)
    path = tmp_path / 'store_tables.json'
    path.write_text(json.dumps(p))
    return path, t


def test_the_runner_records_the_store_and_writes_its_tables_and_namelist(tmp_path):
    path, t = product(tmp_path)
    s = ls.setting(path, 'sandy_banked', 'loam', 'cm1_6km', 'mature', 0.3, G_MOON, 121590.0)
    assert s['floor'] == 'sandy_banked' and s['wetted_fraction'] == 0.3 and 'not settled' in s['wetted_fraction_status']
    assert s['pool']['free_convection_m_s_k3'] == pytest.approx(PC[15]) and s['output']['160'] == 'mrro'
    for bad in (dict(floor='concrete'), dict(surface='rock'), dict(library='earth'), dict(hollows='old'),
                dict(wetted_fraction=0.0)):
        args = dict(floor='clay_floored', surface='loam', library='walnut_gulch', hollows='young', wetted_fraction=0.3)
        args.update(bad)
        with pytest.raises(ValueError):
            ls.setting(path, args['floor'], args['surface'], args['library'], args['hollows'], args['wetted_fraction'],
                       G_MOON, 121590.0)
    nl = ls.namelist_values(s)
    assert nl['NLSTORE'] == '1' and nl['NLSFLOOR'] == '2' and nl['NLSSURF'] == '2' and nl['DLSWET'] == '0.3'
    assert set(k.lower() for k in nl) == set(ls.NAMELIST_NAMES)
    out = ls.write_table(ls.read_product(path), 'cm1_6km', 'mature', tmp_path / 'landstore.dat')
    lines = out.read_text().splitlines()
    assert lines[0].split() == ['1', '2048', '12', '19', '11', '3', '4']
    assert len(lines) == 4 + len(ls.SURFACES) + ls.NUGP
    tab = np.array(lines[4 + 1].split(), float).reshape(ls.NLSW, ls.NLSD).T     # the loam's, depth fastest
    assert tab == pytest.approx(t['table'] * 2.0)
    row = np.array(lines[4 + len(ls.SURFACES) + 5 * 64 + 10].split(), float)
    assert row.size == 1 + 12 * 4 + 12 and row[0] == pytest.approx(t['cap'][0])
    assert row[1 + 12 + 12:1 + 12 + 24] == pytest.approx(t['seep'][0] * 2.0)     # the sandy banks' seepage
    assert row[49:52] == pytest.approx([0.6, 0.001, 0.002])
    assert float(np.abs(np.array(lines[4 + len(ls.SURFACES)].split(), float)).sum()) == 0.0   # no store at the pole
    bad = json.loads(path.read_text())
    bad['schema'] = 'something/1'
    (tmp_path / 'bad.json').write_text(json.dumps(bad))
    with pytest.raises(ValueError):
        ls.read_product(tmp_path / 'bad.json')


def test_the_runner_writes_the_runs_tables_and_refuses_a_changed_product_or_a_store_at_sea(tmp_path):
    path, t = product(tmp_path)
    s = ls.setting(path, 'clay_floored', 'loam', 'walnut_gulch', 'young', 0.3, G_MOON, 121590.0)
    mask = np.zeros((32, 64))
    mask[5, 10] = mask[6, 11] = 1.0
    er._write_sra(tmp_path / 'land.sra', 172, mask)
    out = er.land_store_table(s, tmp_path / 'land.sra', tmp_path / 'inputs' / ls.TABLE_NAME)
    assert out.read_text() == ls.write_table(ls.read_product(path), 'walnut_gulch', 'young', tmp_path / 'again.dat').read_text()
    mask[6, 11] = 0.0                                                   # one store now stands at sea
    er._write_sra(tmp_path / 'sea.sra', 172, mask)
    with pytest.raises(ValueError, match='sea'):
        er.land_store_table(s, tmp_path / 'sea.sra', tmp_path / 'x.dat')
    path.write_text(path.read_text().replace('"basis": {"run": "test"}', '"basis": {"run": "changed"}'))
    with pytest.raises(ValueError, match='changed'):
        er.land_store_table(s, tmp_path / 'land.sra', tmp_path / 'y.dat')


def test_the_runner_wires_the_store_without_touching_other_runs():
    text = (er.HERE / 'exoplasim_run.py').read_text()
    assert "lsf.LAND_STORE_PATCHES" in text and "'NDIAGGP2D', str(lsf.NLSOUT)" in text
    assert "if land_store is not None and not wind10:" in text
    assert "*([inputs['land_store']] if 'land_store' in inputs else [])" in text
    import inspect
    for f in (er.run, er.configuration):
        assert inspect.signature(f).parameters['land_store'].default is None
    assert inspect.signature(er.ensure_patched).parameters['land_store'].default is False
    # the store's codes are the 2-d slots, which no other Terluna field and no PlaSim code takes
    assert set(ls.CODES) == set(range(1, ls.NLSOUT + 1)) and not set(ls.CODES) & set(er.PLASIM_CODES)
    names = [v[0] for v in ls.CODES.values()]
    assert not set(names) & set(er.OUTPUT) and len(set(names)) == len(names)
