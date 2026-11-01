"""Beljaars's gust in PlaSim's surface transfer (climate/gcm/surface_gust.py): its depth is the nonlocal scheme's and its
convective velocity vanishes in stable air; PlaSim's surflx with the patch, compiled with ExoPlaSim's flags against a
stub of PlaSim's module, gives the Python copy's transfer to 1e-10 on lunar and Earth-like columns over water and land,
unstable, calm, windy, stable and neutral, with the gust and with either free-convection constant; with the switch off
it prints every digit the unpatched surflx prints; the patch applies once over the nonlocal scheme's; the patched
fluxmod.f90 compiles with PlaSim's modules; and the runner records the gust and writes its namelist."""
import inspect
import re
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import levels as lv
from climate.gcm import pbl_fortran as pbf
from climate.gcm import pbl_scheme as ps
from climate.gcm import surface_gust as sg
from climate.tests.test_pbl_fortran import (ALV, EARTH, FORTRAN, GCM_FLAGS, LUNAR, apply, num, twilight_fluxmod, undo,
                                            with_free_convection)


def without_gust(text):
    """The shared source as it stood before this patch: its edits undone where it carries them, the similarity over
    water's, which goes over them, first."""
    from climate.gcm import surface_similarity as ssf
    return undo(undo(text, ssf.SIMW_EDITS, ssf.SIMW_MARKER), sg.GUST_EDITS, sg.GUST_MARKER)


def with_gust(text):
    text = apply(with_free_convection(without_gust(text)), pbf.PBL_EDITS, pbf.PBL_MARKER)
    return apply(text, sg.GUST_EDITS, sg.GUST_MARKER)


def surflx(text, last='surflx'):
    start = text.index('      subroutine surflx\n')
    end = text.index(f'      end subroutine {last}\n') + len(f'      end subroutine {last}\n')
    return text[start:end]


# ------------------------------------------------------------------------------------------------------- columns

def column(c, half, ps_pa, theta, zi, lapse, q_bl, wind, dts, dls, dz0, shf, lhf, ust, seed):
    """One of PlaSim's columns as surflx finds it: a mixed layer zi deep at potential temperature theta, stable lapse
    K/m above, vapour q_bl in the mixed layer falling above, the wind at the lowest level wind m/s, the surface dts K
    warmer than the lowest level's air (in potential temperature), and the surface fluxes of the step before: sensible
    and latent heat shf and lhf (W/m2, upward) and u* ust."""
    rng = np.random.default_rng(seed)
    sigma = ps.full_sigma(half)
    akap = c['gascon'] / c['acpd']
    z = -c['gascon'] * 280.0 / c['ga'] * np.log(sigma)
    th = np.where(z < zi, theta, theta + lapse * (z - zi)) + rng.normal(0.0, 0.05, z.size)
    t = th * sigma ** akap
    q = np.where(z < zi, q_bl, q_bl * np.exp(-(z - zi) / 2500.0)) * (1.0 + rng.normal(0.0, 0.01, z.size))
    u = wind + 1.0e-3 * z + rng.normal(0.0, 0.3, z.size)
    u[-1] = wind
    v = rng.normal(0.0, 0.8, z.size)
    v[-1] = 0.0
    ts = th[-1] + dts                                   # theta is referred to the ground, where sigma is 1
    qs = 0.622 * 610.78 * np.exp(17.27 * (ts - 273.16) / (ts - 35.86)) / ps_pa
    rho = ps_pa / (c['gascon'] * ts)
    return dict(dt=np.r_[t, ts], dq=np.r_[q, qs if dls < 1.0 else 0.5 * qs], du=u, dv=v, dz0=dz0, dls=dls, dp=ps_pa,
                dshfl=-shf * max(1.0, 1.0 + c['ADV'] * q[-1]), devap=-lhf / (ALV * 1000.0), dust3=ust ** 3, rho=rho)


def cases():
    """(constants, half levels, columns): lunar on fifteen levels and on forty equal ones, and Earth-like on PlaSim's
    ten; over the sea unstable, calm, windy, stable and with last step's fluxes downward though the air is unstable;
    over land a convective day, a calm stable night with the wind below every floor, and neutral air."""
    half15 = np.asarray(lv.half_levels('l15'), float)
    half10 = np.asarray(lv.half_levels(lv.DEFAULT), float)
    equal40 = np.linspace(1.0, 40.0, 40) / 40.0
    sets = []
    for c, half, ps_pa, zi in ((LUNAR, half15, 121590.0, 4000.0), (LUNAR, equal40, 121590.0, 3000.0),
                               (EARTH, half10, 101325.0, 900.0)):
        rows = [  # theta, q, wind, dts, dls, dz0, shf, lhf, u*
            (299.0, 0.012, 3.0, 2.7, 0.0, 2.0e-4, 9.0, 110.0, 0.10),
            (299.0, 0.013, 0.2, 3.0, 0.0, 1.5e-5, 6.0, 60.0, 0.02),
            (299.0, 0.011, 9.0, 1.5, 0.0, 1.0e-3, 12.0, 220.0, 0.35),
            (301.0, 0.014, 1.5, -1.0, 0.0, 1.0e-4, -3.0, 5.0, 0.05),
            (299.0, 0.012, 0.6, 1.0, 0.0, 1.0e-4, -2.0, -1.0, 0.03),
            (300.0, 0.010, 2.0, 6.0, 1.0, 0.1, 40.0, 150.0, 0.25),
            (300.0, 0.006, 0.05, -5.0, 1.0, 0.1, -5.0, 2.0, 0.01),
            (300.0, 0.008, 4.0, 0.0, 1.0, 0.05, 0.0, 30.0, 0.2),
        ]
        cols = [column(c, half, ps_pa + 200.0 * j, th, zi * (0.7 + 0.1 * j), 1.0e-3, q, w, dts, dls, dz0, shf, lhf, ust,
                       seed=j) for j, (th, q, w, dts, dls, dz0, shf, lhf, ust) in enumerate(rows)]
        sets.append((c, half, cols))
    return sets


STUB = """
      module fluxmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLEM={nlev}-1, NLON={nhor}, NLAT=1, NROOT=0
      real, parameter :: vonkarman=0.4
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), du(NHOR,NLEP), dv(NHOR,NLEP)
      real :: dtsa(NHOR), dz0(NHOR), dls(NHOR), dp(NHOR), dshfl(NHOR), devap(NHOR), dust3(NHOR)
      real :: du10(NHOR), dv10(NHOR), dtransh(NHOR), dtransm(NHOR)
      real :: sigma(NLEV), sigmah(NLEV)
      real :: ga, gascon, akap, rdbrv, acpd, adv
      integer :: ntsa=2, nwind10=0, nstress=0, nshfl=0, nevap=0, ntime=0, nprint=0, nprhor=1, nud=6, mypid=0
      real :: time4tr=0., time4st=0., time4sh=0., time4ev=0.
      real :: zumin=1., vdiff_lamm=160., vdiff_b=5., vdiff_c=5., vdiff_d=5., cfree=0.0016
{declarations}      end module fluxmod

      subroutine mksecond(a,b)
      real a,b
      end subroutine mksecond
      subroutine mpgagp(a,b,n)
      real a(*),b(*)
      integer n
      end subroutine mpgagp
      subroutine mkstress
      end subroutine mkstress
      subroutine mkshfl
      end subroutine mkshfl
      subroutine mkevap
      end subroutine mkevap
"""
HARNESS = """
      program harness
      use fluxmod
      implicit none
      integer :: jhor, jlev
      read(*,*) ga, gascon, akap, rdbrv, acpd, adv
      read(*,*) sigma, sigmah
      read(*,*) zumin, cfree{switch}
      do jhor=1,NHOR
       read(*,*) dp(jhor), dshfl(jhor), devap(jhor), dust3(jhor), dz0(jhor), dls(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEP)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEP)
       read(*,*) (du(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dv(jhor,jlev), jlev=1,NLEV)
      enddo
      call surflx
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') dtransh(jhor), dtransm(jhor), dtsa(jhor)
      enddo
      end program harness
"""


def compile_surflx(folder, code, nhor, nlev, with_gust: bool):
    folder.mkdir(parents=True, exist_ok=True)
    decl = sg.declarations() + pbf.declarations() if with_gust else ''
    files = {'stub.f90': STUB.format(nhor=nhor, nlev=nlev, declarations=decl), 'surflx.f90': code,
             'harness.f90': HARNESS.replace('{switch}', ', ngust, gust_beta, gust_umin' if with_gust else '')}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def harness_input(c, half, cols, cfree, ngust=None):
    sigma = ps.full_sigma(half)
    akap = c['gascon'] / c['acpd']
    switch = '' if ngust is None else f' {int(ngust)} {sg.GUST_BETA!r} {sg.GUST_UMIN!r}'
    lines = [num([c['ga'], c['gascon'], akap, c['rdbrv'], c['acpd'], c['ADV']]), num(np.r_[sigma, half]),
             f'{sg.EXOPLASIM_ZUMIN!r} {cfree!r}' + switch]
    for x in cols:
        lines.append(num([x['dp'], x['dshfl'], x['devap'], x['dust3'], x['dz0'], x['dls']]))
        lines += [num(x[k]) for k in ('dt', 'dq', 'du', 'dv')]
    return '\n'.join(lines) + '\n'


def run_harness(exe, text, ncols):
    rows = subprocess.run([str(exe)], input=text, capture_output=True, text=True, check=True).stdout.strip('\n')
    rows = rows.split('\n')
    assert len(rows) == ncols
    return rows


# ------------------------------------------------------------------------------------------------------- tests

def test_the_gust_takes_the_nonlocal_schemes_depth_and_vanishes_in_stable_air():
    for c, half, cols in cases():
        n = half.size
        sigma = ps.full_sigma(half)
        for x in cols:
            s = pbf.surface_from_plasim(x['dshfl'], x['devap'], x['dust3'], x['dt'][n], x['dp'], x['dq'][n - 1], c)
            g = sg.gust_column(x['dt'][:n], x['dq'][:n], x['du'], x['dv'], sigma, half, s['wth0'], s['wq0'], s['ustar'], c)
            nl = ps.nonlocal_column(x['dt'][:n], x['dq'][:n], x['du'], x['dv'], sigma, half, s['wth0'], s['wq0'],
                                    s['ustar'], np.zeros(n - 1), np.zeros(n - 1), np.ones(n - 1), 1800.0, c)
            assert g['unstable'] == nl['unstable']
            if g['unstable']:
                assert g['h'] == pytest.approx(nl['h'], rel=1e-12) and g['wstar'] == pytest.approx(nl['wstar'], rel=1e-12)
                assert g['wstar'] == pytest.approx((c['ga'] / (x['dt'][n - 1] / sigma[-1] ** (c['gascon'] / c['acpd'])
                                                               * (1.0 + (1.0 / c['rdbrv'] - 1.0) * x['dq'][n - 1]))
                                                    * g['wthv0'] * g['h']) ** (1.0 / 3.0), rel=1e-12)
            else:
                assert g['h'] == 0.0 and g['wstar'] == 0.0


def test_over_the_lunar_sea_the_gust_is_small_and_the_tank_term_goes():
    c, half, cols = cases()[0]
    sea = cols[0]                                                     # 3 m/s, 2.7 K, 9 and 110 W/m2 the step before
    earth, moon = (sg.surface_transfer(sea, ps.full_sigma(half), half, c, cfree=k) for k in (0.0016, 0.00083))
    gust = sg.surface_transfer(sea, ps.full_sigma(half), half, c, ngust=1)
    assert 0.3 < gust['wstar'] < 1.0 and 1000.0 < gust['h'] < 8000.0
    assert gust['absu2'] == pytest.approx(9.0 + gust['wstar'] ** 2)
    neutral = np.sqrt(9.0) * (ps.KARMAN / np.log(earth['znl'] / sea['dz0'] + 1.0)) ** 2
    # the transfer of heat against the neutral: the tank term's at Earth's and the Moon's constants, the gust's
    assert earth['dtransh'] / neutral > moon['dtransh'] / neutral > 1.25 > gust['dtransh'] / neutral > 1.0
    # over land the transfer was Louis's already; with the gust it takes the gust's wind alone
    land = cols[5]
    plain, gusty = (sg.surface_transfer(land, ps.full_sigma(half), half, c, ngust=k) for k in (0, 1))
    assert gusty['dtransh'] > plain['dtransh'] and gusty['wstar'] > 0.0
    # a calm stable night: 1 m/s without the gust, 0.1 m/s with it
    night = cols[6]
    assert sg.surface_transfer(night, ps.full_sigma(half), half, c)['absu2'] == 1.0
    assert sg.surface_transfer(night, ps.full_sigma(half), half, c, ngust=1)['absu2'] == pytest.approx(0.01)


@FORTRAN
def test_the_fortran_surflx_gives_the_python_copys_transfer(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    code = surflx(with_gust(source.read_text()), last='gustcol')
    for i, (c, half, cols) in enumerate(cases()):
        exe = compile_surflx(tmp_path / f'case{i}', code, len(cols), half.size, True)
        for cfree, ngust in ((0.0016, 0), (0.00083, 0), (0.0016, 1)):
            rows = run_harness(exe, harness_input(c, half, cols, cfree, ngust), len(cols))
            for x, row in zip(cols, rows):
                h, m, tsa = (float(v) for v in row.split())
                py = sg.surface_transfer(x, ps.full_sigma(half), half, c, ngust=ngust, cfree=cfree)
                assert h == pytest.approx(py['dtransh'], rel=1e-10) and m == pytest.approx(py['dtransm'], rel=1e-10)
                assert tsa == pytest.approx(py['dtsa'], rel=1e-12)


@FORTRAN
def test_with_the_switch_off_the_patched_surflx_prints_every_digit_the_old_one_prints(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    old, new = surflx(with_free_convection(without_gust(text))), surflx(with_gust(text), last='gustcol')
    for i, (c, half, cols) in enumerate(cases()):
        exe_old = compile_surflx(tmp_path / f'old{i}', old, len(cols), half.size, False)
        exe_new = compile_surflx(tmp_path / f'new{i}', new, len(cols), half.size, True)
        for cfree in (0.0016, 0.00083):
            a = run_harness(exe_old, harness_input(c, half, cols, cfree), len(cols))
            b = run_harness(exe_new, harness_input(c, half, cols, cfree, 0), len(cols))
            assert a == b


def test_the_patch_applies_once_over_the_boundary_layers():
    assert not any(sg.GUST_MARKER in old for old, _ in sg.GUST_EDITS)
    assert all(len(line) <= 132 for _, new in sg.GUST_EDITS for line in new.splitlines())
    # its module edits stand on the free-convection patch's line and the nonlocal scheme's namelist
    assert sg.GUST_EDITS[0][0] == pbf.PBL_EDITS[0][0] and sg.GUST_EDITS[1][0] == pbf.namelist_text()
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = with_gust(source.read_text())
    assert apply(text, sg.GUST_EDITS, sg.GUST_MARKER) == text and without_gust(text) != text
    assert text.count('       call gustnl(zgust2)') == 1 and text.count('      subroutine gustcol(') == 1
    assert text.count('if(dls(jhor) < 1. .and. ngust /= 1) then') == 1
    for name in sg.NAMELIST_NAMES:
        assert len(re.findall(rf'(real|integer) :: {name} = ', text)) == 1
        assert len(re.findall(rf'call mpbc[ri]\({name}\)', text)) == 1
        assert len(re.findall(rf'&\s+,.*\b{name}\b', text)) == 1


@FORTRAN
def test_the_patched_fluxmod_compiles_with_plasims_modules(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = source.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(src / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'fluxmod.f90').write_text(with_gust(source.read_text()))
    for name in ('resmod.f90', 'plasimmod.f90', 'fluxmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def test_the_runner_records_the_gust_and_writes_its_namelist():
    for f in (er.run, er.configuration):
        assert inspect.signature(f).parameters['surface_gust'].default is None
    s = er.surface_gust_setting('beljaars')
    assert s['beta'] == sg.GUST_BETA == 1.0 and s['least_wind_m_s'] == sg.GUST_UMIN == 0.1
    assert s['depth']['parameters'] == pbf.fortran_values()
    for bad in ('earth', 'coare'):
        with pytest.raises(ValueError):
            er.surface_gust_setting(bad)
    text = (er.HERE / 'exoplasim_run.py').read_text()
    assert "'fluxmod_namelist', 'NGUST', '1'" in text
