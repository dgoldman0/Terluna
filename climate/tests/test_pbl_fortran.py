"""The nonlocal boundary layer's Fortran for PlaSim (climate/gcm/pbl_fortran.py): the patch applies once to the shared
twilight source over the free-convection patch; PlaSim's vdiff with it, compiled with ExoPlaSim's flags against a stub of
PlaSim's module, gives the Python copy's tendencies (pbl_scheme.vdiff) to 1e-10 on lunar and Earth-like columns,
unstable, calm, stable and neutral, with either mixing length; it keeps the column's heat, water and momentum; with
the switch off it prints every digit the unpatched vdiff prints; the patched fluxmod.f90 and rainmod.f90 compile with
PlaSim's modules on fifteen levels; the runner records the options and writes their namelist; and the Python dry
adjustment (pbl_scheme.dca) gives the tendencies of PlaSim's own mkdca."""
import math
import re
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import levels as lv
from climate.gcm import pbl_fortran as pbf
from climate.gcm import pbl_scheme as ps

GCM_FLAGS = ('-O3', '-cpp', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffpe-summary=none',
             '-finit-real=zero', '-fdefault-real-8')
FORTRAN = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')
G_MOON = 1.6242188765608128
LUNAR = dict(ga=G_MOON, gascon=288.3767, acpd=288.3767 / 0.286, rdbrv=288.3767 / 461.51, ADV=0.8546)
EARTH = dict(ga=9.80665, gascon=287.05, acpd=1004.6, rdbrv=287.05 / 461.51, ADV=0.8608)
ALV = 2.5008e6
num = lambda xs: ' '.join(repr(float(x)) for x in np.atleast_1d(xs))


def twilight_fluxmod():
    path = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'fluxmod.f90'
    return path if path.exists() else None


def apply(text, edits, marker):
    if marker in text:
        return text
    for old, new in edits:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    return text


def with_free_convection(text):
    return apply(text, er.FREE_CONVECTION_EDITS, er.FREE_CONVECTION_MARKER)


def undo(text, edits, marker):
    """The text as it stood before a patch: its edits undone, the last first, where the text carries it."""
    if marker not in text:
        return text
    for old, new in reversed(edits):
        assert text.count(new) == 1, new[:60]
        text = text.replace(new, old)
    return text


def before_later_patches(text):
    """The shared source as the free-convection patch left it: the similarity over water and the surface gust's patch,
    which go over this one, and this one undone where the source carries them."""
    from climate.gcm import surface_gust as sg
    from climate.gcm import surface_similarity as ssf
    text = undo(text, ssf.SIMW_EDITS, ssf.SIMW_MARKER)
    return undo(undo(text, sg.GUST_EDITS, sg.GUST_MARKER), pbf.PBL_EDITS, pbf.PBL_MARKER)


def unpatched(text):
    """The shared source as it stood before this patch: the free-convection patch on, this one off."""
    return with_free_convection(before_later_patches(text))


def patched(text):
    return apply(with_free_convection(text), pbf.PBL_EDITS, pbf.PBL_MARKER)


def routines(text, last='vdiff'):
    start = text.index('      subroutine vdiff\n')
    end = text.index(f'      end subroutine {last}\n') + len(f'      end subroutine {last}\n')
    return text[start:end]


# ------------------------------------------------------------------------------------------------------- columns

def sounding(c, ps_pa, half, theta, zi, lapse, q_bl, q_scale, wind, shear, seed, noise=0.05):
    """A column top first on half levels half: a mixed layer zi deep at potential temperature theta (referred to the
    surface), stable lapse K/m above, vapour q_bl falling with scale q_scale above it, wind wind + shear z, noise."""
    rng = np.random.default_rng(seed)
    sigma = ps.full_sigma(half)
    akap = c['gascon'] / c['acpd']
    z = -c['gascon'] * 280.0 / c['ga'] * np.log(sigma)
    th = np.where(z < zi, theta, theta + lapse * (z - zi)) + rng.normal(0.0, noise, z.size)
    t = np.maximum(th * sigma ** akap, 160.0)
    q = np.where(z < zi, q_bl, q_bl * np.exp(-(z - zi) / q_scale)) * (1.0 + rng.normal(0.0, 0.01, z.size))
    u = wind + shear * z / 1000.0 + rng.normal(0.0, 0.5, z.size)
    v = rng.normal(0.0, 1.5, z.size)
    return t, np.maximum(q, 1.0e-9), u, v


def cases():
    """(constants, half levels, surface pressure, columns): each column PlaSim's state (dt, dq, du, dv and the
    tendencies so far) and its surface fluxes as surflx leaves them (dshfl and devap downward, dust3), varied."""
    half15 = np.asarray(lv.half_levels('l15'), float)
    half10 = np.asarray(lv.half_levels(lv.DEFAULT), float)
    equal40 = np.linspace(1.0, 40.0, 40) / 40.0
    sets = []
    rng = np.random.default_rng(7)
    for c, half, ps_pa, kind in ((LUNAR, half15, 121500.0, 'lunar'), (LUNAR, equal40, 121500.0, 'lunar'),
                                 (EARTH, half10, 101325.0, 'earth')):
        cols = []
        zi = 4000.0 if kind == 'lunar' else 900.0
        lapse = 1.0e-3 if kind == 'lunar' else 4.0e-3
        q_scale = 3000.0 if kind == 'lunar' else 1500.0
        # (sensible, latent heat W/m2 upward, u*): unstable, calm, windy, stable, neutral, a dry column moistened hard
        for j, (shf, lhf, ust, q_bl) in enumerate([(8.0, 100.0, 0.1, 0.012), (12.0, 150.0, 0.0, 0.014),
                                                   (5.0, 80.0, 0.6, 0.010), (-20.0, -2.0, 0.25, 0.008),
                                                   (0.0, 0.0, 0.2, 0.010), (15.0, 300.0, 0.1, 1.0e-9),
                                                   (3.0, 200.0, 0.05, 0.016), (-4.0, 30.0, 0.15, 0.011)]):
            t, q, u, v = sounding(c, ps_pa + 300.0 * j, half, 300.0 - j, zi * (0.6 + 0.15 * j), lapse,
                                  q_bl, q_scale, 3.0 + j, 1.0e-3 * (j - 3), seed=j)
            n = t.size
            ts = t[-1] / ps.full_sigma(half)[-1] ** (c['gascon'] / c['acpd']) + shf / 10.0
            rho = (ps_pa + 300.0 * j) / (c['gascon'] * ts)
            dq_old = np.maximum(q - rng.normal(0.0, 1.0e-6, n), 1.0e-9)
            cols.append(dict(dp=ps_pa + 300.0 * j, dt=np.r_[t, ts], dq=dq_old,
                             du=u, dv=v, dtdt=rng.normal(0.0, 2.0e-6, n), dqdt=(q - dq_old) / 3600.0,
                             dudt=rng.normal(0.0, 1.0e-5, n), dvdt=rng.normal(0.0, 1.0e-5, n),
                             dshfl=-shf * max(1.0, 1.0 + c['ADV'] * dq_old[-1]), devap=-lhf / (ALV * 1000.0),
                             dust3=ust ** 3, rho=rho))
        sets.append((c, half, cols))
    return sets


STUB = """
      module fluxmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLEM={nlev}-1
      real, parameter :: vonkarman=0.4
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), du(NHOR,NLEP), dv(NHOR,NLEP)
      real :: dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dudt(NHOR,NLEP), dvdt(NHOR,NLEP)
      real :: dp(NHOR), dshfl(NHOR), devap(NHOR), dust3(NHOR)
      real :: sigma(NLEV), sigmah(NLEV), dsigma(NLEV)
      real :: ga, gascon, akap, rdbrv, acpd, adv, deltsec2
      real :: dgp3d(NHOR,NLEV,21)
      real :: dentropy(NHOR,40), dentro(NHOR), dentrot(NHOR,NLEV), dentroq(NHOR,NLEV), dentrop(NHOR)
      real :: dentro3d(NHOR,NLEV,40), denergy(NHOR,40), dener3d(NHOR,NLEV,40)
      integer :: ndheat=1, ndiaggp=1, nentropy=0, nentro3d=0, nenergy=0, nener3d=0
      real :: zumin=1., vdiff_lamm=160., vdiff_b=5., vdiff_c=5., vdiff_d=5.
{declarations}      end module fluxmod
"""
HARNESS = """
      program harness
      use fluxmod
      implicit none
      integer :: jhor, jlev
      real :: zt0(NHOR,NLEP), zq0(NHOR,NLEP), zu0(NHOR,NLEP), zv0(NHOR,NLEP)
      read(*,*) ga, gascon, akap, rdbrv, acpd, adv, deltsec2
      read(*,*) sigma, sigmah, dsigma
      read(*,*) vdiff_lamm, ndheat{switch}
      do jhor=1,NHOR
       read(*,*) dp(jhor), dshfl(jhor), devap(jhor), dust3(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEP)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEV)
       read(*,*) (du(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dv(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dtdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dudt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dvdt(jhor,jlev), jlev=1,NLEV)
      enddo
      zt0(:,:)=dtdt(:,:)
      zq0(:,:)=dqdt(:,:)
      zu0(:,:)=dudt(:,:)
      zv0(:,:)=dvdt(:,:)
      call vdiff
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') (dtdt(jhor,jlev)-zt0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dqdt(jhor,jlev)-zq0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dudt(jhor,jlev)-zu0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dvdt(jhor,jlev)-zv0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dgp3d(jhor,jlev,4), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dgp3d(jhor,jlev,14), jlev=1,NLEV)
      enddo
      end program harness
"""
OUTPUTS = ('dtdt', 'dqdt', 'dudt', 'dvdt', 'dgp4', 'dgp14')


def compile_vdiff(folder, code, nhor, nlev, with_switch: bool):
    folder.mkdir(parents=True, exist_ok=True)
    decl = pbf.declarations() if with_switch else ''
    files = {'stub.f90': STUB.format(nhor=nhor, nlev=nlev, declarations=decl), 'vdiff.f90': code,
             'harness.f90': HARNESS.replace('{switch}', ', nlpbl' if with_switch else '')}
    for name, text in files.items():
        (folder / name).write_text(text)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    return folder / 'harness'


def harness_input(c, half, cols, dtime, lamm, ndheat, nlpbl=None):
    sigma = ps.full_sigma(half)
    akap = c['gascon'] / c['acpd']
    lines = [num([c['ga'], c['gascon'], akap, c['rdbrv'], c['acpd'], c['ADV'], dtime]),
             num(np.r_[sigma, half, np.diff(np.r_[0.0, half])]),
             f"{lamm!r} {int(ndheat)}" + ('' if nlpbl is None else f' {int(nlpbl)}')]
    for x in cols:
        lines.append(num([x['dp'], x['dshfl'], x['devap'], x['dust3']]))
        lines += [num(x[k]) for k in ('dt', 'dq', 'du', 'dv', 'dtdt', 'dqdt', 'dudt', 'dvdt')]
    return '\n'.join(lines) + '\n'


def run_harness(exe, text, ncols):
    rows = subprocess.run([str(exe)], input=text, capture_output=True, text=True, check=True).stdout.strip('\n')
    rows = rows.split('\n')
    assert len(rows) == 6 * ncols
    return rows


def parse(rows, ncols):
    return [{k: np.array(rows[6 * j + i].split(), float) for i, k in enumerate(OUTPUTS)} for j in range(ncols)]


def python_vdiff(c, half, x, dtime, lamm, ndheat, nlpbl):
    n = half.size
    t = x['dt'][:n] + x['dtdt'] * dtime
    q = x['dq'][:n] + x['dqdt'] * dtime
    u = x['du'] + x['dudt'] * dtime
    v = x['dv'] + x['dvdt'] * dtime
    surface = None
    if nlpbl:
        surface = pbf.surface_from_plasim(x['dshfl'], x['devap'], x['dust3'], x['dt'][n], x['dp'], x['dq'][n - 1], c)
    r = ps.vdiff(t, q, u, v, ps.full_sigma(half), half, dtime, c, lamm=lamm, q_heat=x['dq'][:n], ndheat=ndheat,
                 surface=surface)
    dsig = np.diff(np.r_[0.0, half])
    return dict(dtdt=r['dtdt'], dqdt=r['dqdt'], dudt=r['dudt'], dvdt=r['dvdt'], dgp4=r['dtdt_diffusion'],
                dgp14=r['dqdt'] * x['dp'] * dsig / c['ga'] / 1000.0, pbl=r['pbl'])


@FORTRAN
def test_the_fortran_vdiff_gives_the_python_copys_tendencies(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    code = routines(patched(source.read_text()), last='pblcol')
    seen = set()
    for i, (c, half, cols) in enumerate(cases()):
        exe = compile_vdiff(tmp_path / f'set{i}', code, len(cols), half.size, True)
        for lamm in (ps.LAMBDA_EARTH_M, ps.lunar_length(c['ga'])):
            for nlpbl in (0, 1):
                out = parse(run_harness(exe, harness_input(c, half, cols, 3600.0, lamm, True, nlpbl), len(cols)),
                            len(cols))
                py = [python_vdiff(c, half, x, 3600.0, lamm, True, nlpbl) for x in cols]
                # each output's scale over the set, so a column whose tendency is round-off about nil is held to it
                scale = {k: max(np.abs(np.asarray(p[k], float)).max() for p in py) for k in OUTPUTS}
                for j, (x, f, p) in enumerate(zip(cols, out, py)):
                    for k in OUTPUTS:
                        b = np.asarray(p[k], float)
                        assert np.allclose(f[k], b, rtol=1e-10, atol=1e-10 * scale[k]), (i, lamm, nlpbl, j, k)
                    if nlpbl:
                        pb = p['pbl']
                        seen.add((pb['unstable'], pb['ustar'] > 0.0, pb['vapour_kept'],
                                  pb['inversion_interface'] is not None))
                        # with the switch in unstable air the scheme acts: the tendencies differ from the local
                        # diffusion's, the vapour's unless the explicit step would have left it negative
                        if pb['unstable']:
                            local = python_vdiff(c, half, x, 3600.0, lamm, True, 0)
                            assert not np.allclose(local['dtdt'], p['dtdt'], rtol=1e-6, atol=0.0), (i, j)
                            if not pb['vapour_kept']:
                                assert not np.allclose(local['dqdt'], p['dqdt'], rtol=1e-6, atol=0.0), (i, j)
    # the columns reach unstable air with and without u*, stable air, an entrainment flux and a kept vapour profile
    assert {(True, True, False, True), (True, False, False, True), (False, True, False, False)} <= seen
    assert any(s[2] for s in seen)


@FORTRAN
def test_the_fortran_vdiff_keeps_the_columns_heat_water_and_momentum(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    code = routines(patched(source.read_text()), last='pblcol')
    c, half, cols = cases()[0]
    exe = compile_vdiff(tmp_path / 'keep', code, len(cols), half.size, True)
    out = parse(run_harness(exe, harness_input(c, half, cols, 3600.0, ps.lunar_length(c['ga']), True, 1), len(cols)),
                len(cols))
    dsig = np.diff(np.r_[0.0, half])
    n = half.size
    for x, f in zip(cols, out):
        # the column keeps its sums of sigma-thickness times vapour, temperature and wind to round-off of their content
        for k, content in (('dqdt', x['dq'][:n]), ('dgp4', x['dt'][:n]), ('dudt', x['du']), ('dvdt', x['dv'])):
            assert abs(np.sum(dsig * f[k])) <= 1e-12 * np.sum(dsig * np.abs(content)) / 3600.0, k


@FORTRAN
def test_with_the_switch_off_the_patched_vdiff_prints_every_digit_the_old_one_prints(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    old, new = routines(unpatched(text)), routines(patched(text), last='pblcol')
    for i, (c, half, cols) in enumerate(cases()):
        exe_old = compile_vdiff(tmp_path / f'old{i}', old, len(cols), half.size, False)
        exe_new = compile_vdiff(tmp_path / f'new{i}', new, len(cols), half.size, True)
        for lamm in (ps.LAMBDA_EARTH_M, ps.lunar_length(c['ga'])):
            a = run_harness(exe_old, harness_input(c, half, cols, 3600.0, lamm, True), len(cols))
            b = run_harness(exe_new, harness_input(c, half, cols, 3600.0, lamm, True, 0), len(cols))
            assert a == b


def test_the_patch_applies_once_over_the_free_convection_patch():
    source = twilight_fluxmod()
    assert not any(pbf.PBL_MARKER in old for old, _ in pbf.PBL_EDITS)
    assert not any(pbf.DIAGNOSTIC_MARKER in old for old, _ in pbf.DIAGNOSTIC_EDITS)
    assert all(len(line) <= 132 for _, new in pbf.PBL_EDITS + pbf.DIAGNOSTIC_EDITS for line in new.splitlines())
    # every edit adds text after what it finds, but for the namelist line it continues
    changed = [old for old, new in pbf.PBL_EDITS + pbf.DIAGNOSTIC_EDITS if not new.startswith(old)]
    assert changed == ['     &                ,zumin,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d,cfree\n']
    # the module edit's anchor is the line the free-convection patch writes
    assert pbf.PBL_EDITS[0][0] == er.FREE_CONVECTION_EDITS[0][1].split('\n')[1] + '\n'
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = patched(source.read_text())
    assert apply(text, pbf.PBL_EDITS, pbf.PBL_MARKER) == text
    assert text.count('       call pblnl(zt,zq,zu,zv,zkdiffm,zkdiffh)') == 1
    assert text.count('      subroutine pblcol(') == 1 and text.count('      end subroutine pblcol') == 1
    for name in pbf.namelist_names():
        assert len(re.findall(rf'(real|integer) :: {name} = ', text)) == 1
        assert len(re.findall(rf'call mpbc[ri]\({name}\)', text)) == 1
        assert len(re.findall(rf'&\s+,.*\b{name}\b', text)) == 1


@FORTRAN
def test_the_patched_sources_compile_with_plasims_modules(tmp_path):
    source = twilight_fluxmod()
    if source is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = source.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(src / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'fluxmod.f90').write_text(patched(source.read_text()))
    rain = (src / 'rainmod.f90').read_text()
    from climate.gcm import massflux_fortran as mff
    rain = apply(apply(rain, er.SBM_EDITS, er.SBM_MARKER), mff.MASSFLUX_EDITS, mff.MASSFLUX_MARKER)
    rain = apply(rain, pbf.DIAGNOSTIC_EDITS, pbf.DIAGNOSTIC_MARKER)
    assert rain.count(pbf.DIAGNOSTIC_MARKER) == 2
    (tmp_path / 'rainmod.f90').write_text(rain)
    for name in ('resmod.f90', 'plasimmod.f90', 'fluxmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def test_the_runner_records_the_options_and_their_namelist():
    import inspect
    for f in (er.run, er.configuration):
        params = inspect.signature(f).parameters
        assert params['mixing_length'].default is None and params['boundary_layer'].default is None
    m = er.mixing_length_setting('lunar', G_MOON)
    assert m['momentum_m'] == pytest.approx(160.0 * 9.80665 / G_MOON) and m['heat_m'] == pytest.approx(
        m['momentum_m'] * math.sqrt(7.5))
    b = er.boundary_layer_setting('nonlocal')
    assert b['parameters'] == pbf.fortran_values() and set(b['parameters']) == set(pbf.namelist_names()[1:])
    for bad in (('earth', G_MOON), ('mars', G_MOON), ('lunar', 0.0)):
        with pytest.raises(ValueError):
            er.mixing_length_setting(*bad)
    for bad in ('local', 'eddy'):
        with pytest.raises(ValueError):
            er.boundary_layer_setting(bad)
    text = (er.HERE / 'exoplasim_run.py').read_text()
    assert "'fluxmod_namelist', 'VDIFF_LAMM'" in text and "'fluxmod_namelist', 'NLPBL', '1'" in text
    # the vapour by process, written only by runs with the options
    assert set(er.BUDGET) == {30, 33, 34, 35, 36, 37, 39} and not set(er.BUDGET) & (set(er.HEATING) | set(er.MOISTENING))


DCA_STUB = """
      module rainmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLEM={nlev}-1
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV)
      real :: akap, deltsec2, ga, acpd, adv
      real :: dgp3d(NHOR,NLEV,21)
      real :: dentropy(NHOR,40), dentro(NHOR), dentrot(NHOR,NLEV), dentroq(NHOR,NLEV), dentrop(NHOR)
      real :: dentro3d(NHOR,NLEV,40), denergy(NHOR,40), dener3d(NHOR,NLEV,40)
      integer :: ndiaggp=1, nentropy=0, nentro3d=0, nenergy=0, nener3d=0
      end module rainmod
"""
DCA_HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: jhor, jlev
      real :: zt0(NHOR,NLEP), zq0(NHOR,NLEP)
      read(*,*) akap, deltsec2, ga, acpd, adv
      read(*,*) sigma, dsigma
      do jhor=1,NHOR
       read(*,*) dp(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dtdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqdt(jhor,jlev), jlev=1,NLEV)
      enddo
      zt0(:,:)=dtdt(:,:)
      zq0(:,:)=dqdt(:,:)
      call mkdca
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') (dtdt(jhor,jlev)-zt0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dqdt(jhor,jlev)-zq0(jhor,jlev), jlev=1,NLEV)
      enddo
      end program harness
"""


@FORTRAN
def test_the_python_dry_adjustment_is_plasims_mkdca(tmp_path):
    source = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'rainmod.f90'
    if not source.exists():
        pytest.skip('no ExoPlaSim source at hand')
    text = source.read_text()
    code = text[text.index('      subroutine mkdca\n'):text.index('      end subroutine mkdca\n')] + \
        '      end subroutine mkdca\n'
    rng = np.random.default_rng(11)
    for c, half, cols in cases()[:2]:
        sigma = ps.full_sigma(half)
        akap = c['gascon'] / c['acpd']
        dsig = np.diff(np.r_[0.0, half])
        n = half.size
        for x in cols:                               # warm the lowest layers so the columns adjust
            k = int(rng.integers(n - 4, n))
            x['dtdt'] = x['dtdt'] + np.where(np.arange(n) >= k, rng.uniform(0.0, 2.0e-3), 0.0)
        # and one smooth stable column, which the adjustment leaves alone
        z = -c['gascon'] * 280.0 / c['ga'] * np.log(sigma)
        stable = dict(cols[0], dt=np.r_[(300.0 + 3.0e-3 * z) * sigma ** akap, 302.0], dtdt=np.zeros(n),
                      dqdt=np.zeros(n))
        cols = cols + [stable]
        folder = tmp_path / f'dca{n}'
        folder.mkdir()
        files = {'stub.f90': DCA_STUB.format(nhor=len(cols), nlev=n), 'mkdca.f90': code, 'harness.f90': DCA_HARNESS}
        for name, body in files.items():
            (folder / name).write_text(body)
        subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
        lines = [num([akap, 3600.0, c['ga'], c['acpd'], c['ADV']]), num(np.r_[sigma, dsig])]
        for x in cols:
            lines += [num(x['dp'])] + [num(x[key][:n]) for key in ('dt', 'dq', 'dtdt', 'dqdt')]
        rows = subprocess.run([str(folder / 'harness')], input='\n'.join(lines) + '\n', capture_output=True, text=True,
                              check=True).stdout.strip('\n').split('\n')
        adjusted = 0
        for j, x in enumerate(cols):
            t = x['dt'][:n] + x['dtdt'] * 3600.0
            q = x['dq'][:n] + x['dqdt'] * 3600.0
            dtd, dqd, passes = ps.dca(t, q, sigma, half, akap, 3600.0)
            ft, fq = (np.array(rows[2 * j + i].split(), float) for i in (0, 1))
            assert np.allclose(ft, dtd, rtol=1e-10, atol=1e-14) and np.allclose(fq, dqd, rtol=1e-10, atol=1e-18), j
            adjusted += passes > 0
        assert adjusted == len(cols) - 1             # every warmed column adjusts, the stable one is left alone
        # PlaSim leaves the round-off of T / sigma^kappa * sigma^kappa in an untouched column's temperature
        assert np.abs(np.array(rows[-2].split(), float)).max() < 1.0e-15
        assert not np.any(np.array(rows[-1].split(), float))
