"""The simplified Betts-Miller scheme patched into PlaSim: the patch applies once to ExoPlaSim's rainmod.f90 with Kuo's
path kept, the Fortran routine compiled alone against a stub of PlaSim's module gives the Python copy's tendencies
column by column, and both keep the column's moist enthalpy, relax toward their references and rain only from a
buoyant column moister than its reference."""
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import sbm

SIGMA = np.array([0.038, 0.119, 0.211, 0.317, 0.437, 0.567, 0.699, 0.823, 0.924, 0.983])
EDGES = [0.0]
for s in SIGMA:
    EDGES.append(2.0 * s - EDGES[-1])
DSIGMA = np.diff(EDGES)
PS, STEP, AKAP = 121490.0, 1800.0, 287.05 / 1004.6


def sounding(warm_top=0.0, moist=1.0):
    """The Near-side Sea's mean column (gcm_site, A28R_301k_cold_moist years 0-5), top first, every level above the
    lowest warmed by warm_top K and its vapour scaled by moist."""
    t = np.array([165.97, 194.74, 221.45, 243.81, 260.84, 272.61, 281.35, 289.48, 295.52, 299.16])
    q = np.array([1e-8, 2.7e-6, 6.02e-5, 4.808e-4, 1.3599e-3, 2.9967e-3, 4.7614e-3, 6.7683e-3, 1.13541e-2, 1.32844e-2])
    t = t + warm_top * (SIGMA < 0.95)
    return t, q * moist


def run(t, q, **kw):
    return sbm.column(t, q, SIGMA, DSIGMA, PS, step=STEP, akap=AKAP, **kw)


def enthalpy(out, q):
    cp = sbm.PLASIM['acpd'] * (1.0 + sbm.PLASIM['ADV'] * q)
    m = PS * DSIGMA / sbm.PLASIM['ga']
    return float(((cp * out['dtdt'] + sbm.PLASIM['ALV'] * out['dqdt']) * m).sum())


def test_a_moist_unstable_column_rains_its_vapour_and_keeps_its_enthalpy():
    t, q = sounding(moist=1.15)
    out = run(t, q, tau=7200.0, rh=0.5)
    assert out['top'] is not None and out['top'] <= 3 and not out['shallow']
    rain = out['rain'].sum()
    m = PS * DSIGMA / sbm.PLASIM['ga']
    assert rain > 0.0 and rain == pytest.approx(-(out['dqdt'] * m).sum(), rel=1e-12)
    assert abs(enthalpy(out, q)) < 1e-9 * sbm.PLASIM['ALV'] * rain
    # a longer relaxation time rains proportionally less from the same column
    slow = run(t, q, tau=4 * 7200.0, rh=0.5)
    assert slow['rain'].sum() == pytest.approx(rain / 4.0, rel=1e-9)


def test_a_dry_or_stable_column_does_not_convect_and_one_too_dry_to_rain_goes_shallow():
    for t, q in (sounding(moist=0.4), sounding(warm_top=10.0)):
        out = run(t, q, rh=0.5)
        assert out['top'] is None and out['rain'].sum() == 0.0 and np.all(out['dtdt'] == 0.0)
        assert not out['shallow']
    t, q = sounding()                                   # buoyant, but drier than 70% of its reference's saturation
    out = run(t, q, rh=0.7)
    assert out['top'] is not None and out['shallow'] and out['rain'].sum() == 0.0 and np.all(out['dqdt'] == 0.0)


def test_mixing_cools_the_reference_and_the_inhibition_limit_holds_convection_back():
    t, q = sounding(moist=1.15)
    plain, mixed = run(t, q, rh=0.5), run(t, q, rh=0.5, eps=0.5)
    # through the plain parcel's convecting layers the mixed one is cooler, so its top lies no higher
    layers = slice(plain['top'], -1)
    assert np.all(mixed['reference_t'][layers] < plain['reference_t'][layers])
    assert mixed['reference_t'][-1] == plain['reference_t'][-1]
    assert mixed['top'] is None or mixed['top'] >= plain['top']
    held = run(t, q, rh=0.5, cin_limit=1e-6)
    if plain['cin'] > 1e-6:
        assert held['rain'].sum() == 0.0 and held['shallow']


def test_the_patch_applies_once_and_keeps_kuos_path():
    exoplasim = pytest.importorskip('exoplasim')
    from pathlib import Path
    src = Path(exoplasim.__file__).parent / 'plasim' / 'src' / 'rainmod.f90'
    if not src.exists():
        pytest.skip(f'{src} is not at hand')
    text = src.read_text()
    if er.SBM_MARKER not in text:
        for old, new in er.SBM_EDITS:
            assert text.count(old) == 1, old
            text = text.replace(old, new)
    assert er.SBM_MARKER in text and text.count('subroutine mksbm') == 2
    assert '      if(nprc == 1) call kuo\n      if(nprc == 2) call mksbm' in text
    assert text.count('call kuo') == 1 and 'namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm' in text
    assert ('rainmod.f90', er.SBM_MARKER, er.SBM_EDITS) == er.PATCHES[-1]


STUB = """
      module rainmod
      integer, parameter :: NHOR=1, NLEV=10, NLEM=9, NLEP=11
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dqt(NHOR,NLEV), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV), dgp3d(NHOR,NLEV,20)
      integer :: icclev(NHOR,NLEV), icctot(NHOR)
      real :: akap, ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon, deltsec2, pdeep
      real :: tau_bm, rh_bm, eps_bm, cin_bm
      integer :: nshallow = 0, ndiaggp = 0
      end module rainmod

      subroutine mkshallow(pt,pq,kshallow,klift,ktop,pdtdt,pdqdt)
      use rainmod
      integer kshallow(NHOR), klift(NHOR), ktop(NHOR)
      real pt(NHOR,NLEV), pq(NHOR,NLEV), pdtdt(NHOR,NLEV), pdqdt(NHOR,NLEV)
      pdtdt(:,:)=0.
      pdqdt(:,:)=0.
      end subroutine mkshallow
"""
HARNESS = """
      program harness
      use rainmod
      integer :: k
      read(*,*) sigma, dsigma, dp, deltsec2, akap
      read(*,*) ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon
      read(*,*) tau_bm, rh_bm, eps_bm, cin_bm
      read(*,*) (dt(1,k), k=1,NLEV)
      read(*,*) (dq(1,k), k=1,NLEV)
      pdeep=999999.
      dtdt(:,:)=0.
      dqdt(:,:)=0.
      dqt(:,:)=0.
      call mksbm
      write(*,'(10es24.16)') (dtdt(1,k), k=1,NLEV)
      write(*,'(10es24.16)') (dqdt(1,k), k=1,NLEV)
      write(*,'(10es24.16)') (dprcl(1,k)*dp(1)/ga, k=1,NLEV)
      end program harness
"""


@pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')
def test_the_fortran_routine_gives_the_python_copys_tendencies(tmp_path):
    (tmp_path / 'stub.f90').write_text(STUB)
    (tmp_path / 'mksbm.f90').write_text(er.SBM_CODE)
    (tmp_path / 'harness.f90').write_text(HARNESS)
    subprocess.run(['gfortran', '-fdefault-real-8', '-O2', '-o', 'harness', 'stub.f90', 'mksbm.f90', 'harness.f90'],
                   cwd=tmp_path, check=True, capture_output=True)
    c = sbm.PLASIM
    for warm, moist, kw in ((0.0, 1.15, {}), (0.0, 1.15, dict(eps=0.5)), (0.0, 1.3, dict(tau=3600.0, rh=0.6)),
                            (0.0, 0.4, {}), (25.0, 1.0, {})):
        t, q = sounding(warm_top=warm, moist=moist)
        args = dict(tau=7200.0, rh=0.5, eps=0.0, cin_limit=0.0) | kw
        num = lambda xs: ' '.join(repr(float(x)) for x in xs)          # plain numbers, as Fortran reads them
        lines = [num(SIGMA) + ' ' + num(DSIGMA) + f' {PS!r} {STEP!r} {AKAP!r}',
                 ' '.join(repr(c[k]) for k in ('ga', 'acpd', 'ADV', 'ALV', 'ALS', 'TMELT', 'rdbrv', 'ra1', 'ra2',
                                                'ra4', 'gascon')),
                 f"{args['tau']!r} {args['rh']!r} {args['eps']!r} {args['cin_limit']!r}",
                 num(t), num(q)]
        done = subprocess.run(['./harness'], cwd=tmp_path, input='\n'.join(lines) + '\n', capture_output=True,
                              text=True, check=True)
        dtdt, dqdt, rain = (np.array(row.split(), float) for row in done.stdout.split('\n')[:3])
        out = run(t, q, **args)
        assert np.allclose(dtdt, out['dtdt'], rtol=1e-10, atol=1e-18), (warm, moist, kw)
        assert np.allclose(dqdt, out['dqdt'], rtol=1e-10, atol=1e-20), (warm, moist, kw)
        assert np.allclose(rain, out['rain'], rtol=1e-10, atol=1e-18), (warm, moist, kw)


def test_the_configuration_records_the_scheme_and_refuses_what_it_cannot_run():
    assert er.convection_setting() == dict(scheme='sbm', nprc=2, tau_s=7200.0, rh=0.7, eps=0.0, cin_j_kg=0.0)
    lunar = er.convection_setting(tau_s=21600, rh=0.5, eps=0.3, cin_j_kg=None)
    assert lunar['tau_s'] == 21600.0 and lunar['rh'] == 0.5 and lunar['eps'] == 0.3 and lunar['cin_j_kg'] == 0.0
    for bad in (dict(scheme='kuo'), dict(rh=1.2), dict(tau_s=0.0), dict(eps=-0.1), dict(width=3.0)):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
