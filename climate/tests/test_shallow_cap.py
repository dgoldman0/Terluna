"""The shallow cap with the simplified Betts-Miller scheme (exoplasim_run.SBM_CAP_PATCHES, 2026-10-26): with the option
off the routine prints every digit it printed before; with it on, a shallow cloud's top stands no higher than the
depth test's threshold above its own cloud base, as the Python copy puts it, and PlaSim's absolute 700 hPa cap gives
way to it. The patch applies once to both installed sources after the depth patch, with or without the CAPE trigger."""
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
GCM_FLAGS = ('-O3', '-cpp', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffpe-summary=none',
             '-finit-real=zero', '-fdefault-real-8')
FORTRAN = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')
T0 = np.array([165.97, 194.74, 221.45, 243.81, 260.84, 272.61, 281.35, 289.48, 295.52, 299.16])
Q0 = np.array([1e-8, 2.7e-6, 6.02e-5, 4.808e-4, 1.3599e-3, 2.9967e-3, 4.7614e-3, 6.7683e-3, 1.13541e-2, 1.32844e-2])

STUB = """
      module rainmod
      integer, parameter :: NHOR=1, NLEV=10, NLEM=9, NLEP=11
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dqt(NHOR,NLEV), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV), dgp3d(NHOR,NLEV,20)
      integer :: icclev(NHOR,NLEV), icctot(NHOR)
      real :: akap, ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon, deltsec2, pdeep
      real :: tau_bm, rh_bm, eps_bm, cin_bm, cape_bm, epst_bm, dpdeep_bm
      integer :: nsub_bm = 1, nshcap_bm = 0
      integer :: nshallow = 0, ndiaggp = 0
      integer :: kshal(NHOR), ktopo(NHOR)
      end module rainmod

      subroutine mkshallow(pt,pq,kshallow,klift,ktop,pdtdt,pdqdt)
      use rainmod
      integer kshallow(NHOR), klift(NHOR), ktop(NHOR)
      real pt(NHOR,NLEV), pq(NHOR,NLEV), pdtdt(NHOR,NLEV), pdqdt(NHOR,NLEV)
      kshal(:)=kshallow(:)
      ktopo(:)=ktop(:)
      pdtdt(:,:)=0.
      pdqdt(:,:)=0.
      end subroutine mkshallow
"""
HARNESS = """
      program harness
      use rainmod
      integer :: k, jcase, ncase
      read(*,*) sigma, dsigma, dp, deltsec2, akap
      read(*,*) ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon
      read(*,*) ncase
      pdeep=999999.
      nshallow=1
      do jcase=1,ncase
       read(*,*) tau_bm, rh_bm, eps_bm, cin_bm, cape_bm, epst_bm, dpdeep_bm, nsub_bm, nshcap_bm
       read(*,*) (dt(1,k), k=1,NLEV)
       read(*,*) (dq(1,k), k=1,NLEV)
       dtdt(:,:)=0.
       dqdt(:,:)=0.
       dqt(:,:)=0.
       kshal(:)=0
       ktopo(:)=0
       call mksbm
       write(*,'(10es24.16)') (dtdt(1,k), k=1,NLEV)
       write(*,'(10es24.16)') (dqdt(1,k), k=1,NLEV)
       write(*,'(10es24.16)') (dprcl(1,k)*dp(1)/ga, k=1,NLEV)
       write(*,'(2i6)') kshal(1), ktopo(1)
      enddo
      end program harness
"""
ARGS = dict(tau=7200.0, rh=0.7, eps=0.0, cin_limit=0.0, cape_limit=0.0, eps_trigger=0.0, deep_pa=20000.0, steps=1,
            cap=0)
# Columns: the sea's mean column with its vapour scaled, warmed above or capped by an inversion, at the scheme's
# settings; a moist reference (rh 0.9-0.95) leaves a deep, buoyant column too dry to rain, which goes to the shallow
# scheme with its parcel's own top far above its base.
CASES = ((1.15, 0.0, 0.0, dict(rh=0.5)), (1.15, 0.0, 0.0, dict(rh=0.95)), (1.3, 0.0, 0.0, dict(rh=0.9)),
         (1.3, 0.0, 6.0, {}), (1.15, 0.0, 0.0, dict(rh=0.95, deep_pa=50000.0)), (1.15, 0.0, 0.0, dict(rh=0.95, eps=0.5)),
         (1.15, 25.0, 0.0, {}), (0.4, 0.0, 0.0, {}), (1.3, 0.0, 0.0, dict(rh=0.95, eps=2.4, steps=16)),
         (1.15, 0.0, 0.0, dict(rh=0.95, cape_limit=20.0)))


def column(moist, warm_top, inversion):
    t = T0 + warm_top * (SIGMA < 0.95) + inversion * (SIGMA < 0.88)
    return t, Q0 * moist


def routine(cap=True, trigger=True):
    """The scheme's routine with the trigger's, the depth patch's and the cap's edits on its text, as the runner
    applies them (the module's own edits go into the stub)."""
    code = er.SBM_CODE
    for edits in [er.SBM_TRIGGER_EDITS] * trigger + [er.SBM_DEPTH_EDITS] + [er.SBM_CAP_EDITS] * cap:
        for old, new in edits:
            if old in code:
                assert code.count(old) == 1
                code = code.replace(old, new)
    return code


def build(folder, code):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'stub.f90').write_text(STUB)
    (folder / 'mksbm.f90').write_text(code)
    (folder / 'harness.f90').write_text(HARNESS)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', 'stub.f90', 'mksbm.f90', 'harness.f90'], cwd=folder,
                   check=True, capture_output=True)
    return folder / 'harness'


def harness_run(exe, cases) -> str:
    c = sbm.PLASIM
    num = lambda xs: ' '.join(repr(float(x)) for x in xs)
    lines = [num(SIGMA) + ' ' + num(DSIGMA) + f' {PS!r} {STEP!r} {AKAP!r}',
             ' '.join(repr(c[k]) for k in ('ga', 'acpd', 'ADV', 'ALV', 'ALS', 'TMELT', 'rdbrv', 'ra1', 'ra2', 'ra4',
                                            'gascon')), str(len(cases))]
    for t, q, args in cases:
        a = ARGS | args
        lines += [' '.join(repr(float(a[k])) for k in ('tau', 'rh', 'eps', 'cin_limit', 'cape_limit', 'eps_trigger',
                                                        'deep_pa')) + f" {int(a['steps'])} {int(a['cap'])}",
                  num(t), num(q)]
    return subprocess.run([str(exe)], cwd=exe.parent, input='\n'.join(lines) + '\n', capture_output=True, text=True,
                          check=True).stdout


def parsed(printed, n):
    rows = printed.strip('\n').split('\n')
    assert len(rows) == 4 * n
    return [dict(dtdt=np.array(rows[4 * i].split(), float), dqdt=np.array(rows[4 * i + 1].split(), float),
                 shallow=bool(int(rows[4 * i + 3].split()[0])), top=int(rows[4 * i + 3].split()[1])) for i in range(n)]


def python(t, q, args):
    a = ARGS | args
    return sbm.column(t, q, SIGMA, DSIGMA, PS, tau=a['tau'], rh=a['rh'], eps=a['eps'], cin_limit=a['cin_limit'],
                      step=STEP, akap=AKAP, cape_limit=a['cape_limit'], eps_trigger=a['eps_trigger'],
                      deep_pa=a['deep_pa'], steps=a['steps'], shallow_cap=bool(a['cap']))


def test_the_python_copy_caps_a_deep_column_too_dry_to_rain_within_the_threshold_of_its_base():
    t, q = column(1.15, 0.0, 0.0)
    plain = python(t, q, dict(rh=0.95))
    capped = python(t, q, dict(rh=0.95, cap=1))
    p = PS * SIGMA
    assert plain['shallow'] and capped['shallow'] and plain['depth'] > ARGS['deep_pa']
    assert capped['parcel_top'] == plain['top'] < capped['top']
    assert p[capped['top']] >= plain['base'] - ARGS['deep_pa'] > p[capped['top'] - 1]
    thin = python(*column(1.3, 0.0, 6.0), dict(cap=1))                # sent by the depth test: its top already within
    assert thin['shallow'] and thin['top'] == python(*column(1.3, 0.0, 6.0), {})['top'] and 'parcel_top' not in thin
    deep = python(t, q, dict(rh=0.5, cap=1))                          # a raining column is not touched
    assert not deep['shallow'] and np.array_equal(deep['dtdt'], python(t, q, dict(rh=0.5))['dtdt'])


@FORTRAN
def test_with_the_cap_off_the_routine_is_the_scheme_as_it_was(tmp_path):
    cases = [(*column(m, w, i), kw) for m, w, i, kw in CASES]
    for trigger in (False, True):
        use = cases if trigger else [c for c in cases if not c[2].get('cape_limit')]
        before = harness_run(build(tmp_path / f'before{trigger}', routine(cap=False, trigger=trigger)), use)
        after = harness_run(build(tmp_path / f'after{trigger}', routine(cap=True, trigger=trigger)), use)
        assert before == after                                       # every printed digit, so every bit


@FORTRAN
def test_with_the_cap_on_the_fortran_gives_the_python_copys_tops(tmp_path):
    cases = [(*column(m, w, i), dict(kw, cap=1)) for m, w, i, kw in CASES]
    printed = parsed(harness_run(build(tmp_path / 'cap', routine()), cases), len(cases))
    capped = 0
    for (t, q, kw), fortran in zip(cases, printed):
        py = python(t, q, kw)
        assert fortran['shallow'] == py['shallow'], kw
        assert np.allclose(fortran['dtdt'], py['dtdt'], rtol=1e-10, atol=1e-18), kw
        if py['shallow']:
            assert fortran['top'] == py['top'] + 1, kw                # Fortran's levels count from 1
            capped += 'parcel_top' in py and py['top'] > py['parcel_top']
    assert capped >= 3


def installed_sources() -> list:
    found = [er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'rainmod.f90']
    tail = 'climate/gcm/.venv/lib/python*/site-packages/exoplasim/plasim/src/rainmod.f90'
    found += sorted(er.ROOT.glob(tail)) + sorted(er.ROOT.parent.glob('*/' + tail))
    return [p for p in dict.fromkeys(found) if p.exists()]


def apply(text, edits, marker):
    if marker in text:
        return text
    for old, new in edits:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    return text


def test_the_cap_patch_applies_once_to_both_sources_after_the_depth_patch():
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    marker = er.SBM_CAP_MARKER
    assert er.SBM_CAP_PATCHES == [('rainmod.f90', marker, er.SBM_CAP_EDITS)]
    assert not any(marker in old for old, _ in er.SBM_CAP_EDITS) and any(marker in new for _, new in er.SBM_CAP_EDITS)
    for _, new in er.SBM_CAP_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines())
    for path in sources:
        text = apply(path.read_text(), er.SBM_EDITS, er.SBM_MARKER)
        for trigger in (False, True):
            base = apply(text, er.SBM_TRIGGER_EDITS, er.SBM_TRIGGER_MARKER) if trigger else text
            if not trigger and er.SBM_TRIGGER_MARKER in text:
                continue
            patched = apply(apply(base, er.SBM_DEPTH_EDITS, er.SBM_DEPTH_MARKER), er.SBM_CAP_EDITS, marker)
            assert patched.count('namelist/rainmod_nl/nshcap_bm') == 1 and patched.count('call mpbci(nshcap_bm)') == 1
            assert patched.count('if(nshcap_bm == 1 .and. dpdeep_bm > 0.) then') == 1
            assert patched.count('ktop(:)=MAX(ktop(:),kktop(:))') == 1             # still once, now behind the switch
            assert patched.count('if(nprc /= 2 .or. nshcap_bm /= 1) ktop(:)=MAX(ktop(:),kktop(:))') == 1
            assert apply(patched, er.SBM_CAP_EDITS, marker) == patched


def test_the_configuration_records_the_cap_only_when_set_and_only_with_the_depth_test():
    base = er.convection_setting(deep_hpa=200.0)
    assert er.convection_setting(deep_hpa=200.0, shallow_cap=None) == base
    assert er.convection_setting(deep_hpa=200.0, shallow_cap=False) == base
    capped = er.convection_setting(deep_hpa=200.0, shallow_cap=True)
    assert capped == dict(base, shallow_cap=True)
    assert er.convection_namelist(capped)[-1] == ('NSHCAP_BM', '1')
    assert er.convection_namelist(base)[-1] == ('DPDEEP_BM', '20000.0')
    with pytest.raises(ValueError):
        er.convection_setting(shallow_cap=True)                      # it caps at the depth test's threshold
    assert er.uses_cap_patch(dict(deep_hpa=200.0, shallow_cap=True))
    assert not er.uses_cap_patch(None) and not er.uses_cap_patch(dict(deep_hpa=200.0, shallow_cap=None))
