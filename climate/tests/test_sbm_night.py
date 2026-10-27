"""The night switch of the simplified Betts-Miller scheme, a diagnostic and never a configuration
(exoplasim_run.SBM_NIGHT_PATCHES, 2026-10-27). With the switch off, or on where the Sun is up, the routine prints every
digit it printed before. With it on where the geometric cosine of the Sun's zenith is at most 0, every column takes the
branch a column held back by the CAPE trigger takes, PlaSim's shallow scheme with the shallow cap, so it prints every
digit the routine prints with a CAPE limit no column meets. The scheme reads the cosine through sbmnight, a routine of
radmod.f90 that hands on the twilight patch's gmuz. The patch applies once over the CAPE trigger's, the depth test's and
the cap's patches, the patched sources compile with PlaSim's modules, and the runner records the switch, labelled a
diagnostic, only for a run that sets it."""
import shutil
import signal
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import pbl_fortran as pbf
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
      integer :: nsub_bm = 1, nshcap_bm = 0, ndiagnight_bm = 0
      integer :: nshallow = 0, ndiaggp = 0
      integer :: kshal(NHOR), ktopo(NHOR)
      real :: zmuz_test(NHOR)
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

      subroutine sbmnight(pmuz)
      use rainmod
      real pmuz(NHOR)
      pmuz(:)=zmuz_test(:)
      end subroutine sbmnight
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
       read(*,*) tau_bm, rh_bm, eps_bm, cin_bm, cape_bm, epst_bm, dpdeep_bm, nsub_bm, nshcap_bm, &
     &           ndiagnight_bm, zmuz_test(1)
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
            cap=1, night=0, muz=0.5)
# The shallow cap's test columns (test_shallow_cap.py): the sea's mean column with its vapour scaled, warmed above or
# capped by an inversion, at the run's settings (the depth test at 200 hPa and the shallow cap) and others.
CASES = ((1.15, 0.0, 0.0, dict(rh=0.5)), (1.15, 0.0, 0.0, dict(rh=0.95)), (1.3, 0.0, 0.0, dict(rh=0.9)),
         (1.3, 0.0, 6.0, {}), (1.15, 0.0, 0.0, dict(rh=0.95, deep_pa=50000.0)),
         (1.15, 0.0, 0.0, dict(rh=0.95, eps=0.5)),
         (1.15, 25.0, 0.0, {}), (0.4, 0.0, 0.0, {}), (1.3, 0.0, 0.0, dict(rh=0.95, eps=2.4, steps=16)),
         (1.15, 0.0, 0.0, dict(rh=0.95, cape_limit=20.0)), (1.3, 0.0, 0.0, {}), (1.15, 0.0, 0.0, {}),
         (1.3, 0.0, 0.0, dict(cap=0)), (1.15, 0.0, 0.0, dict(rh=0.5, cap=0)), (1.3, 0.0, 0.0, dict(rh=0.5)),
         (1.15, 0.0, 0.0, dict(rh=0.6)), (1.3, 0.0, 0.0, dict(rh=0.6, deep_pa=10000.0)))


def column(moist, warm_top, inversion):
    t = T0 + warm_top * (SIGMA < 0.95) + inversion * (SIGMA < 0.88)
    return t, Q0 * moist


def routine(night=True):
    """The scheme's routine with the trigger's, the depth test's, the cap's and (night) the night switch's edits on its
    text, as the runner applies them (the module's own edits go into the stub)."""
    code = er.SBM_CODE
    for edits in [er.SBM_TRIGGER_EDITS, er.SBM_DEPTH_EDITS, er.SBM_CAP_EDITS] + [er.SBM_NIGHT_EDITS] * night:
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
                                                        'deep_pa'))
                  + f" {int(a['steps'])} {int(a['cap'])} {int(a['night'])} {float(a['muz'])!r}", num(t), num(q)]
    return subprocess.run([str(exe)], cwd=exe.parent, input='\n'.join(lines) + '\n', capture_output=True, text=True,
                          check=True).stdout


def parsed(printed, n):
    rows = printed.strip('\n').split('\n')
    assert len(rows) == 4 * n
    return [dict(dtdt=np.array(rows[4 * i].split(), float), shallow=bool(int(rows[4 * i + 3].split()[0])),
                 top=int(rows[4 * i + 3].split()[1])) for i in range(n)]


def cases(**over):
    return [(*column(m, w, i), dict(kw, **over)) for m, w, i, kw in CASES]


@FORTRAN
def test_with_the_switch_off_or_the_sun_up_the_routine_prints_every_digit_it_printed_before(tmp_path):
    before = harness_run(build(tmp_path / 'before', routine(night=False)), cases())
    after = build(tmp_path / 'after', routine())
    assert harness_run(after, cases(night=0, muz=-0.5)) == before            # off: the cosine is not read
    for muz in (1e-6, 0.3, 1.0):                                             # on, the Sun up
        assert harness_run(after, cases(night=1, muz=muz)) == before


@FORTRAN
def test_with_the_sun_down_every_column_takes_the_branch_of_a_column_the_cape_trigger_holds_back(tmp_path):
    held = harness_run(build(tmp_path / 'before', routine(night=False)), cases(cape_limit=1e30))
    after = build(tmp_path / 'after', routine())
    for muz in (0.0, -1e-6, -0.6):
        assert harness_run(after, cases(night=1, muz=muz)) == held          # every printed digit, so every bit
    day = parsed(harness_run(after, cases(night=1, muz=0.5)), len(CASES))
    night = parsed(held, len(CASES))
    deep = [i for i, d in enumerate(day) if not d['shallow'] and np.any(d['dtdt'])]
    assert len(deep) >= 3                                                    # the cases convect deep by day
    for i in deep:
        assert night[i]['shallow'] and not np.any(night[i]['dtdt'])         # and hand those columns on at night


def installed_sources(name='rainmod.f90') -> list:
    found = [er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / name]
    tail = f'climate/gcm/.venv/lib/python*/site-packages/exoplasim/plasim/src/{name}'
    found += sorted(er.ROOT.glob(tail)) + sorted(er.ROOT.parent.glob('*/' + tail))
    return [p for p in dict.fromkeys(found) if p.exists()]


def apply(text, edits, marker):
    if marker in text:
        return text
    for old, new in edits:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    return text


def patched_rainmod(text):
    """A rainmod.f90 source with the patches a night run's source carries, in the runner's order."""
    from climate.gcm import massflux_fortran as mff
    for edits, marker in ((er.SBM_EDITS, er.SBM_MARKER), (er.SBM_TRIGGER_EDITS, er.SBM_TRIGGER_MARKER),
                          (er.SBM_DEPTH_EDITS, er.SBM_DEPTH_MARKER), (er.SBM_CAP_EDITS, er.SBM_CAP_MARKER),
                          (mff.MASSFLUX_EDITS, mff.MASSFLUX_MARKER), (pbf.DIAGNOSTIC_EDITS, pbf.DIAGNOSTIC_MARKER),
                          (er.SBM_NIGHT_EDITS, er.SBM_NIGHT_MARKER)):
        text = apply(text, edits, marker)
    return text


def test_the_night_patch_applies_once_over_the_trigger_depth_and_cap_patches():
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    marker = er.SBM_NIGHT_MARKER
    assert er.SBM_NIGHT_PATCHES == [('rainmod.f90', marker, er.SBM_NIGHT_EDITS),
                                    ('radmod.f90', marker, er.SBM_NIGHT_RADMOD_EDITS)]
    for edits in (er.SBM_NIGHT_EDITS, er.SBM_NIGHT_RADMOD_EDITS):
        assert not any(marker in old for old, _ in edits) and any(marker in new for _, new in edits)
        assert all(len(line) <= 132 for _, new in edits for line in new.splitlines())
    assert 'diagnostic' in marker
    for path in sources:
        text = patched_rainmod(path.read_text())
        assert text.count(marker) == 4
        assert text.count('namelist/rainmod_nl/ndiagnight_bm') == 1 and text.count('call mpbci(ndiagnight_bm)') == 1
        assert text.count('if(ndiagnight_bm == 1) call sbmnight(zmuz)') == 1
        assert text.count('.or. (ndiagnight_bm == 1 .and. zmuz(jhor) <= 0.)) then') == 1
        assert apply(text, er.SBM_NIGHT_EDITS, marker) == text                     # once only
    rad = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'radmod.f90'
    if rad.exists():
        text = rad.read_text()
        assert 'gmuz(jhor) = zmuz' in text                                       # the twilight patch's cosine
        text = apply(text, er.SBM_NIGHT_RADMOD_EDITS, marker)
        assert text.count('subroutine sbmnight(pmuz)') == 1 and text.count('pmuz(:)=gmuz(:)') == 1
        assert text.index('end subroutine solang') < text.index('subroutine sbmnight(pmuz)')


@FORTRAN
def test_the_patched_sources_compile_with_plasims_modules_and_the_getter_hands_on_the_cosine(tmp_path):
    sources = installed_sources()
    rad = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'radmod.f90'
    if not sources or not rad.exists():
        pytest.skip('no ExoPlaSim source at hand')
    src = rad.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    for name in ('plasimmod.f90', 'specblock.f90'):
        shutil.copy(src / name, tmp_path / name)
    (tmp_path / 'rainmod.f90').write_text(patched_rainmod(sources[0].read_text()))
    (tmp_path / 'radmod.f90').write_text(apply(rad.read_text(), er.SBM_NIGHT_RADMOD_EDITS, er.SBM_NIGHT_MARKER))
    for name in ('resmod.f90', 'plasimmod.f90', 'specblock.f90', 'rainmod.f90', 'radmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)
    # the getter on a module holding the cosine, below the horizon too
    folder = tmp_path / 'getter'
    folder.mkdir()
    (folder / 'radmod.f90').write_text('      module radmod\n      integer, parameter :: NHOR=4\n'
                                       '      real :: gmuz(NHOR)\n      end module radmod\n'
                                       + er.SBM_NIGHT_GETTER)
    (folder / 'harness.f90').write_text('      program harness\n      use radmod\n      real :: z(NHOR)\n'
                                        '      gmuz(:)=(/0.7, 0.0, -0.25, -1.0/)\n      z(:)=9.\n'
                                        '      call sbmnight(z)\n      write(*,*) z\n      end program harness\n')
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', 'radmod.f90', 'harness.f90'], cwd=folder, check=True,
                   capture_output=True)
    out = subprocess.run([str(folder / 'harness')], capture_output=True, text=True, check=True).stdout
    assert np.allclose(np.array(out.split(), float), [0.7, 0.0, -0.25, -1.0])


def test_the_runner_records_the_switch_only_when_set_and_labels_it_a_diagnostic():
    base = er.convection_setting(deep_hpa=200.0, shallow_cap=True)
    assert er.convection_setting(deep_hpa=200.0, shallow_cap=True, diagnostic_night=None) == base
    assert er.convection_setting(deep_hpa=200.0, shallow_cap=True, diagnostic_night=False) == base
    night = er.convection_setting(deep_hpa=200.0, shallow_cap=True, diagnostic_night=True)
    assert night == dict(base, diagnostic_night=er.SBM_NIGHT_RECORD)
    assert 'diagnostic' in night['diagnostic_night']['status']
    assert er.convection_namelist(night)[-1] == ('NDIAGNIGHT_BM', '1')
    assert ('NDIAGNIGHT_BM', '1') not in er.convection_namelist(base)
    assert er.uses_night_switch(dict(deep_hpa=200.0, diagnostic_night=True))
    assert not er.uses_night_switch(None) and not er.uses_night_switch(dict(deep_hpa=200.0, diagnostic_night=None))
    with pytest.raises(ValueError):
        er.convection_setting(scheme='massflux', diagnostic_night=True)       # the Betts-Miller scheme's alone
    # refused without --convection sbm (main sets the runner's signal handlers before it checks, so the test puts the
    # process's own back)
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        with pytest.raises(SystemExit):
            er.main(['A28R', '--sbm-diagnostic-night'])
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
    text = (er.HERE / 'exoplasim_run.py').read_text()
    assert "'--sbm-diagnostic-night'" in text and "help='a diagnostic, never a configuration" in text
    # the runner's chain: the switch over the trigger, depth and cap patches, with the convection's vapour diagnostic,
    # and the vapour by process for a run with it and the moisture diagnostics
    assert '(SBM_TRIGGER_PATCHES if trigger or night else [])' in text
    assert '(SBM_CAP_PATCHES if cap or night else [])' in text
    assert '(pbf.DIAGNOSTIC_PATCHES if boundary_layer or night else [])' in text
    assert text.index('(pbf.DIAGNOSTIC_PATCHES if boundary_layer or night') < text.index('(SBM_NIGHT_PATCHES if night')
    assert 'boundary_layer is not None or uses_night_switch(convection):' in text
