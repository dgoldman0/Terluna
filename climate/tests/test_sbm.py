"""The simplified Betts-Miller scheme patched into PlaSim: the patch applies once to ExoPlaSim's rainmod.f90 with Kuo's
path kept, the Fortran routine compiled alone against a stub of PlaSim's module gives the Python copy's tendencies
column by column, and both keep the column's moist enthalpy, relax toward their references and rain only from a
buoyant column moister than its reference; the CAPE trigger, off by default, leaves the routine as it was and
holds back a column whose mixed parcel has too little CAPE. The cloud-depth patch's two options, each off by default,
leave the routine as it was; with them on, the Fortran gives the Python copy's tendencies and shallow decisions, the
depth test sends a cloud thinner than its threshold to the shallow scheme, and the parcel's mixing within layers
converges as its steps are refined. The patch applies to both installed sources, with or without the trigger."""
import math
import shutil
import subprocess
from pathlib import Path

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
# ExoPlaSim's own compiler options for PlaSim (most_compiler's MOST_F90_OPTS and MOST_PREC), less the fixed-form width
GCM_FLAGS = ('-O3', '-cpp', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffpe-summary=none',
             '-finit-real=zero', '-fdefault-real-8')
FORTRAN = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')


def sounding(warm_top=0.0, moist=1.0, inversion=0.0, warm_ground=0.0):
    """The Near-side Sea's mean column (gcm_site, A28R_301k_cold_moist years 0-5), top first, every level above the
    lowest warmed by warm_top K and its vapour scaled by moist; inversion warms the levels above sigma 0.88, which
    caps the parcel just above the lowest layer, and warm_ground the lowest level alone."""
    t = np.array([165.97, 194.74, 221.45, 243.81, 260.84, 272.61, 281.35, 289.48, 295.52, 299.16])
    q = np.array([1e-8, 2.7e-6, 6.02e-5, 4.808e-4, 1.3599e-3, 2.9967e-3, 4.7614e-3, 6.7683e-3, 1.13541e-2, 1.32844e-2])
    t = t + warm_top * (SIGMA < 0.95) + inversion * (SIGMA < 0.88) + warm_ground * (SIGMA > 0.95)
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
    if er.SBM_TRIGGER_MARKER not in text:                       # the CAPE trigger, a patch on the scheme's text
        for old, new in er.SBM_TRIGGER_EDITS:
            assert text.count(old) == 1, old
            text = text.replace(old, new)
    assert er.SBM_TRIGGER_MARKER in text and 'namelist/rainmod_nl/cape_bm,epst_bm' in text
    assert text.count('call mpbcr(cape_bm)') == 1 and text.count('zcapet < cape_bm') == 1


# The routine compiled alone: a stub of PlaSim's module, whose shallow scheme records the columns it is handed, and a
# harness that runs the routine on each of a list of columns and prints its tendencies, rain and shallow decision.
STUB = """
      module rainmod
      integer, parameter :: NHOR=1, NLEV=10, NLEM=9, NLEP=11
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dqt(NHOR,NLEV), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV), dgp3d(NHOR,NLEV,20)
      integer :: icclev(NHOR,NLEV), icctot(NHOR)
      real :: akap, ga, acpd, ADV, ALV, ALS, TMELT, rdbrv, ra1, ra2, ra4, gascon, deltsec2, pdeep
      real :: tau_bm, rh_bm, eps_bm, cin_bm, cape_bm, epst_bm, dpdeep_bm
      integer :: nsub_bm = 1
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
       read(*,*) tau_bm, rh_bm, eps_bm, cin_bm, cape_bm, epst_bm, dpdeep_bm, nsub_bm
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
ARGS = dict(tau=7200.0, rh=0.5, eps=0.0, cin_limit=0.0, cape_limit=0.0, eps_trigger=0.0, deep_pa=0.0, steps=1)


def routine(trigger=False, depth=False, depth_first=False):
    """The scheme's routine with the CAPE trigger's and the cloud-depth patch's edits on it, in the order the runner
    would apply them (the module's own edits go into the stub)."""
    code = er.SBM_CODE
    order = [(er.SBM_TRIGGER_EDITS, 2)] * trigger + [(er.SBM_DEPTH_EDITS, 4)] * depth
    for edits, expected in (order[::-1] if depth_first else order):
        applied = 0
        for old, new in edits:
            if old in code:
                assert code.count(old) == 1
                code, applied = code.replace(old, new), applied + 1
        assert applied == expected
    return code


def trigger_code():
    return routine(trigger=True)


def build(folder, code, flags=GCM_FLAGS):
    """The harness compiled with the routine's code: the executable's path."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'stub.f90').write_text(STUB)
    (folder / 'mksbm.f90').write_text(code)
    (folder / 'harness.f90').write_text(HARNESS)
    subprocess.run(['gfortran', *flags, '-o', 'harness', 'stub.f90', 'mksbm.f90', 'harness.f90'], cwd=folder,
                   check=True, capture_output=True)
    return folder / 'harness'


def harness_input(cases):
    """Each case is (t, q, args) with args as ARGS; plain numbers, as Fortran reads them."""
    c = sbm.PLASIM
    num = lambda xs: ' '.join(repr(float(x)) for x in xs)
    lines = [num(SIGMA) + ' ' + num(DSIGMA) + f' {PS!r} {STEP!r} {AKAP!r}',
             ' '.join(repr(c[k]) for k in ('ga', 'acpd', 'ADV', 'ALV', 'ALS', 'TMELT', 'rdbrv', 'ra1', 'ra2', 'ra4',
                                            'gascon')),
             str(len(cases))]
    for t, q, args in cases:
        a = ARGS | args
        lines += [' '.join(repr(float(a[k])) for k in ('tau', 'rh', 'eps', 'cin_limit', 'cape_limit', 'eps_trigger',
                                                        'deep_pa')) + f" {int(a['steps'])}", num(t), num(q)]
    return '\n'.join(lines) + '\n'


def harness_run(exe, cases) -> str:
    return subprocess.run([str(exe)], cwd=exe.parent, input=harness_input(cases), capture_output=True, text=True,
                          check=True).stdout


def parsed(printed, cases):
    """The harness's printed rows by case: tendencies, rain, and the shallow flag and top its shallow scheme saw."""
    rows = printed.strip('\n').split('\n')
    assert len(rows) == 4 * len(cases)
    out = []
    for n in range(len(cases)):
        dtdt, dqdt, rain = (np.array(rows[4 * n + i].split(), float) for i in range(3))
        flag, top = (int(v) for v in rows[4 * n + 3].split())
        out.append(dict(dtdt=dtdt, dqdt=dqdt, rain=rain, shallow=bool(flag), top=top))
    return out


def assert_matches(fortran, python, label):
    assert np.allclose(fortran['dtdt'], python['dtdt'], rtol=1e-10, atol=1e-18), label
    assert np.allclose(fortran['dqdt'], python['dqdt'], rtol=1e-10, atol=1e-20), label
    assert np.allclose(fortran['rain'], python['rain'], rtol=1e-10, atol=1e-18), label
    assert fortran['shallow'] == python['shallow'], label
    if python['shallow']:
        assert fortran['top'] == python['top'] + 1, label             # Fortran's levels count from 1


CASES = ((0.0, 1.15, {}), (0.0, 1.15, dict(eps=0.5)), (0.0, 1.3, dict(tau=3600.0, rh=0.6)), (0.0, 0.4, {}),
         (25.0, 1.0, {}), (0.0, 1.15, dict(cape_limit=50.0)), (0.0, 1.15, dict(cape_limit=1.0e4)),
         (0.0, 1.3, dict(cape_limit=5.0, eps_trigger=1.0, rh=0.6)), (0.0, 1.15, dict(cape_limit=20.0, eps=0.5,
                                                                                     eps_trigger=2.0)))


@FORTRAN
def test_the_fortran_routine_gives_the_python_copys_tendencies(tmp_path):
    cases = [(*sounding(warm_top=warm, moist=moist), kw) for warm, moist, kw in CASES]
    printed = parsed(harness_run(build(tmp_path / 'trigger', trigger_code()), cases), cases)
    for (t, q, kw), fortran in zip(cases, printed):
        assert_matches(fortran, run(t, q, **(ARGS | kw)), kw)


@FORTRAN
def test_with_the_trigger_off_the_fortran_routine_is_the_scheme_as_it_was(tmp_path):
    cases = [(*sounding(warm_top=warm, moist=moist), kw) for warm, moist, kw in CASES[:5]]
    assert (harness_run(build(tmp_path / 'trigger', trigger_code()), cases)
            == harness_run(build(tmp_path / 'plain', er.SBM_CODE), cases))


# Columns for the cloud-depth patch: (sounding's settings, the scheme's settings); a capped moist column convects
# only through the lowest layers, a saturated one starts its cloud at the lowest level, a dry one with a warm lowest
# level is buoyant without condensing below its top.
DEPTH_CASES = (
    (dict(moist=1.15), dict(deep_pa=20000.0)),                              # deep, its cloud 896 hPa thick
    (dict(moist=1.3, inversion=6.0), dict(deep_pa=20000.0)),                # capped: 65 hPa, to the shallow scheme
    (dict(moist=1.3, inversion=6.0), dict(deep_pa=5000.0)),                 # the same cloud above a 50 hPa test
    (dict(moist=1.15), dict(deep_pa=90000.0)),                              # deep but thinner than 900 hPa
    (dict(moist=1.4), dict(deep_pa=20000.0)),                               # saturated at the lowest level
    (dict(moist=0.3, warm_ground=8.0), dict(deep_pa=20000.0)),              # unsaturated up to its top
    (dict(moist=0.4), dict(deep_pa=20000.0)),                               # buoyant nowhere
    (dict(moist=1.15), dict(eps=2.4, steps=16)),                            # the mixing within layers alone
    (dict(moist=1.3), dict(eps=4.8, steps=4)),
    (dict(moist=1.15), dict(eps=0.5, steps=2, rh=0.6, tau=3600.0)),
    (dict(moist=1.15), dict(eps=2.4, steps=16, deep_pa=20000.0)),           # both
    (dict(moist=1.3, inversion=2.0), dict(eps=2.4, steps=16, deep_pa=20000.0)),
    (dict(moist=1.3, inversion=6.0), dict(eps=4.8, steps=16, deep_pa=20000.0)),
    (dict(moist=1.15), dict(eps=0.5, deep_pa=20000.0)),                     # the depth test with one-step mixing
    (dict(moist=1.15), dict(eps=4.8, steps=2, deep_pa=20000.0, cape_limit=20.0)),   # and with the CAPE trigger
    (dict(moist=1.3), dict(eps=2.4, steps=8, cape_limit=5.0, eps_trigger=1.0, rh=0.6)),
)


@FORTRAN
def test_with_both_new_options_off_the_routine_is_the_scheme_as_it_was(tmp_path):
    # every earlier setting, the CAPE trigger's included, with the depth test at 0 and one mixing step a layer
    cases = [(*sounding(warm_top=warm, moist=moist), kw) for warm, moist, kw in CASES]
    cases += [(*sounding(**col), {k: v for k, v in kw.items() if k not in ('deep_pa', 'steps')})
              for col, kw in DEPTH_CASES]
    plain_cases = [case for case in cases if not case[2].get('cape_limit')]
    for trigger in (False, True):
        before = harness_run(build(tmp_path / f'before{trigger}', routine(trigger=trigger)),
                             cases if trigger else plain_cases)
        after = harness_run(build(tmp_path / f'after{trigger}', routine(trigger=trigger, depth=True)),
                            cases if trigger else plain_cases)
        assert before == after                                      # every printed digit, so every bit
    # the patch on a source without the trigger, then the trigger on top, gives the same routine's arithmetic
    lines = harness_run(build(tmp_path / 'depth_first', routine(trigger=True, depth=True, depth_first=True)), cases)
    assert lines == before


@FORTRAN
def test_the_fortran_routine_gives_the_python_copys_tendencies_with_each_option_on(tmp_path):
    cases = [(*sounding(**col), kw) for col, kw in DEPTH_CASES]
    for name, code in (('both', routine(trigger=True, depth=True)),
                       ('depth_first', routine(trigger=True, depth=True, depth_first=True))):
        printed = parsed(harness_run(build(tmp_path / name, code), cases), cases)
        for (t, q, kw), fortran in zip(cases, printed):
            assert_matches(fortran, run(t, q, **(ARGS | kw)), (name, kw))
    plain = routine(depth=True)                                     # without the trigger in the source
    cases = [case for case in cases if not case[2].get('cape_limit')]
    printed = parsed(harness_run(build(tmp_path / 'no_trigger', plain), cases), cases)
    for (t, q, kw), fortran in zip(cases, printed):
        assert_matches(fortran, run(t, q, **(ARGS | kw)), ('no_trigger', kw))
    # the cases do what they were chosen for
    shallow = [run(t, q, **(ARGS | kw))['shallow'] for t, q, kw in cases[:7]]
    assert shallow == [False, True, False, True, False, True, False]


def test_the_cape_trigger_holds_back_a_column_whose_mixed_parcel_has_too_little_cape():
    t, q = sounding(moist=1.15)
    plain = run(t, q, rh=0.5)
    held = run(t, q, rh=0.5, cape_limit=1.0e4)
    assert held['trigger_cape'] < 1.0e4 and held['shallow'] and held['rain'].sum() == 0.0
    passed = run(t, q, rh=0.5, cape_limit=1.0)
    assert passed['trigger_cape'] > 1.0 and not passed['shallow']
    assert np.array_equal(passed['dtdt'], plain['dtdt']) and np.array_equal(passed['rain'], plain['rain'])
    capes = [run(t, q, rh=0.5, cape_limit=1.0, eps_trigger=e)['trigger_cape'] for e in (0.0, 0.5, 2.0)]
    assert capes[0] > capes[1] > capes[2] >= 0.0                        # mixing takes the trigger parcel's CAPE away
    assert 'trigger_cape' not in plain


def test_the_configuration_records_the_scheme_and_refuses_what_it_cannot_run():
    assert er.convection_setting() == dict(scheme='sbm', nprc=2, tau_s=7200.0, rh=0.7, eps=0.0, cin_j_kg=0.0)
    lunar = er.convection_setting(tau_s=21600, rh=0.5, eps=0.3, cin_j_kg=None)
    assert lunar['tau_s'] == 21600.0 and lunar['rh'] == 0.5 and lunar['eps'] == 0.3 and lunar['cin_j_kg'] == 0.0
    for bad in (dict(scheme='kuo'), dict(rh=1.2), dict(tau_s=0.0), dict(eps=-0.1), dict(width=3.0)):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    # the CAPE trigger joins the record only when set
    assert er.convection_setting(cape_j_kg=None, trigger_eps=None) == er.convection_setting()
    assert er.convection_setting(cape_j_kg=0.0) == er.convection_setting()
    trig = er.convection_setting(cape_j_kg=100.0)
    assert trig['cape_j_kg'] == 100.0 and trig['trigger_eps'] == 0.0
    assert er.convection_setting(cape_j_kg=20.0, trigger_eps=2.0)['trigger_eps'] == 2.0
    for bad in (dict(cape_j_kg=-1.0), dict(trigger_eps=1.0), dict(cape_j_kg=10.0, trigger_eps=-1.0)):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)


# ---------------------------------------------------------------------------------- the cloud-depth patch's options

def test_the_depth_test_sends_a_thin_cloud_to_the_shallow_scheme_and_leaves_a_deep_one():
    t, q = sounding(moist=1.3, inversion=6.0)
    capped = run(t, q, rh=0.5)
    assert capped['top'] == 8 and not capped['shallow'] and capped['rain'].sum() > 0.0     # topped at sigma 0.924
    held = run(t, q, rh=0.5, deep_pa=20000.0)
    assert held['shallow'] and held['rain'].sum() == 0.0 and np.all(held['dtdt'] == 0.0)
    assert held['base'] == pytest.approx(118720.0, rel=1e-3) and held['depth'] < 20000.0
    assert held['depth'] == pytest.approx(held['base'] - PS * SIGMA[8], rel=1e-12)
    assert np.array_equal(run(t, q, rh=0.5, deep_pa=5000.0)['rain'], capped['rain'])
    t, q = sounding(moist=1.15)
    deep, tested = run(t, q, rh=0.5), run(t, q, rh=0.5, deep_pa=20000.0)
    assert tested['depth'] > 80000.0 and not tested['shallow']
    for k in ('dtdt', 'dqdt', 'rain'):
        assert np.array_equal(deep[k], tested[k])
    assert run(t, q, rh=0.5, deep_pa=tested['depth'])['shallow']    # a cloud exactly the threshold's depth is shallow


def test_the_condensation_level_is_where_the_lifted_air_saturates():
    c, p0 = sbm.PLASIM, PS * SIGMA[-1]
    for moist in (0.9, 1.0, 1.15, 1.3):
        t, q = sounding(moist=moist)
        base = sbm.condensation_level(t[-1], q[-1], p0, PS * SIGMA[3], AKAP)
        assert PS * SIGMA[3] < base < p0
        # 24 halvings of a bracket at most 3.3 wide in log pressure: within 2e-7 of it, 1e-6 of the vapour
        assert sbm.saturation(t[-1] * (base / p0) ** AKAP, base, c) == pytest.approx(q[-1], rel=1e-6)
    t, q = sounding(moist=1.4)                                      # saturated where it starts
    assert sbm.condensation_level(t[-1], q[-1], p0, PS * SIGMA[3], AKAP) == p0
    t, q = sounding(moist=0.3)                                      # still unsaturated at the top it is given
    assert sbm.condensation_level(t[-1], q[-1], p0, PS * SIGMA[8], AKAP) == PS * SIGMA[8]


def old_parcel(t, q, sigma, p, eps, akap, c=sbm.PLASIM):
    """The scheme's parcel as it was before the mixing within layers (sbm.py at 467a146)."""
    tp = np.zeros(t.size)
    tn, qn = sbm.condense(t[-1], q[-1], p[-1], c)
    tp[-1] = tn
    for k in range(t.size - 2, -1, -1):
        tn *= (sigma[k] / sigma[k + 1]) ** akap
        if eps > 0.0:
            mix = 1.0 - math.exp(-eps * math.log(sigma[k + 1] / sigma[k]))
            tn, qn = (1.0 - mix) * tn + mix * t[k], (1.0 - mix) * qn + mix * q[k]
        tn, qn = sbm.condense(tn, qn, p[k], c)
        tp[k] = tn
    return tp


def test_the_mixing_within_layers_converges_and_one_step_is_the_scheme_as_before():
    p = PS * SIGMA
    for moist in (1.0, 1.15):
        t, q = sounding(moist=moist)
        for eps in (0.0, 0.5, 2.4):
            assert np.array_equal(sbm.parcel(t, q, SIGMA, p, eps, AKAP, steps=1, ps=PS),
                                  old_parcel(t, q, SIGMA, p, eps, AKAP))
        assert np.array_equal(sbm.parcel(t, q, SIGMA, p, 0.0, AKAP, steps=16, ps=PS), old_parcel(t, q, SIGMA, p, 0.0,
                                                                                                 AKAP))
        for eps in (2.4, 4.8):
            limit = sbm.parcel(t, q, SIGMA, p, eps, AKAP, steps=1024, ps=PS)
            errors = [np.abs(sbm.parcel(t, q, SIGMA, p, eps, AKAP, steps=n, ps=PS) - limit)[2:]
                      for n in (4, 8, 16, 32, 64)]
            worst = [float(e.max()) for e in errors]
            # first order: each doubling of the steps halves the error below sigma 0.2, the convecting troposphere
            assert all(1.7 < a / b < 2.3 for a, b in zip(worst, worst[1:])), worst
            assert worst[2] < 0.12                                  # 16 steps a layer: within about 0.1 K


def installed_sources() -> list:
    """The installed ExoPlaSim rainmod.f90 sources: the twilight copy every run takes, and ExoPlaSim's own (the GCM
    venv's, here or in a sibling checkout of the repository)."""
    found = [er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / 'rainmod.f90']
    tail = 'climate/gcm/.venv/lib/python*/site-packages/exoplasim/plasim/src/rainmod.f90'
    found += sorted(er.ROOT.glob(tail)) + sorted(er.ROOT.parent.glob('*/' + tail))
    return [p for p in dict.fromkeys(found) if p.exists()]


def apply(text, edits, marker):
    """A patch applied as ensure_patched applies it: each old text found exactly once, unless the marker is there."""
    if marker in text:
        return text
    for old, new in edits:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    return text


def test_the_depth_patch_applies_to_both_sources_with_or_without_the_trigger():
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    marker = er.SBM_DEPTH_MARKER
    assert any(marker in new for _, new in er.SBM_DEPTH_EDITS)
    assert not any(marker in old for old, _ in er.SBM_DEPTH_EDITS)
    assert er.SBM_DEPTH_PATCHES == [('rainmod.f90', marker, er.SBM_DEPTH_EDITS)]
    for _, new in er.SBM_DEPTH_EDITS:
        assert all(len(line) <= 132 for line in new.splitlines())
    for path in sources:
        text = apply(path.read_text(), er.SBM_EDITS, er.SBM_MARKER)
        orders = {'trigger, then depth': apply(apply(text, er.SBM_TRIGGER_EDITS, er.SBM_TRIGGER_MARKER),
                                               er.SBM_DEPTH_EDITS, marker)}
        if er.SBM_TRIGGER_MARKER not in text:
            alone = apply(text, er.SBM_DEPTH_EDITS, marker)
            orders['depth alone'] = alone
            orders['depth, then trigger'] = apply(alone, er.SBM_TRIGGER_EDITS, er.SBM_TRIGGER_MARKER)
        for label, patched in orders.items():
            assert patched.count(marker) == 2 and patched.count('function sbmlcl') == 2, (path, label)
            assert patched.count('namelist/rainmod_nl/dpdeep_bm,nsub_bm') == 1, (path, label)
            assert patched.count('call mpbcr(dpdeep_bm)') == 1 and patched.count('call mpbci(nsub_bm)') == 1
            assert patched.count('do jsub=1,nsub_bm') == 1 and patched.count('if(dpdeep_bm > 0.) then') == 1
            assert apply(patched, er.SBM_DEPTH_EDITS, marker) == patched           # once only
            assert er.SBM_TRIGGER_MARKER in patched or label == 'depth alone'
        if 'depth, then trigger' in orders:
            # either order yields the same lines, the trigger's and the depth patch's declarations interleaved
            assert sorted(orders['depth, then trigger'].split('\n')) == sorted(orders['trigger, then depth'].split('\n'))


def test_the_configuration_records_the_depth_patchs_options_only_when_set():
    base = er.convection_setting()
    assert er.convection_setting(deep_hpa=None, mixing_steps=None) == base
    assert er.convection_setting(deep_hpa=0.0, mixing_steps=1) == base
    assert er.convection_setting(deep_hpa=200)['deep_hpa'] == 200.0
    mixed = er.convection_setting(eps=2.4, mixing_steps=16, deep_hpa=200.0)
    assert mixed == dict(base, eps=2.4, deep_hpa=200.0, mixing_steps=16)
    for bad in (dict(deep_hpa=-1.0), dict(deep_hpa=1500.0), dict(mixing_steps=16), dict(eps=2.4, mixing_steps=0),
                dict(eps=2.4, mixing_steps=2.5)):
        with pytest.raises(ValueError):
            er.convection_setting(**bad)
    # the namelist: as before for the earlier settings, the depth in Pa and the steps as a whole number
    assert er.convection_namelist(base) == [('NPRC', '2'), ('TAU_BM', '7200.0'), ('RH_BM', '0.7'), ('EPS_BM', '0.0'),
                                            ('CIN_BM', '0.0')]
    assert er.convection_namelist(er.convection_setting(cape_j_kg=100.0))[-2:] == [('CAPE_BM', '100.0'),
                                                                                   ('EPST_BM', '0.0')]
    assert er.convection_namelist(mixed)[-2:] == [('DPDEEP_BM', '20000.0'), ('NSUB_BM', '16')]
    assert not er.uses_depth_patch(None) and not er.uses_depth_patch(dict(rh=0.7, deep_hpa=None, mixing_steps=None))
    assert er.uses_depth_patch(dict(deep_hpa=200.0)) and er.uses_depth_patch(dict(eps=2.4, mixing_steps=4))
