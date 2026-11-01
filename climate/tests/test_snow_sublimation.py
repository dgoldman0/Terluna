"""PlaSim's snow from its flakes (climate/gcm/snow_sublimation.py, nsnowsub 1): the scheme's closure gives the validated
closure_columns' flux under every matching and both table sets; the column keeps its water, sublimates only air
unsaturated over ice and only below the layers that form snow, splits graupel off by kind and melts all in a warm layer;
a column colder than 0 C to the ground takes its freezing level below the ground by the ladder's lapse; mkrain with
nsnowsub 0 is mkrain with gamma; compiled with ExoPlaSim's flags against a stub of PlaSim's modules, the patched mkrain
gives the Python copy's tendencies, precipitation and snow diagnostics to 1e-10, shadows on and off and with either top,
and stays finite in dry air and heavy snow; the patched rainmod.f90 and plasimmod.f90 compile with PlaSim's modules; and
the output patches apply once on lines no other patch touches (outmod.f90 compiles in the model's own build)."""
import math
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import levels as lv
from climate.gcm import rain_drops as rd
from climate.gcm import snow_sublimation as ss
from climate.tests.test_condensate import FORTRAN, GCM_FLAGS, apply, num, source, undo
from climate.tests.test_rain_drops import C, copy, mkrain_code

SNOW = dict(ss.PARAMETERS)


def snow_columns(n=8, seed=21, cold_ground=False, warm_layer=False, dry=False, plasim=True):
    """Lunar columns (top first) on the fifteen levels, freezing near 24 km: snow forming in one or two layers of the
    cold upper air at 0.01-3 mm/h, clear and 70-99% humid over ice in the layers between them and the freezing level,
    rain forming below in some; with cold_ground half the columns' ground and lowest layers below 0 C, and with
    warm_layer a warm layer standing over a cold one near the ground in some. With plasim, as PlaSim's condensation
    and convection make it, all of it forms as rain (dprscl and dprsll nought) and mkrain freezes it; without, the
    snow is formed as snow."""
    rng = np.random.default_rng(seed)
    half = np.asarray(lv.half_levels('l15'), float)
    nlev = half.size
    edges = np.r_[0.0, half]
    sigma = 0.5 * (edges[1:] + edges[:-1])
    dsigma = np.diff(edges)
    ps = 121500.0 + rng.uniform(-3000.0, 1000.0, n)
    z = -C['gascon'] * 285.0 / C['ga'] * np.log(sigma)
    t = (300.5 - 1.15e-3 * z)[None, :] + rng.uniform(-1.0, 1.0, (n, nlev))
    t = np.maximum(t, 190.0)
    p = ps[:, None] * sigma[None, :]
    ground = t[:, -1] + 1.0
    if cold_ground:
        cold = rng.uniform(size=n) < 0.5
        t[cold, -3:] = np.array([268.0, 266.0, 265.0])[None, :] + rng.uniform(-1.0, 1.0, (cold.sum(), 3))
        ground = np.where(cold, 263.0, ground)
        if warm_layer:
            some = cold & (rng.uniform(size=n) < 0.6)
            t[some, -3] = 276.0
    esi = C['ra1'] * np.exp(ss.MURRAY_ICE_A * (t - C['TMELT']) / (t - ss.MURRAY_ICE_B))
    esw = C['ra1'] * np.exp(C['ra2'] * (t - C['TMELT']) / (t - C['ra4']))
    e_sat = np.where(t <= C['TMELT'], esi, esw)
    rh = rng.uniform(0.05, 0.3, (n, nlev)) if dry else rng.uniform(0.7, 0.99, (n, nlev))
    q = C['rdbrv'] * rh * e_sat / p
    cold_air = t <= C['TMELT']
    conv, large = np.zeros((n, nlev)), np.zeros((n, nlev))
    for j in range(n):
        upper = [k for k in range(nlev) if cold_air[j, k] and z[k] > 30000.0 and z[k] < 80000.0]
        for k in rng.choice(upper, size=min(len(upper), int(rng.integers(1, 3))), replace=False):
            rate = math.exp(rng.uniform(math.log(0.01), math.log(3.0))) / 3600.0       # kg/m2/s
            each = rate / (ps[j] / C['ga'])
            share = rng.uniform(0.0, 1.0)
            conv[j, k], large[j, k] = share * each, (1.0 - share) * each
            q[j, k] = C['rdbrv'] * 1.01 * esw[j, k] / p[j, k]                      # the forming layer saturated
        if rng.uniform() < 0.5:
            warmk = [k for k in range(nlev) if not cold_air[j, k]]
            if warmk:
                k = int(rng.choice(warmk))
                conv[j, k] += rng.uniform(0.0, 4e-10)
    if n > 2:
        conv[1], large[1] = 0.0, 0.0                                                # one column without snow
    snow_c, snow_l = (np.zeros_like(conv), np.zeros_like(large)) if plasim else \
        (np.where(cold_air, conv, 0.0), np.where(cold_air, large, 0.0))
    return dict(dt=np.c_[t, ground], dq=q, dp=ps, sigma=sigma, dsigma=dsigma,
                dtdt=rng.uniform(-2e-5, 1e-5, (n, nlev)), dqdt=rng.uniform(-2e-9, 2e-9, (n, nlev)),
                dprcl=conv, dprscl=snow_c, dprll=large, dprsll=snow_l)


def frozen(s):
    """The snow each layer forms, convective and large-scale: all that a cold layer forms, the snow formed as snow in
    a warm one (where it melts)."""
    cold = s['dt'][:, :-1] <= C['TMELT']
    return np.where(cold, s['dprcl'], s['dprscl']), np.where(cold, s['dprll'], s['dprsll'])


def run_python(s, snow=SNOW):
    py = copy(s)
    out = rd.mkrain(py, C, rd.PARAMETERS, snow)
    return py, out


# ------------------------------------------------------------------------------------------------------ Python

def _one_column(h1, top, f1, prm):
    st = ss.snowprep(dict(dq=np.zeros((1, 2)), dt=np.full((1, 3), 250.0), sigma=np.array([0.3, 0.8]),
                          dsigma=np.array([0.6, 0.4]), dp=np.array([1e5])), C, prm)
    st['zh'][0, 1] = h1
    st['topw'][0], st['topz'][0], st['high'][0] = 1.0, top, top
    st['mass'] = np.array([1.0])
    st['fl'][0, :, 0] = f1
    ss._start(st, 0, 1, prm)
    return st


@pytest.mark.parametrize('h1,top,rate_mm_h,sub', [(4100.0, 24000.0, 0.3, 0.05), (3957.0, 16280.0, 0.6, 0.053),
                                                  (8254.0, 30000.0, 0.05, 0.04), (6000.0, 40000.0, 0.01, 0.02),
                                                  (1200.0, 9000.0, 2.0, 0.08)])
def test_the_schemes_closure_gives_closure_columns_flux_under_every_matching_and_both_tables(h1, top, rate_mm_h, sub):
    prm = dict(SNOW, snowstep=500.0)
    f1 = rate_mm_h / 3600.0
    fams = ss.load()['families']
    f10 = ss.load(sensitivity='10cm')['families']
    edges = np.array([[h1, h1 - 0.3 * h1, h1 - 0.65 * h1, 0.0]])
    st = _one_column(h1, top, f1, prm)
    for a, b in zip(edges[0, :-1], edges[0, 1:]):
        ss._integrate(st, 0, a, b, sub, prm)
    want = lambda fam, sigma: ss.closure_columns(fam, np.array([h1]), np.array([f1]), np.array([top]), edges,
                                                 np.full((1, 3), sub), sigma=sigma, step=500.0, exact=True)[0, -1] / f1
    mu = max(prm['snowlnc'] + prm['snowbeta'] * math.log(f1), math.log(f1))
    spread = lambda fam: sum(w * want(fam, f1 / math.exp(mu + prm['snowspread'] * x))
                             for x, w in zip(ss.GAUSS_X, ss.GAUSS_W))
    assert ss._kept(st, 0, 0) == pytest.approx(spread(fams), rel=1e-12)
    assert ss._kept(st, 0, 1) == pytest.approx(want(fams, 1.0), rel=1e-12)
    assert ss._kept(st, 0, 2) == pytest.approx(want(fams, prm['snowfixed']), rel=1e-12)
    assert ss._kept(st, 0, 3) == pytest.approx(spread(f10), rel=1e-12)


def test_the_constants_come_from_the_tests_and_the_products():
    assert ss.PARAMETERS['snowgconv'] == 0.75 and ss.PARAMETERS['snowglarge'] == 0.4
    assert (ss.FRACTION_LN_C, ss.FRACTION_BETA, ss.FRACTION_SPREAD) == pytest.approx((-2.8933, 0.6273, 1.0468),
                                                                                       abs=1e-4)
    assert ss.FIXED_FRACTION == pytest.approx(0.5 * (0.5130966962394068 + 0.5822050450211864))
    assert sum(ss.GAUSS_W) == pytest.approx(1.0) and sum(w * x * x for x, w in zip(ss.GAUSS_X, ss.GAUSS_W)) == \
        pytest.approx(1.0)
    t = ss.TABLES
    assert t['tops'][0] * t['lapse'] == pytest.approx(15.0, abs=0.01) and len(t['refs']) == 11
    assert [t['refs'][r] for r in t['sets'][1]][2] == 'capl_ref_q1_f1'
    text = ss.declarations()
    for value in (19.3, 0.37, 400.0, 0.86, 0.28, 21.8745584, 7.66, ss.RHOSU, ss.GAMMA_4B, ss.GAMMA_VENT):
        assert repr(value) in text, value
    assert max(len(line) for line in (text + ss.ROUTINES).splitlines()) <= 132    # gfortran's free-form limit


def test_the_column_keeps_its_water_and_sublimates_only_air_unsaturated_over_ice():
    s = snow_columns()
    py, out = run_python(s)
    d = out['snow']
    fc, fl = frozen(s)
    formed = (fc + fl).sum(axis=1) * s['dp'] / C['ga'] / 1000.0
    assert d[:, 0] == pytest.approx(formed, rel=1e-12)
    for v in range(ss.NV):                          # a warm ground: all of it sublimates or melts
        sub, melt = ss.SLOTS[v]
        lost = d[:, list(sub)].sum(axis=1) + d[:, list(melt)].sum(axis=1) if v == 0 else d[:, sub] + d[:, melt]
        if ss.VARIANTS[v][1] == 2:                  # graupel alone: its share of the snow formed
            share = fc * ss.GRAUPEL_SHARE_CONVECTIVE + fl * ss.GRAUPEL_SHARE_LARGE_SCALE
            formed_v = share.sum(axis=1) * s['dp'] / C['ga'] / 1000.0
        else:
            formed_v = formed
        assert lost == pytest.approx(formed_v, rel=1e-9, abs=1e-20), ss.VARIANTS[v][0]
    assert np.all(out['dprs'] == 0.0)
    assert np.all(out['snow_sublimation'] >= 0.0) and out['snow_sublimation'].sum() > 0.0
    column = (out['snow_sublimation'] * s['dsigma']).sum(axis=1) * s['dp'] / C['ga'] / 1000.0
    assert column == pytest.approx(d[:, 1] + d[:, 2], rel=1e-12, abs=1e-25)
    saturated = copy(s)
    qsi = np.array([[ss.ice_saturation(t, s['dp'][j] * s['sigma'][k], C)[0]
                     for k, t in enumerate(s['dt'][j, :-1])] for j in range(len(s['dp']))])
    saturated['dq'] = np.maximum(s['dq'], 1.02 * qsi)
    saturated['dqdt'] = np.zeros_like(s['dqdt'])
    _, wet = run_python(saturated)
    assert np.all(wet['snow'][:, 1:3] == 0.0) and np.all(wet['snow_sublimation'] == 0.0)


def test_the_layers_that_form_snow_do_not_sublimate_and_the_graupel_follows_its_kind():
    s = snow_columns(seed=22)
    py, out = run_python(s)
    fc, fl = frozen(s)
    forming = (fc + fl) > 0.0
    assert np.all(out['snow_sublimation'][forming] == 0.0)
    assert np.any(out['snow_sublimation'][~forming] > 0.0)
    assert out['snow'][:, 1].sum() > 0.0 and out['snow'][:, 2].sum() > 0.0       # flakes and graupel both sublimate
    cold = s['dt'][:, :-1] <= C['TMELT']
    only_large = copy(s)                                  # the cold layers' convective snow made large-scale
    only_large['dprll'] = s['dprll'] + np.where(cold, s['dprcl'], 0.0)
    only_large['dprcl'] = np.where(cold, 0.0, s['dprcl'])
    _, a = run_python(only_large)
    _, b = run_python(s)
    assert a['snow'][:, 0] == pytest.approx(b['snow'][:, 0])          # the same snow formed, less of it graupel
    assert a['snow'][:, 4].sum() + a['snow'][:, 2].sum() < b['snow'][:, 4].sum() + b['snow'][:, 2].sum()


def test_a_column_colder_than_freezing_to_the_ground_takes_its_freezing_level_below_the_ground():
    s = snow_columns(seed=23, cold_ground=True)
    cold = np.arange(len(s['dp'])) % 2 == 0
    s['dt'][cold] = np.minimum(s['dt'][cold], 272.0)
    s['dt'][cold, -1] = 263.0
    st = ss.snowprep(s, C, SNOW)
    want = -(C['TMELT'] - s['dt'][cold, -1]) / ss.TABLES['lapse']
    assert st['z0'][cold] == pytest.approx(want)
    assert np.all(st['zh'][cold, -1] > 0.0)
    _, out = run_python(s)
    fc, fl = frozen(s)
    assert np.all(out['dprs'][cold & ((fc + fl).sum(axis=1) > 0.0)] > 0.0)


def test_with_nsnowsub_0_mkrain_is_mkrain_with_gamma():
    s = snow_columns(seed=24, cold_ground=True)
    a, oa = run_python(s, None)
    b, ob = run_python(s, dict(SNOW, nsnowsub=0))
    for key in ('dtdt', 'dqdt'):
        assert np.array_equal(a[key], b[key])
    for key in ('dprc', 'dprl', 'dprs'):
        assert np.array_equal(oa[key], ob[key])
    assert ob['snow'] is None


def test_the_shadows_never_feed_back():
    s = snow_columns(seed=25, cold_ground=True, warm_layer=True)
    a, oa = run_python(s, dict(SNOW, nsnowshad=0))
    b, ob = run_python(s, dict(SNOW, nsnowshad=1))
    for key in ('dtdt', 'dqdt'):
        assert np.array_equal(a[key], b[key])
    assert np.array_equal(oa['snow'][:, :5], ob['snow'][:, :5])
    assert np.all(oa['snow'][:, 5:19] == 0.0) and np.any(ob['snow'][:, 5:19] > 0.0)


# ------------------------------------------------------------------------------------------------------ Fortran

STUB = """
      module rainmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLEM={nlev}-1, NLON=1, NLAT=1, NROOT=0
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV)
      real :: deltsec, deltsec2, ga, gascon, acpd, adv, alv, als, rdbrv, ra1, ra2, ra4, tmelt, gamma
      real :: dprll(NHOR,NLEV), dprsll(NHOR,NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV)
      real :: dprc(NHOR), dprl(NHOR), dprs(NHOR)
      real :: dgp3d(NHOR,NLEV,25)
      real :: dentropy(NHOR,40), dentro(NHOR), dentrot(NHOR,NLEV), dentroq(NHOR,NLEV), dentrop(NHOR)
      real :: dentro3d(NHOR,NLEV,40), denergy(NHOR,40), dener3d(NHOR,NLEV,40)
      integer :: ndiaggp=1, nentropy=0, nentro3d=0, nenergy=0, nener3d=0, nprc=2
      integer :: nevapprec= 1
      integer, parameter :: NSNOWD = {ndiag}
      integer :: nsnowout = 0, nsnacc = 0
      real :: asnow(NHOR,NSNOWD) = 0.
{declarations}
      end module rainmod
"""
HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: jhor, jlev, jd
      real :: zt0(NHOR,NLEP), zq0(NHOR,NLEP)
      read(*,*) deltsec, deltsec2, ga, gascon, acpd, adv, alv, als, rdbrv, ra1, ra2, ra4, tmelt, gamma
      read(*,*) nevapprec, {rain}
      read(*,*) {snow}
      read(*,*) sigma, dsigma
      do jhor=1,NHOR
       read(*,*) dp(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEP)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dtdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprcl(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprscl(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprll(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprsll(jhor,jlev), jlev=1,NLEV)
      enddo
      zt0(:,:)=dtdt(:,:)
      zq0(:,:)=dqdt(:,:)
      call mkrain
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') (dtdt(jhor,jlev)-zt0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dqdt(jhor,jlev)-zq0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') dprc(jhor), dprl(jhor), dprs(jhor)
       write(*,'(*(es25.16e3))') (asnow(jhor,jd), jd=1,NSNOWD)
      enddo
      end program harness
"""


def run_fortran(folder, code, s, snow):
    n, nlev = s['dq'].shape
    folder.mkdir()
    decl = rd.declarations(rd.PARAMETERS) + ss.declarations(snow)
    files = {'stub.f90': STUB.format(nhor=n, nlev=nlev, ndiag=ss.NDIAG, declarations=decl), 'rain.f90': code,
             'harness.f90': HARNESS.replace('{rain}', ', '.join(rd.NAMELIST_NAMES))
             .replace('{snow}', ', '.join(ss.NAMELIST_NAMES))}
    for name, body in files.items():
        (folder / name).write_text(body)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    lines = [num([C[k] for k in ('deltsec', 'deltsec2', 'ga', 'gascon', 'acpd', 'ADV', 'ALV', 'ALS', 'rdbrv', 'ra1',
                                 'ra2', 'ra4', 'TMELT', 'gamma')]),
             '2 ' + ' '.join(repr(rd.PARAMETERS[k]) for k in rd.NAMELIST_NAMES),
             ' '.join(repr(snow[k]) for k in ss.NAMELIST_NAMES), num(np.r_[s['sigma'], s['dsigma']])]
    for j in range(n):
        lines += [num(s['dp'][j])] + [num(s[k][j]) for k in ('dt', 'dq', 'dtdt', 'dqdt', 'dprcl', 'dprscl', 'dprll',
                                                               'dprsll')]
    text = subprocess.run([str(folder / 'harness')], input='\n'.join(lines) + '\n', capture_output=True, text=True,
                          check=True).stdout
    rows = [np.array(r.split(), float) for r in text.strip('\n').split('\n')]
    pick = lambda i: np.array([rows[4 * j + i] for j in range(n)])
    return dict(dtdt=pick(0), dqdt=pick(1), ground=pick(2), snow=pick(3))


def patched_code():
    """The patched mkrain with the rain's and the snow's routines after it."""
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = apply(path.read_text(), rd.RAIN_DROP_EDITS, rd.MARKER)
    start = text.index('      subroutine mkrain\n')
    end = '      end subroutine snowsubl\n'
    return text[start:text.index(end) + len(end)]


def compare(got, py, s, out, label):
    scale = lambda a: max(np.abs(a).max(), 1e-30)
    for key, want in (('dtdt', py['dtdt'] - s['dtdt']), ('dqdt', py['dqdt'] - s['dqdt']),
                      ('ground', np.c_[out['dprc'], out['dprl'], out['dprs']]), ('snow', out['snow'])):
        assert np.isfinite(got[key]).all() and np.isfinite(want).all(), (label, key)
        if key == 'snow':
            for d in range(want.shape[1]):
                assert np.abs(got[key][:, d] - want[:, d]).max() <= 1e-10 * scale(want[:, d]), (label, key, d)
        else:
            assert np.abs(got[key] - want).max() <= 1e-10 * scale(want), (label, key)


@FORTRAN
@pytest.mark.parametrize('seed,cold,warm,top,shadows,plasim', [(31, False, False, 1, 1, True),
                                                                (32, True, False, 1, 1, True),
                                                                (33, True, True, 1, 1, True),
                                                                (34, True, False, 2, 1, True),
                                                                (35, False, False, 1, 0, True),
                                                                (37, True, True, 1, 1, False)])
def test_the_fortran_gives_the_python_copys_snow(tmp_path, seed, cold, warm, top, shadows, plasim):
    code = patched_code()
    snow = dict(SNOW, nsnowtop=top, nsnowshad=shadows)
    s = snow_columns(seed=seed, cold_ground=cold, warm_layer=warm, plasim=plasim)
    got = run_fortran(tmp_path / f'snow{seed}', code, s, snow)
    py, out = run_python(s, snow)
    compare(got, py, s, out, seed)


@FORTRAN
def test_in_dry_air_and_heavy_snow_the_fortran_stays_finite_and_follows_the_copy(tmp_path):
    code = patched_code()
    s = snow_columns(seed=36, cold_ground=True, warm_layer=True, dry=True)
    for key in ('dprcl', 'dprscl', 'dprll', 'dprsll'):
        s[key] = s[key] * 30.0
    for step in (200.0, 3000.0):
        snow = dict(SNOW, snowstep=step)
        got = run_fortran(tmp_path / f'dry{step:g}', code, s, snow)
        py, out = run_python(s, snow)
        compare(got, py, s, out, step)


@FORTRAN
def test_the_patched_modules_compile_with_plasims_own(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = path.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    (tmp_path / 'plasimmod.f90').write_text(apply((src / 'plasimmod.f90').read_text(), ss.OUTPUT_MOD_EDITS,
                                                  ss.OUTPUT_MARKER))
    shutil.copy(src / 'restartmod.f90', tmp_path / 'restartmod.f90')
    (tmp_path / 'rainmod.f90').write_text(apply(path.read_text(), rd.RAIN_DROP_EDITS, rd.MARKER))
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90', 'restartmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def test_the_output_patches_apply_once_and_undo_on_lines_no_other_patch_touches():
    for name, marker, edits in ss.OUTPUT_PATCHES:
        path = source(name)
        if path is None:
            pytest.skip('no ExoPlaSim source at hand')
        plain = undo(path.read_text(), edits, marker)
        patched = apply(plain, edits, marker)
        assert apply(patched, edits, marker) == patched and undo(patched, edits, marker) == plain
    import climate.gcm.condensate as cnd
    import climate.gcm.massflux_fortran as mff
    others = [e for _, _, edits in er.PATCHES + er.SBM_DEPTH_PATCHES + er.SBM_CAP_PATCHES + mff.MASSFLUX_PATCHES
              + cnd.CLOUD_WATER_PATCHES + er.HUS_MEAN_PATCHES + er.WIND10_PATCHES + er.PLASIM_SWITCH_PATCHES
              + [x for x in rd.RAIN_DROP_PATCHES if x not in ss.OUTPUT_PATCHES] for e in edits]
    for _, _, edits in ss.OUTPUT_PATCHES:
        for old, _ in edits:
            assert all(old not in o and old not in n for o, n in others), old
