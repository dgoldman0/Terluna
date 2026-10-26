"""The mass-flux scheme's Fortran routine for PlaSim (massflux_fortran.py): the patch applies once to both installed
ExoPlaSim sources over the Betts-Miller patch, with or without that scheme's options, and compiles; the column routine
compiled alone gives the Python copy's tendencies, rain, snow and mass fluxes (massflux.column with PlaSim's phase
rule) to 1e-10 on lunar and Earth-like columns, and keeps the column's energy and water; deep, shallow and mid-level
columns convect as designed; the PlaSim routine hands mkrain its water, which the patched mkrain carries to the
ground without loss and without a negative running total, and with nprc 1 or 2 the patched mkrain is the old one bit
for bit; the routine's defaults are the Python copy's parameters; and the runner records the scheme with every
parameter, writes its namelist and refuses the Betts-Miller scheme's options beside it."""
import math
import re
import shutil
import signal
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import massflux as mf
from climate.gcm import massflux_fortran as mff

GCM_FLAGS = ('-O3', '-cpp', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffpe-summary=none',
             '-finit-real=zero', '-fdefault-real-8')
FORTRAN = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')
G_MOON = 1.6242188765608128
LUNAR = dict(mf.PLASIM, ga=G_MOON, gascon=288.3767, acpd=288.3767 / 0.286, rdbrv=288.3767 / 461.51, ADV=0.8546)
EARTH = dict(mf.PLASIM, ga=9.80665)
CONSTANTS = ('ALV', 'ALS', 'acpd', 'ADV', 'TMELT', 'ra1', 'ra2', 'ra4', 'rdbrv', 'gascon')
num = lambda xs: ' '.join(repr(float(x)) for x in xs)


# ------------------------------------------------------------------------------------------------------- columns

def sounding(g, ps, t0, n, top_m, lapse_bl, lapse, t_min, rh_bl, rh, stable=(0.0, 0.0), seed=0, c=LUNAR,
             levels=None, noise=0.05):
    """A column top first on n sigma levels (levels: half levels in sigma, the top layer's bottom first, or n equal
    in height to top_m): a well-mixed boundary layer at lapse_bl K/m from t0 at the ground, its vapour rh_bl of
    saturation at the ground, up to its condensation level; above, the temperature falls lapse K/m to t_min at rh
    of saturation, with stable = (K/m, m) warming the first metres above the condensation level by that much more;
    noise (K, and a fifth of it as a share of the vapour) from seed. Hydrostatic at gravity g."""
    rng = np.random.default_rng(seed)
    fine = np.linspace(0.0, top_m, 20001)
    temp, lnp = np.empty(fine.size), np.empty(fine.size)
    lnp[0], temp[0] = math.log(ps), t0
    q0 = rh_bl * mf.qsat(t0, ps, c)
    lcl = None
    for i in range(1, fine.size):
        dz = fine[i] - fine[i - 1]
        lnp[i] = lnp[i - 1] - dz * g / (c['gascon'] * temp[i - 1])
        if lcl is None:
            temp[i] = temp[i - 1] - lapse_bl * dz
            if mf.qsat(temp[i], math.exp(lnp[i]), c) <= q0:
                lcl = fine[i]
        elif fine[i] < lcl + stable[1]:
            temp[i] = temp[i - 1] - (lapse - stable[0]) * dz
        else:
            temp[i] = max(temp[i - 1] - lapse * dz, t_min)
    if levels is None:                              # equal layers in height, the top half level at nil pressure
        ph = np.r_[0.0, np.exp(np.interp(np.linspace(0.0, top_m * 0.999, n + 1), fine, lnp))[::-1][1:]]
    else:
        ph = np.r_[0.0, np.asarray(levels, float) * ps]
    pf = 0.5 * (ph[:-1] + ph[1:])                   # PlaSim's full levels, midway in sigma
    zf = np.interp(-np.log(pf), -lnp, fine)
    t = np.interp(zf, fine, temp) + rng.normal(0.0, noise, pf.size)
    q = np.where(zf < lcl, q0, [rh * mf.qsat(a, b, c) for a, b in zip(t, pf)])
    q = q * (1.0 + rng.normal(0.0, 0.2 * noise, pf.size))
    return t, np.maximum(q, 1.0e-9), pf, ph


def lunar_columns():
    """Lunar columns on PlaSim's fifteen levels and on 40 equal layers: deep, capped, dry and moist."""
    half15 = mf.grid_sets()['15']
    out = []
    for j, (kw, lev) in enumerate([
            (dict(t0=300.0, rh_bl=0.75, rh=0.8), half15),
            (dict(t0=302.0, rh_bl=0.8, rh=0.85), half15),
            (dict(t0=300.0, rh_bl=0.8, rh=0.6, stable=(2.0e-3, 6000.0)), half15),
            (dict(t0=296.0, rh_bl=0.6, rh=0.5), half15),
            (dict(t0=300.0, rh_bl=0.75, rh=0.8), None),
            (dict(t0=301.0, rh_bl=0.8, rh=0.7, stable=(1.5e-3, 6000.0)), None),
            (dict(t0=280.0, rh_bl=0.9, rh=0.9), half15),
            (dict(t0=300.0, rh_bl=0.2, rh=0.1, stable=(3.0e-3, 9000.0)), half15)]):
        t, q, pf, ph = sounding(G_MOON, 121500.0, n=15 if lev is not None else 40, top_m=150000.0, lapse_bl=1.6e-3,
                                lapse=1.0e-3, t_min=200.0, seed=j, c=LUNAR, levels=lev, **kw)
        out.append((t, q, pf, ph, G_MOON, LUNAR))
    return out


def earth_columns():
    """Earth-like columns at Earth's gravity on PlaSim's ten levels and on 30 equal layers."""
    half10 = mf.grid_sets()['10']
    out = []
    for j, (kw, lev) in enumerate([(dict(t0=300.0, rh_bl=0.8, rh=0.75), half10),
                                   (dict(t0=302.0, rh_bl=0.85, rh=0.8), None),
                                   (dict(t0=298.0, rh_bl=0.8, rh=0.5, stable=(8.0e-3, 2000.0)), None),
                                   (dict(t0=290.0, rh_bl=0.5, rh=0.4), half10)]):
        t, q, pf, ph = sounding(9.80665, 101325.0, n=10 if lev is not None else 30, top_m=30000.0, lapse_bl=9.8e-3,
                                lapse=6.5e-3, t_min=200.0, seed=10 + j, c=EARTH, levels=lev, **kw)
        out.append((t, q, pf, ph, 9.80665, EARTH))
    return out


def forcing(n, seed, land, g):
    """Surface fluxes and the rest of a column's inputs, varied by seed: (shf, lhf, bl_dtdt, bl_dqdt, omega, wind,
    land, the temperature offset of mkrain's phase temperature)."""
    rng = np.random.default_rng(100 + seed)
    shf, lhf = rng.uniform(5.0, 120.0), rng.uniform(20.0, 250.0)
    bt = rng.normal(0.0, 1.0e-5 * g / G_MOON, n)
    bq = rng.normal(0.0, 1.0e-8 * g / G_MOON, n)
    om = rng.normal(-0.02, 0.05, n) * (g / G_MOON)
    wind = np.abs(rng.normal(4.0, 2.0, n))
    return shf, lhf, bt, bq, om, wind, land, rng.normal(0.0, 0.3)


def mid_level_column():
    """A dry ground under a saturated layer in large-scale ascent, on the fifteen levels."""
    half15 = mf.grid_sets()['15']
    t, q, pf, ph = sounding(G_MOON, 121500.0, t0=290.0, n=15, top_m=150000.0, lapse_bl=1.6e-3, lapse=1.0e-3,
                            t_min=200.0, rh_bl=0.3, rh=0.99, levels=half15)
    om = np.where(pf > 0.5 * ph[-1], -0.05, 0.0)
    return dict(t=t, q=q, pf=pf, ph=ph, g=G_MOON, c=LUNAR, dt=3600.0, land=0, shf=-5.0, lhf=0.0, bt=np.zeros(15),
                bq=np.zeros(15), om=om, w=np.full(15, 3.0), tph=t)


def all_cases():
    cases = [mid_level_column()]
    for j, (t, q, pf, ph, g, c) in enumerate(lunar_columns() + earth_columns()):
        shf, lhf, bt, bq, om, wind, land, dtph = forcing(t.size, j, j % 2, g)
        cases.append(dict(t=t, q=q, pf=pf, ph=ph, g=g, c=c, dt=3600.0, shf=shf, lhf=lhf, bt=bt, bq=bq, om=om, w=wind,
                          land=land, tph=t + dtph))
        if j < 3:                                   # the same column under sinking air and cooled from below
            cases.append(dict(cases[-1], om=np.abs(om), shf=-10.0, lhf=5.0, land=1 - land))
    return cases


# --------------------------------------------------------------------------------------------- the column harness

COLUMN_STUB = """
      module rainmod
      real :: ALV, ALS, acpd, ADV, TMELT, ra1, ra2, ra4, rdbrv, gascon
"""
COLUMN_HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: kn, kland, ncase, jc, kinfo(4)
      real :: zg, zdt, zshf, zlhf, zinfo(10)
      real, allocatable :: zt(:), zq(:), zpf(:), zph(:), zbt(:), zbq(:), zom(:), zw(:), ztph(:)
      real, allocatable :: zdtdt(:), zdqdt(:), zr(:), zs(:), zmu(:), zmd(:)
      read(*,*) ncase
      do jc = 1, ncase
       read(*,*) ALV, ALS, acpd, ADV, TMELT, ra1, ra2, ra4, rdbrv, gascon
       read(*,*) kn, kland, zg, zdt, zshf, zlhf
       allocate(zt(kn), zq(kn), zpf(kn), zph(0:kn), zbt(kn), zbq(kn), zom(kn), zw(kn), ztph(kn))
       allocate(zdtdt(kn), zdqdt(kn), zr(kn), zs(kn), zmu(0:kn), zmd(0:kn))
       read(*,*) zt
       read(*,*) zq
       read(*,*) zpf
       read(*,*) zph
       read(*,*) zbt
       read(*,*) zbq
       read(*,*) zom
       read(*,*) zw
       read(*,*) ztph
       call mfcolumn(kn,zt,zq,zpf,zph,zg,zdt,zshf,zlhf,zbt,zbq,zom,zw,kland,ztph,zdtdt,zdqdt,zr,zs,zmu,zmd,kinfo,zinfo)
       write(*,'(4i6)') kinfo
       write(*,'(*(es25.16e3))') zinfo
       write(*,'(*(es25.16e3))') zdtdt
       write(*,'(*(es25.16e3))') zdqdt
       write(*,'(*(es25.16e3))') zr
       write(*,'(*(es25.16e3))') zs
       write(*,'(*(es25.16e3))') zmu
       write(*,'(*(es25.16e3))') zmd
       deallocate(zt, zq, zpf, zph, zbt, zbq, zom, zw, ztph, zdtdt, zdqdt, zr, zs, zmu, zmd)
      enddo
      end program harness
"""
OUTPUTS = ('dtdt', 'dqdt', 'rain', 'snow', 'mu', 'md')


def compile_in(folder, files: dict, name='harness'):
    folder.mkdir(parents=True, exist_ok=True)
    for fname, text in files.items():
        (folder / fname).write_text(text)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', name, *files], cwd=folder, check=True, capture_output=True)
    return folder / name


def column_harness(folder):
    stub = COLUMN_STUB + mff.declarations() + '      end module rainmod\n'
    return compile_in(folder, {'stub.f90': stub, 'mfcolumn.f90': mff.COLUMN_CODE, 'harness.f90': COLUMN_HARNESS})


def run_columns(exe, cases) -> list:
    lines = [str(len(cases))]
    for x in cases:
        lines.append(num([x['c'][k] for k in CONSTANTS]))
        lines.append(f"{x['t'].size} {x['land']} {x['g']!r} {x['dt']!r} {x['shf']!r} {x['lhf']!r}")
        lines += [num(x[k]) for k in ('t', 'q', 'pf', 'ph', 'bt', 'bq', 'om', 'w', 'tph')]
    rows = subprocess.run([str(exe)], input='\n'.join(lines) + '\n', capture_output=True, text=True,
                          check=True).stdout.strip('\n').split('\n')
    assert len(rows) == 8 * len(cases)
    out = []
    for j in range(len(cases)):
        r = rows[8 * j:8 * j + 8]
        res = dict(kinfo=[int(v) for v in r[0].split()], info=np.array(r[1].split(), float))
        res.update({k: np.array(r[2 + i].split(), float) for i, k in enumerate(OUTPUTS)})
        out.append(res)
    return out


def python_column(x):
    return mf.column(x['t'], x['q'], x['pf'], x['ph'], x['g'], x['dt'], x['shf'], x['lhf'], x['bt'], x['bq'], x['om'],
                     x['w'], bool(x['land']), c=x['c'], melt='plasim', t_phase=x['tph'])


TYPES = {None: 0, 'deep': 1, 'shallow': 2, 'mid': 3}


@FORTRAN
def test_the_fortran_column_gives_the_python_copys_tendencies(tmp_path):
    cases = all_cases()
    fortran = run_columns(column_harness(tmp_path / 'column'), cases)
    kinds = set()
    for j, (x, f) in enumerate(zip(cases, fortran)):
        p = python_column(x)
        kinds.add((p['type'], x['g'] > 5.0))
        assert f['kinfo'][0] == TYPES[p['type']], j
        if p['type'] is None:
            assert all(not np.any(f[k]) for k in OUTPUTS), j
            continue
        assert f['kinfo'][1] == p['up']['base_half'] and f['kinfo'][2] == p['up']['top_half'], j
        assert f['kinfo'][3] == int(p['rooted']), j
        assert f['info'][0] == pytest.approx(p['m_base'], rel=1e-10), j
        assert f['info'][2] == pytest.approx(p['z_top'], rel=1e-10)
        assert f['info'][5] == pytest.approx(p['pcape'], rel=1e-10, abs=1e-12)
        for k, pk in zip(OUTPUTS, ('dtdt', 'dqdt', 'net_rain', 'net_snow', 'mu', 'md')):
            b = np.asarray(p[pk], float)
            assert np.allclose(f[k], b, rtol=1e-10, atol=1e-10 * np.abs(b).max()), (j, k)
    # the cases reach every type on the Moon and deep convection at Earth's gravity
    assert {('deep', False), ('shallow', False), ('mid', False), (None, False), ('deep', True)} <= kinds


@FORTRAN
def test_the_fortran_column_keeps_its_energy_and_water(tmp_path):
    cases = all_cases()
    fortran = run_columns(column_harness(tmp_path / 'column'), cases)
    rained = 0
    for x, f in zip(cases, fortran):
        if f['kinfo'][0] == 0:
            continue
        c = x['c']
        m = np.diff(x['ph']) / x['g']
        heat = float(np.sum(c['acpd'] * (1.0 + c['ADV'] * x['q']) * f['dtdt'] * m))
        water = float(np.sum(f['dqdt'] * m))
        rain, snow = float(f['rain'].sum()), float(f['snow'].sum())
        assert water == pytest.approx(-(rain + snow), rel=1e-9, abs=1e-15)
        assert heat == pytest.approx(c['ALV'] * rain + c['ALS'] * snow, rel=1e-9, abs=1e-9)
        rained += rain + snow > 0.0
    assert rained >= 4


@FORTRAN
def test_deep_shallow_and_midlevel_columns_convect_as_designed(tmp_path):
    half15 = mf.grid_sets()['15']
    base = dict(g=G_MOON, c=LUNAR, dt=3600.0, land=0)
    cols = {}
    # deep: a moist, conditionally unstable column heated from below
    t, q, pf, ph = sounding(G_MOON, 121500.0, t0=301.0, n=15, top_m=150000.0, lapse_bl=1.6e-3, lapse=1.0e-3,
                            t_min=200.0, rh_bl=0.8, rh=0.8, levels=half15)
    cols['deep'] = dict(base, t=t, q=q, pf=pf, ph=ph, shf=50.0, lhf=150.0, bt=np.zeros(15), bq=np.zeros(15),
                        om=np.zeros(15), w=np.full(15, 3.0), tph=t)
    # shallow: a dry-adiabatic boundary layer under a layer 0.4 K/km more stable for 6 km above its condensation
    # level, on layers 1 km deep
    t, q, pf, ph = sounding(G_MOON, 121500.0, t0=300.0, n=150, top_m=150000.0, lapse_bl=G_MOON / LUNAR['acpd'],
                            lapse=1.0e-3, t_min=220.0, rh_bl=0.75, rh=0.8, stable=(0.4e-3, 6000.0), noise=0.0)
    bq = np.zeros(150)
    bq[-1] = 2.0e-8
    cols['shallow'] = dict(base, t=t, q=q, pf=pf, ph=ph, shf=10.0, lhf=150.0, bt=np.zeros(150), bq=bq,
                           om=np.zeros(150), w=np.full(150, 3.0), tph=t)
    # mid-level: a dry ground under a saturated layer in large-scale ascent
    cols['mid'] = mid_level_column()
    cols['dry'] = dict(cols['mid'], q=0.05 * cols['mid']['q'], om=np.zeros(15))
    out = dict(zip(cols, run_columns(column_harness(tmp_path / 'column'), list(cols.values()))))
    assert [out[k]['kinfo'][0] for k in ('deep', 'shallow', 'mid', 'dry')] == [1, 2, 3, 0]
    deep, shallow, mid = out['deep'], out['shallow'], out['mid']
    # deep: rooted at the ground, its cloud over 200 hPa thick, raining and heating the column aloft
    assert deep['kinfo'][3] == 1 and deep['info'][3] - deep['info'][4] > mf.PARAMETERS['deep_dp']
    assert deep['rain'].sum() + deep['snow'].sum() > 0.0 and deep['dtdt'][:10].max() > 0.0
    # shallow: its top within the stable layer, its cloud thinner than 200 hPa
    zh_base, z_top = shallow['info'][1], shallow['info'][2]
    assert z_top - zh_base < 6000.0 + 3000.0 and shallow['info'][3] - shallow['info'][4] <= mf.PARAMETERS['deep_dp']
    assert shallow['info'][0] > 0.0
    # mid-level: its base mass flux the large-scale ascent's at its departure level, -omega / g
    assert mid['info'][0] == pytest.approx(0.05 / G_MOON, rel=1e-12) or mid['info'][8] == 1.0
    assert mid['kinfo'][3] == 0
    assert not np.any(out['dry']['dtdt'])


# ------------------------------------------------------------------------------------ the patch and PlaSim's routines

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


def patched(text, optional=(), first=False):
    """A source with the Betts-Miller patch, and the mass-flux patch before or after the optional patches."""
    text = apply(text, er.SBM_EDITS, er.SBM_MARKER)
    if first:
        text = apply(text, mff.MASSFLUX_EDITS, mff.MASSFLUX_MARKER)
    for edits, marker in optional:
        text = apply(text, edits, marker)
    return apply(text, mff.MASSFLUX_EDITS, mff.MASSFLUX_MARKER)


OPTIONAL = ((er.SBM_TRIGGER_EDITS, er.SBM_TRIGGER_MARKER), (er.SBM_DEPTH_EDITS, er.SBM_DEPTH_MARKER),
            (er.SBM_CAP_EDITS, er.SBM_CAP_MARKER))


def test_the_patch_applies_once_to_both_sources_with_or_without_the_other_options():
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    marker = mff.MASSFLUX_MARKER
    assert not any(marker in old for old, _ in mff.MASSFLUX_EDITS)
    # every edit adds text after what it finds, but for mkrain's two evaporations, whose arithmetic under nprc 1 and 2
    # test_with_nprc_1_or_2_mkrain_is_as_it_was checks
    changed = [old for old, new in mff.MASSFLUX_EDITS if not new.startswith(old)]
    assert changed == [old for old in changed if ',zprc(:))' in old or ',zprsc(:))' in old] and len(changed) == 2
    assert all(len(line) <= 132 for _, new in mff.MASSFLUX_EDITS for line in new.splitlines())
    for path in sources:
        plain = path.read_text()
        versions = [patched(plain)]
        if er.SBM_TRIGGER_MARKER not in plain:
            versions += [patched(plain, OPTIONAL), patched(plain, OPTIONAL, first=True)]
        for text in versions:
            assert apply(text, mff.MASSFLUX_EDITS, marker) == text                      # once only
            assert text.count('      if(nprc == 3) call mkmassflux') == 1
            assert '      if(nprc == 1) call kuo\n      if(nprc == 2) call mksbm' in text
            assert text.count('subroutine mfcolumn') == 2 and text.count('subroutine mkmassflux') == 2
            assert text.count('call mpbcr(epsbg_mf)') == 1 and text.count('zrem(:)=MAX(zfloor(:,jlev)') == 1
            assert text.count(',zprc(:),zrem(:))') == 1 and text.count(',zprsc(:),zrem(:))') == 1
            for name, _, _, _ in mff.PARAMETER_NAMES:                                 # declared, read, broadcast
                assert len(re.findall(rf'(real|integer) :: {name} = ', text)) == 1
                assert len(re.findall(rf'namelist/rainmod_nl/.*\b{name}\b', text)) == 1
                assert len(re.findall(rf'call mpbc[ri]\({name}\)', text)) == 1
        if len(versions) == 3:                      # either order gives the same lines
            assert sorted(versions[1].split('\n')) == sorted(versions[2].split('\n'))


@FORTRAN
def test_the_patched_source_compiles_with_plasims_modules(tmp_path):
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    src = sources[0].parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    shutil.copy(src / 'plasimmod.f90', tmp_path / 'plasimmod.f90')
    (tmp_path / 'rainmod.f90').write_text(patched(sources[0].read_text(), OPTIONAL))
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def mkrain_text(text):
    return text[text.index('      subroutine mkrain\n'):text.index('      end subroutine mkrain\n')] + \
        '      end subroutine mkrain\n'


# PlaSim's module as mkmassflux and mkrain read it, for NHOR columns on NLEV levels; mkdttgp hands back a test array.
PLASIM_STUB = """
      module rainmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLEM={nlev}-1
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dqt(NHOR,NLEP), dp(NHOR)
      real :: du(NHOR,NLEP), dv(NHOR,NLEP), dw(NHOR,NLEV), dls(NHOR), dshfl(NHOR), dlhfl(NHOR)
      real :: sigma(NLEV), sigmah(NLEV), dsigma(NLEV)
      real :: dprcl(NHOR,NLEV), dprscl(NHOR,NLEV), dprll(NHOR,NLEV), dprsll(NHOR,NLEV), dprc(NHOR), dprl(NHOR)
      real :: dprs(NHOR), dgp3d(NHOR,NLEV,20), zdtt_test(NHOR,NLEV)
      real :: dentropy(NHOR,40), dentro(NHOR), dentrot(NHOR,NLEV), dentroq(NHOR,NLEV), dentrop(NHOR)
      real :: dentro3d(NHOR,NLEV,40), denergy(NHOR,40), dener3d(NHOR,NLEV,40)
      integer :: icclev(NHOR,NLEV), icctot(NHOR)
      integer :: nprc=3, nevapprec=1, ndiaggp=0, nentropy=0, nentro3d=0, nenergy=0, nener3d=0
      real :: ga, deltsec2, gamma, ALV, ALS, acpd, ADV, TMELT, ra1, ra2, ra4, rdbrv, gascon
{declarations}      end module rainmod

      subroutine mkdttgp(pdtt)
      use rainmod
      real :: pdtt(NHOR,NLEV)
      pdtt(:,:)=zdtt_test(:,:)
      end subroutine mkdttgp
"""
PLASIM_HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: jhor, jlev, jcall
      real :: zdtdt0(NHOR,NLEP), zdqdt0(NHOR,NLEP)
      read(*,*) ga, deltsec2, gamma, ALV, ALS, acpd, ADV, TMELT, ra1, ra2, ra4, rdbrv, gascon
      read(*,*) sigma, sigmah, dsigma
      read(*,*) nprc
      do jhor=1,NHOR
       read(*,*) dp(jhor), dls(jhor), dshfl(jhor), dlhfl(jhor), dt(jhor,NLEP)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dtdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (zdtt_test(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dw(jhor,jlev), jlev=1,NLEV)
       read(*,*) (du(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dv(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprcl(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprscl(jhor,jlev), jlev=1,NLEV)
      enddo
      dprll(:,:)=0.
      dprsll(:,:)=0.
      zdtdt0(:,:)=dtdt(:,:)
      zdqdt0(:,:)=dqdt(:,:)
      do jcall=1,2
       if(jcall == 1 .and. nprc == 3) call mkmassflux
       if(jcall == 2) call mkrain
       do jhor=1,NHOR
        write(*,'(*(es25.16e3))') (dtdt(jhor,jlev)-zdtdt0(jhor,jlev), jlev=1,NLEV)
        write(*,'(*(es25.16e3))') (dqdt(jhor,jlev)-zdqdt0(jhor,jlev), jlev=1,NLEV)
        write(*,'(*(es25.16e3))') (dprcl(jhor,jlev), jlev=1,NLEV)
        write(*,'(*(es25.16e3))') (dprscl(jhor,jlev), jlev=1,NLEV)
        write(*,'(*(es25.16e3))') dprc(jhor), dprs(jhor)
        write(*,'(*(i4))') icctot(jhor), (icclev(jhor,jlev), jlev=1,NLEV)
       enddo
      enddo
      end program harness
"""


def plasim_harness(folder, rainmod_text, nhor, nlev, routines=True):
    stub = PLASIM_STUB.format(nhor=nhor, nlev=nlev, declarations=mff.declarations())
    files = {'stub.f90': stub}
    if routines:
        files['massflux.f90'] = mff.WRAPPER_CODE + mff.COLUMN_CODE
    else:                                           # the harness's call, never made with nprc 1 or 2
        files['massflux.f90'] = '      subroutine mkmassflux\n      end subroutine mkmassflux\n'
    files['mkrain.f90'] = mkrain_text(rainmod_text)
    files['harness.f90'] = PLASIM_HARNESS
    return compile_in(folder, files)


def plasim_columns(nhor, seed=0):
    """PlaSim's state for nhor columns on the fifteen levels, from the lunar test columns: temperature, vapour, the
    tendencies so far and the dynamics', winds, fluxes (PlaSim's signs, downward), and handed water for mkrain."""
    half = mf.grid_sets()['15']
    sigma = mff_full(half)
    rng = np.random.default_rng(seed)
    cols = []
    for j in range(nhor):
        kw = [dict(t0=300.0, rh_bl=0.78, rh=0.8), dict(t0=302.0, rh_bl=0.8, rh=0.85),
              dict(t0=300.0, rh_bl=0.8, rh=0.6, stable=(2.0e-3, 6000.0)), dict(t0=296.0, rh_bl=0.4, rh=0.4)][j % 4]
        t, q, pf, ph = sounding(G_MOON, 121500.0 + 500.0 * j, n=15, top_m=150000.0, lapse_bl=1.6e-3, lapse=1.0e-3,
                                t_min=200.0, seed=20 + j, c=LUNAR, levels=half, **kw)
        n = t.size
        cols.append(dict(ps=ph[-1], land=float(j % 2), shf=-rng.uniform(10.0, 60.0), lhf=-rng.uniform(50.0, 200.0),
                         tsurf=t[-1] + 2.0, t=t, q=q, dtdt=rng.normal(0.0, 1.0e-5, n), dqdt=rng.normal(0.0, 1.0e-9, n),
                         dqt=rng.normal(0.0, 1.0e-9, n), dtt=rng.normal(0.0, 1.0e-5, n),
                         dw=rng.normal(-0.02, 0.05, n), du=rng.normal(3.0, 2.0, n), dv=rng.normal(0.0, 2.0, n),
                         prcl=np.abs(rng.normal(0.0, 1.0e-7, n)), prscl=np.zeros(n)))
        cols[-1]['prscl'] = cols[-1]['prcl'] * (t < LUNAR['TMELT'])
    return sigma, half, cols


def mff_full(half):
    return 0.5 * (np.r_[0.0, half[:-1]] + half)


def plasim_input(sigma, half, cols, nprc, gamma=0.03, c=LUNAR):
    lines = [num([G_MOON, 3600.0, gamma] + [c[k] for k in CONSTANTS]),
             num(np.r_[sigma, half, np.diff(np.r_[0.0, half])]), str(nprc)]
    for x in cols:
        lines.append(num([x['ps'], x['land'], x['shf'], x['lhf'], x['tsurf']]))
        lines += [num(x[k]) for k in ('t', 'q', 'dtdt', 'dqdt', 'dqt', 'dtt', 'dw', 'du', 'dv', 'prcl', 'prscl')]
    return '\n'.join(lines) + '\n'


def parse_plasim(text, nhor):
    rows = text.strip('\n').split('\n')
    assert len(rows) == 12 * nhor
    out = []
    for call in range(2):
        calls = []
        for j in range(nhor):
            r = rows[(call * nhor + j) * 6:(call * nhor + j) * 6 + 6]
            v = [np.array(r[i].split(), float) for i in range(5)]
            flags = [int(a) for a in r[5].split()]
            calls.append(dict(dtdt=v[0], dqdt=v[1], prcl=v[2], prscl=v[3], prc=v[4][0], prs=v[4][1], icctot=flags[0],
                              icclev=np.array(flags[1:])))
        out.append(calls)
    return out


@FORTRAN
def test_plasims_routine_hands_mkrain_water_that_reaches_the_ground_without_loss(tmp_path):
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    nhor = 8
    sigma, half, cols = plasim_columns(nhor)
    exe = plasim_harness(tmp_path / 'plasim', patched(sources[0].read_text()), nhor, sigma.size)
    for gamma in (0.03, 1.0):
        printed = subprocess.run([str(exe)], input=plasim_input(sigma, half, cols, 3, gamma), capture_output=True,
                                 text=True, check=True).stdout
        after_scheme, after_rain = parse_plasim(printed, nhor)
        convected = 0
        for x, a, b in zip(cols, after_scheme, after_rain):
            m = x['ps'] * np.diff(np.r_[0.0, half]) / G_MOON
            # the routine's own column: mfcolumn as PlaSim's state gives it
            zt = x['t'] + x['dtdt'] * 3600.0
            zq = x['q'] + x['dqdt'] * 3600.0
            p = mf.column(zt, zq, sigma * x['ps'], np.r_[0.0, half * x['ps']], G_MOON, 3600.0, -x['shf'], -x['lhf'],
                          x['dtdt'] + x['dtt'], x['dqdt'] + x['dqt'], x['dw'], np.hypot(x['du'], x['dv']),
                          x['land'] > 0.5, c=LUNAR, melt='plasim', t_phase=x['t'])
            assert np.allclose(a['dtdt'], p['dtdt'], rtol=1e-10, atol=1e-10 * max(np.abs(p['dtdt']).max(), 1e-30))
            assert np.allclose(a['prcl'] * x['ps'] / G_MOON, p['net_rain'] + p['net_snow'], rtol=1e-10,
                               atol=1e-10 * max(np.abs(p['net_rain'] + p['net_snow']).max(), 1e-30))
            if p['type'] is None:
                assert a['icctot'] == 0 and not np.any(a['prcl'])
                continue
            convected += 1
            cloud = np.flatnonzero(a['icclev'])
            assert a['icctot'] == cloud.size >= 1
            assert cloud[-1] + 1 == p['up']['base_half'] and cloud[0] + 1 == p['up']['top_half']
            # mkrain carries the water: what the scheme and mkrain take from the vapour reaches the ground
            vapour = float(np.sum(b['dqdt'] * m))
            ground = (b['prc'] + b['prs']) * 1000.0                       # m/s of water to kg/m2/s
            handed = float(np.abs(a['prcl']).sum()) * x['ps'] / G_MOON
            assert ground >= -1e-12 * handed and b['prs'] >= -1e-12 * handed       # never below nil but round-off
            assert vapour == pytest.approx(-ground, rel=1e-9, abs=1e-15)
            # and the heat: the latent heat of what reaches the ground, rain at ALV and snow at ALS, the scheme's
            # heating over the moist heat capacity of its own vapour, mkrain's over the step's starting vapour's
            cp = LUNAR['acpd']
            heat = float(np.sum(cp * (1.0 + LUNAR['ADV'] * zq) * a['dtdt'] * m)
                         + np.sum(cp * (1.0 + LUNAR['ADV'] * x['q']) * (b['dtdt'] - a['dtdt']) * m))
            snow = b['prs'] * 1000.0
            assert heat == pytest.approx(LUNAR['ALV'] * (ground - snow) + LUNAR['ALS'] * snow, rel=1e-9, abs=1e-12)
        assert convected >= 3


@FORTRAN
def test_with_nprc_1_or_2_mkrain_is_as_it_was(tmp_path):
    sources = installed_sources()
    if not sources:
        pytest.skip('no ExoPlaSim source at hand')
    nhor = 8
    sigma, half, cols = plasim_columns(nhor, seed=3)
    before = apply(sources[0].read_text(), er.SBM_EDITS, er.SBM_MARKER)
    old = plasim_harness(tmp_path / 'old', before, nhor, sigma.size, routines=False)
    new = plasim_harness(tmp_path / 'new', patched(sources[0].read_text()), nhor, sigma.size)
    for nprc in (1, 2):
        for gamma in (0.01, 0.5):
            text = plasim_input(sigma, half, cols, nprc, gamma)
            a = subprocess.run([str(old)], input=text, capture_output=True, text=True, check=True).stdout
            b = subprocess.run([str(new)], input=text, capture_output=True, text=True, check=True).stdout
            assert a == b                                                 # every printed digit, so every bit


def test_the_routines_defaults_are_the_python_copys_parameters():
    values = mff.fortran_values()
    declared = dict(re.findall(r'(?:real|integer) :: (\w+) = (\S+)', mff.declarations()))
    assert set(declared) == set(values)
    for name, value in values.items():
        assert float(declared[name]) == value
    assert values['tstar_mf'] == 0.0 and values['nelev_mf'] == 1
    assert values['zsrc_mf'] * 9.80665 / G_MOON == pytest.approx(60.38, abs=0.01)    # 10 m on Earth, 60 m here
    for key, wrong in (('rh_dependence', True), ('delta_form', 'ifs'), ('source_dp', 5000.0),
                       ('excess_scale', 'mixed_layer')):
        with pytest.raises(ValueError):
            mff.fortran_values(dict(mf.PARAMETERS, **{key: wrong}))


def test_the_runner_records_and_writes_the_mass_flux_scheme_only_for_runs_that_set_it():
    c = er.convection_setting(scheme='massflux')
    assert c == dict(scheme='massflux', nprc=3, parameters=mff.fortran_values())
    entries = er.convection_namelist(c)
    assert entries[0] == ('NPRC', '3') and len(entries) == 1 + len(mff.PARAMETER_NAMES)
    assert dict(entries)['ZSRC_MF'] == '10.0' and dict(entries)['NELEV_MF'] == '1'
    for name, text in entries[1:]:                  # Fortran reads each back to the routine's default
        assert float(text) == mff.fortran_values()[name.lower()]
    assert er.convection_setting(scheme='massflux', rh=None, shallow_cap=None) == c
    with pytest.raises(ValueError):
        er.convection_setting(scheme='massflux', rh=0.7)
    with pytest.raises(ValueError):
        er.convection_setting(scheme='kuo')
    assert er.uses_massflux(dict(scheme='massflux')) and not er.uses_massflux(None)
    assert not er.uses_massflux(dict(rh=0.7, deep_hpa=200.0)) and not er.uses_depth_patch(dict(scheme='massflux'))
    assert er.convection_setting() == dict(scheme='sbm', nprc=2, tau_s=7200.0, rh=0.7, eps=0.0, cin_j_kg=0.0)
    # the command line refuses the Betts-Miller scheme's options with the mass-flux scheme (main sets the runner's
    # signal handlers before it checks, so the test puts the process's own back)
    saved = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        for argv in (['A28R', '--convection', 'massflux', '--sbm-rh', '0.7'],
                     ['A28R', '--convection', 'massflux', '--sbm-shallow-cap']):
            with pytest.raises(SystemExit):
                er.main(argv)
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)
