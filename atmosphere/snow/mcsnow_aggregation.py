"""Snowflakes from measured crystals with McSnow's aggregation: an Earth check, then the lunar column.

Plan step 4.  McSnow (Welss et al. 2024; Brdar and Seifert 2018), pinned at the
revision every snow calculation uses, follows superparticles that settle and
collide in a one-dimensional column.  Its collision kernel is the geometric
encounter area at the difference of fall speeds with Boehm's collision
efficiency, and its sticking efficiency is one of its published choices.  Here
vapour growth is switched off by holding the whole column at ice saturation, so
the crystals injected near the cloud top are the only input: on Earth they are
crystals with Takahashi et al.'s (1991) measured -14.4 C mass-span law, on the
Moon the frozen component's crystals from plan step 3.  Aggregates take
McSnow's established mixed-aggregate geometry (Mitchell 1996).

Three changes to the pinned source, each declared and patch-recorded:
* the monomer mass and area laws read from the namelist (defaults equal the
  source's plate constants, so default runs are unchanged);
* a fall-speed multiplier on vterm_snow (default 1), which carries the crystal
  model's fall-speed band into settling and collision speeds alike;
* for the Moon, lunar gravity, the freezing-level pressure and the 301 K
  ladder's lapse rate in the idealized atmosphere.
The random seed is fixed, as in the existing author-case wrapper.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from shared.constants import MOON_SURFACE_GRAVITY
from atmosphere.snow import welss_aggregation_build as builder
from atmosphere.snow import welss_aggregation as author
from atmosphere.snow import flakes

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = ROOT / 'research/runs/snow_agent/mcsnow_agg'
SOURCE = (Path('/media/projectspace/terluna-research/atmosphere-snow/tmp-copy-2026-10-24') / 'tmp'
          / 'lunar_snow_earth_population' / 'mcsnow_source')
PRODUCT = HERE / 'results/mcsnow_aggregation.json'
EARTH_GRAVITY = 9.81            # the pinned source's value
WATER_DENSITY = 1000.
# Takahashi et al. (1991) crystals within 0.5 C of -14.4 C (21 crystals, 0.27-3.77 mm): least-squares
# fit of ln mass on ln span (data_sized_check crystals); McSnow's own sector-branched crystal area law
# (Mitchell 1996, P1b, quoted in the pinned mo_mass2diam.f90) converted to SI.
EARTH_MONOMER = dict(alf=0.01237, bet=2.207, sig=0.55 * 10 ** (2 * 1.97 - 4), gam=1.97)
# McSnow's temperature-only sticking options, as its mo_colleffi.f90 code defines them (the namelist comment
# calls option 3 an upper limit; the code halves option 2, the lower range of Connolly's Fig. 14).
STICKING = {1: 'Pruppacher and Klett step law (Mitchell 1988)', 2: 'Connolly et al. (2012), Fig. 14, peak 0.64 at -15 C',
            3: 'half of option 2, the lower range of Connolly et al. (2012)', 4: 'Karrer et al. (2021)'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lunar_column():
    column, binding = flakes.moon_column()
    return dict(ps=float(column.p_freeze), tlr=float(column.lapse), depth_15c=float(column.height_of(258.15) - column.z_freeze))


# --- source patches ----------------------------------------------------------------

def patched_sources(world):
    src = SOURCE / 'src'
    out = {}
    nml = (src / 'mo_sp_nml.f90').read_text()
    nml = _replace_once(nml, '    sp_kern_sig       ! added velocity variance\n  LOGICAL     ::   &',
                        '    sp_kern_sig,      &  ! added velocity variance\n'
                        '    band_speed,       &  ! [1] fall-speed multiplier (Terluna error band)\n'
                        '    mono_alf, mono_bet, mono_sig, mono_gam  ! monomer m-D and A-D laws (Terluna)\n'
                        '  LOGICAL     ::   &')
    nml = _replace_once(nml, 'coll_kern_type, sp_kern_sig, &\n            prolate_intp\n',
                        'coll_kern_type, sp_kern_sig, &\n            prolate_intp, band_speed, mono_alf, mono_bet, mono_sig, mono_gam\n')
    nml = _replace_once(nml, '  sp_frag, coll_kern_type, sp_kern_sig, prolate_intp\n',
                        '  sp_frag, coll_kern_type, sp_kern_sig, prolate_intp, band_speed, mono_alf, mono_bet, mono_sig, mono_gam\n')
    nml = _replace_once(nml, '  sp_kern_sig    = 0.0_wp\n',
                        '  sp_kern_sig    = 0.0_wp\n  band_speed     = 1.0_wp\n'
                        '  mono_alf       = 0.7453_wp\n  mono_bet       = 2.47_wp\n'
                        '  mono_sig       = 0.6495_wp\n  mono_gam       = 2.00_wp\n')
    out['mo_sp_nml.f90'] = nml
    vel = (src / 'mo_velocity.f90').read_text()
    vel = _replace_once(vel, '  USE mo_sp_nml,        ONLY: sp_vterm_id, rm, Vm, sp_habit, prolate_intp\n',
                        '  USE mo_sp_nml,        ONLY: sp_vterm_id, rm, Vm, sp_habit, prolate_intp, band_speed\n')
    vel = _replace_once(vel, '    ENDIF\n  ENDIF\nEND FUNCTION vterm_snow',
                        '    ENDIF\n  ENDIF\n  vt = vt * band_speed\nEND FUNCTION vterm_snow')
    out['mo_velocity.f90'] = vel
    m2d = (src / 'mo_mass2diam.f90').read_text()
    m2d = _replace_once(m2d, '  USE mo_sp_nml,        ONLY: snowgeo_type, agggeo_type, sp_vterm_id, rm, Vm, sp_habit\n',
                        '  USE mo_sp_nml,        ONLY: snowgeo_type, agggeo_type, sp_vterm_id, rm, Vm, sp_habit, &\n'
                        '                              mono_alf, mono_bet, mono_sig, mono_gam\n')
    m2d = _replace_once(m2d, '      dunr = (sp%m_i / mono_alf_plates)**(1./mono_bet_plates)\n'
                             '      aunr = mono_sig_plates*dunr**mono_gam_plates\n',
                        '      dunr = (sp%m_i / mono_alf)**(1./mono_bet)\n'
                        '      aunr = mono_sig*dunr**mono_gam\n')
    out['mo_mass2diam.f90'] = m2d
    if world == 'moon':
        col = lunar_column()
        atm = (src / 'mo_atmo_types.f90').read_text()
        atm = _replace_once(atm, '    ps    = 101325._wp,     &', f'    ps    = {col["ps"]:.10g}_wp,     &')
        atm = _replace_once(atm, '    grav  = 9.81_wp,        &', f'    grav  = {MOON_SURFACE_GRAVITY:.17g}_wp,        &')
        atm = _replace_once(atm, '    tlr   = 0.0062_wp,      &', f'    tlr   = {col["tlr"]:.10g}_wp,      &')
        out['mo_atmo_types.f90'] = atm
    return out


def _replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Expected one occurrence of: ' + old[:70])
    return text.replace(old, new)


def build(world):
    path = CACHE / f'build_{world}'
    if path.exists():
        return verify_build(world)
    manifest = json.loads(builder.MANIFEST.read_text())
    for item in manifest['files']:
        builder.verify((SOURCE / item['path']).read_bytes(), item)
    path.mkdir(parents=True)
    tree = path / 'source'
    shutil.copytree(SOURCE / 'src', tree / 'src')
    diff = ''
    for name, body in patched_sources(world).items():
        original = (SOURCE / 'src' / name).read_text()
        (tree / 'src' / name).write_text(body)
        diff += ''.join(difflib.unified_diff(original.splitlines(True), body.splitlines(True),
                                             fromfile='pinned/' + name, tofile=f'{world}/' + name))
    main = (SOURCE / 'src/main.f90').read_text()
    (path / 'main.f90').write_text(builder.deterministic_main(main))
    (path / 'aggregation.patch').write_text(diff)
    compiler = str(Path(shutil.which('gfortran')).resolve())
    files = sorted(p for p in (tree / 'src').rglob('*') if p.suffix.lower() == '.f90' and p.name != 'main.f90') + [path / 'main.f90']
    order = builder.compilation_order(files, compiler, tree / 'src')
    objects, commands = [], []
    with (path / 'build.log').open('w') as log:
        for i, p in enumerate(order):
            obj = path / f'{i:02d}_{p.stem}.o'
            objects.append(str(obj))
            cmd = [compiler, *builder.FLAGS, '-I' + str(tree / 'src'), '-I' + str(path), '-J' + str(path), '-c', str(p), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=path, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *builder.FLAGS, *objects, '-o', str(path / 'MCSNOW')]
        commands.append(cmd)
        subprocess.run(cmd, cwd=path, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.snow-mcsnow-aggregation-build/1', world=world,
                  revision=manifest['revision'], source_dir=str(SOURCE),
                  source_manifest_sha256=sha(builder.MANIFEST), patch_sha256=sha(path / 'aggregation.patch'),
                  compiler=subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
                  compiler_sha256=sha(compiler), flags=builder.FLAGS, commands=commands,
                  executable_sha256=sha(path / 'MCSNOW'))
    (path / 'build.json').write_text(json.dumps(record, indent=1) + '\n')
    return record


def verify_build(world):
    path = CACHE / f'build_{world}'
    record = json.loads((path / 'build.json').read_text())
    if sha(path / 'MCSNOW') != record['executable_sha256']:
        raise ValueError('Executable changed')
    for name, body in patched_sources(world).items():
        if (path / 'source/src' / name).read_text() != body:
            raise ValueError('Patched source differs: ' + name)
    return record


# --- runs ----------------------------------------------------------------------------

def configuration(world, *, monomer_mass, water_mm_h, monomer=None, stick=1, band_speed=1., depth=None,
                  source_fraction=.1, duration=None, dz=20., dtc=5., xi=None, nu_gam=10., kern_sig=0., max_sps=200000):
    """Aggregation-only column: ice saturation everywhere, crystals injected in the top fraction."""
    monomer = monomer or EARTH_MONOMER
    if depth is None:
        depth = 15. / .0062 if world == 'earth' else lunar_column()['depth_15c']
    nz = int(round(depth / dz))
    depth = nz * dz
    belt = source_fraction * depth
    number_flux = water_mm_h / 3600. / monomer_mass          # real crystals per m2 per s
    if xi is None:
        # About 20 new superparticles per collision step in the 5 m2 column.
        xi = max(1, int(round(number_flux * 5. * dtc / 20.)))
    # McSnow floors the superparticles created per step, so make that count an exact integer;
    # the realised flux then differs from the target only through the integer multiplicity.
    per_step = max(1, int(round(number_flux * 5. * dtc / xi)))
    nclrate = (per_step + 1e-9) / (dtc * (belt - 1.) * 5.)    # superparticles per m3 per s in the belt
    if duration is None:
        duration = 4 * 3600. if world == 'earth' else 72 * 3600.
    duration = math.ceil(duration / 600.) * 600.
    return dict(general_nml=dict(msg_level=6, out_type=1, tcase_id=111, ini_type=2,
        snowgeo_type=2, agggeo_type=2, rime_type=0, time_end=float(duration), dt_1dprof=600., dt_2dprof=7200.,
        interpol_atmos=True, bnd_type=0, nucl_type=1, melt_type=1, shedding_type=1, diffusion_type=1,
        hydrobreakup_type=1, collbreakup_type=3, multiplication_type=0, mergeRatioLimit=0.,
        mergeNoTresh=1000000., has_asol=False, perf_test=False),
        sp_nml=dict(iwc=1e-5, lwc0=0., nrp0=10000, xi0=int(xi), maxTotalNoSps=int(max_sps), dtc=float(dtc), idtgm=1,
        sigma0=.62, nclrate=float(nclrate), nclmass=float(monomer_mass), rm=1e-5, sp_vterm_id=3, agg_hist=-1,
        sp_habit=0, coll_kern_type=0, sp_kernel_id=2, prolate_intp=True, sp_stick_effi=int(stick), sp_frag=0,
        sp_kern_sig=float(kern_sig), mu_gam=1., nu_gam=float(nu_gam), band_speed=float(band_speed),
        mono_alf=float(monomer['alf']), mono_bet=float(monomer['bet']),
        mono_sig=float(monomer['sig']), mono_gam=float(monomer['gam'])),
        grid_nml=dict(dom_top=float(depth), dom_bottom=0., box_area=5., nz=nz),
        # The belt ends a metre below the top: a crystal placed exactly at dom_top indexes past the grid.
        # A zero-depth riming zone (h2 = h1) crashes McSnow; a 1 m zone at the bottom, where water and ice
        # saturation coincide at 0 C, with no liquid water and riming off, leaves the column at ice saturation.
        atmo_nml=dict(atmo_type=1, IGF_id=1, ssat=0., h0=0., h1=0., h2=1., nh1=float(depth - belt),
                      nh2=float(depth - 1.), area_factor=1.1, atmofile='ecmwf_profile.txt'))


def run(name, world, config, seed=20261024):
    identity = verify_build(world)
    path = CACHE / 'runs' / name
    if (path / 'run.json').exists():
        record = json.loads((path / 'run.json').read_text())
        if record['configuration'] != config:
            raise ValueError('Existing run has different inputs: ' + name)
        return record
    path.mkdir(parents=True, exist_ok=True)
    (path / 'input.nml').write_text(author.namelist(config))
    for item in json.loads(builder.MANIFEST.read_text())['files']:
        if item['path'].startswith('input/'):
            shutil.copyfile(SOURCE / item['path'], path / Path(item['path']).name)
    start = time.monotonic()
    with (path / 'run.log').open('w') as log:
        proc = subprocess.run([str(CACHE / f'build_{world}/MCSNOW'), str(seed)], cwd=path, stdout=log,
                              stderr=subprocess.STDOUT)
    sp, grid, atmo = config['sp_nml'], config['grid_nml'], config['atmo_nml']
    per_step = math.floor(sp['nclrate'] * sp['dtc'] * (atmo['nh2'] - atmo['nh1']) * grid['box_area'])
    number_flux = per_step * sp['xi0'] / (sp['dtc'] * grid['box_area'])
    record = dict(name=name, world=world, configuration=config, seed=seed, returncode=proc.returncode,
                  wall_s=time.monotonic() - start, executable_sha256=identity['executable_sha256'],
                  realised_number_flux_m2_s=number_flux,
                  realised_water_mm_h=number_flux * sp['nclmass'] * 3600.)
    (path / 'run.json').write_text(json.dumps(record, indent=1) + '\n')
    return record


# Observed snowflake size distributions in melted diameter D: N(D) = N0 exp(-Lambda D), R in mm/h of water.
# Gunn and Marshall (1958): Lambda = 25.5 R^-0.48 cm^-1, fitted above D = 1 mm.
# Sekhon and Srivastava (1970): Lambda = 22.9 R^-0.45 cm^-1.  Median-volume melted diameter D0 = 3.67/Lambda.
OBSERVED = {'gunn_marshall_1958': (2.55, -.48), 'sekhon_srivastava_1970': (2.29, -.45)}   # mm^-1


def observed_d0_mm(water_mm_h):
    return {k: 3.67 / (a * water_mm_h ** b) for k, (a, b) in OBSERVED.items()}


def weighted_quantile(values, weights, q):
    order = np.argsort(values)
    v, w = values[order], weights[order]
    c = np.cumsum(w)
    return float(np.interp(q * c[-1], c, v))


def arriving_snow(name, layer_m=300.):
    """Snow in the lowest layer of the final snapshot, weighted by number and mass flux."""
    path = CACHE / 'runs' / name / 'mass2fr.dat'
    rows = []
    with path.open() as f:
        header = [h.strip() for h in f.readline().split(',')]
        for line in f:
            if ',' in line:
                rows.append([float(x) for x in line.split(',')])
    a = np.asarray(rows)
    col = {h: i for i, h in enumerate(header)}
    column = dict(column_superparticles=int(len(a)), column_largest_mm=float(a[:, col['dia']].max() * 1e3) if len(a) else None,
                  column_superparticles_over_10cm=int((a[:, col['dia']] > .1).sum()) if len(a) else 0)
    low = a[a[:, col['sHeight']] < layer_m] if len(a) else a
    if len(low) < 5:
        return dict(superparticles=int(len(low)), layer_m=layer_m, **column)
    m, v, d, xi, mono = (low[:, col[k]] for k in ('mTot', 'vel', 'dia', 'sMult', 'sNmono'))
    run = json.loads((CACHE / 'runs' / name / 'run.json').read_text())
    area = run['configuration']['grid_nml']['box_area']
    melted = (6 * m / (math.pi * WATER_DENSITY)) ** (1 / 3)
    wn, wm = xi * v, xi * m * v
    flux = float(wm.sum() / (layer_m * area))                     # kg m-2 s-1
    return dict(superparticles=int(len(low)), layer_m=layer_m, **column, arriving_water_mm_h=flux * 3600.,
                realised_injection_mm_h=run['realised_water_mm_h'],
                melted_d0_mm=weighted_quantile(melted, wm, .5) * 1e3,
                max_dimension_d0_mm=weighted_quantile(d, wm, .5) * 1e3,
                max_dimension_d90_mm=weighted_quantile(d, wm, .9) * 1e3,
                max_dimension_d99_mm=weighted_quantile(d, wm, .99) * 1e3,
                largest_mm=float(d.max() * 1e3),
                mean_monomers_by_mass=float((wm * mono).sum() / wm.sum()),
                aggregate_mass_fraction=float(wm[mono > 1.01].sum() / wm.sum()),
                number_flux_m2_s=float(wn.sum() / (layer_m * area)))


EARTH_MONOMER_MASS = 2.8e-8   # kg: Takahashi's fitted -14.4 C history at 20 minutes, 2.75 mm across
EARTH_RATES = (.3, 1., 3., 9.)  # mm/h; 9 mm/h is the CASP II storm of Lawson et al. (1993, 1998)


def earth_matrix():
    rows = []
    for rate in EARTH_RATES:
        for stick in STICKING:
            name = f'earth_r{rate:g}_s{stick}'
            cfg = configuration('earth', monomer_mass=EARTH_MONOMER_MASS, water_mm_h=rate, stick=stick)
            r = run(name, 'earth', cfg)
            if r['returncode'] != 0:
                raise RuntimeError('McSnow failed: ' + name)
            rows.append(dict(name=name, water_mm_h=rate, sticking=stick, wall_s=r['wall_s'],
                             observed_melted_d0_mm=observed_d0_mm(rate), **arriving_snow(name)))
    for seed in (1, 2):
        name = f'earth_r1_s2_seed{seed}'
        cfg = configuration('earth', monomer_mass=EARTH_MONOMER_MASS, water_mm_h=1., stick=2)
        r = run(name, 'earth', cfg, seed=20261024 + seed)
        if r['returncode'] != 0:
            raise RuntimeError('McSnow failed: ' + name)
        rows.append(dict(name=name, water_mm_h=1., sticking=2, seed_offset=seed, wall_s=r['wall_s'],
                         observed_melted_d0_mm=observed_d0_mm(1.), **arriving_snow(name)))
    (CACHE / 'earth_matrix.json').write_text(json.dumps(rows, indent=1) + '\n')
    return rows


# The frozen component's crystals arriving at the lunar freezing level (plan step 3, 29 baseline crystals
# of 0.39-52 mm): least-squares fits of ln mass and ln projected area on ln span.
LUNAR_MONOMER = dict(alf=0.01049, bet=2.0346, sig=0.35225, gam=1.9354)
LUNAR_REFERENCE_FLUX = 10232.4266     # crystals m-2 s-1, the earlier lunar runs' incoming flux
LUNAR_CASES = ((3., 1.), (3., .1), (3., .01), (.3, 1.), (.3, .1), (1., 1.))   # (mm/h, fraction of reference flux)


def lunar_config(water, flux_factor, *, stick=2, band_speed=1., size_factor=1., source_fraction=.1, sp_per_step=2.):
    """size_factor scales the monomer span at fixed mass (the step-3 band's crystal-size edges)."""
    monomer = dict(LUNAR_MONOMER)
    monomer['alf'] = LUNAR_MONOMER['alf'] / size_factor ** LUNAR_MONOMER['bet']
    mass = water / 3600. / (LUNAR_REFERENCE_FLUX * flux_factor)
    # About two new superparticles per 10 s step keep some 30 per 20 m cell at lunar residence times.
    number_flux = water / 3600. / mass
    xi = max(1, int(round(number_flux * 5. * 10. / sp_per_step)))
    return configuration('moon', monomer_mass=mass, water_mm_h=water, monomer=monomer, stick=stick,
                         band_speed=band_speed, source_fraction=source_fraction, duration=72 * 3600., dtc=10., xi=xi)


def lunar_plan(band=None):
    band = band or json.loads((HERE / 'results/earth_error_band.json').read_text())['band']
    sizes = json.loads((HERE / 'results/lunar_size_range.json').read_text())
    ref = next(r for r in sizes['grid'] if r['water_mm_h'] == 3. and r['flux_factor'] == 1.)
    small = ref['combined_small']['diameter_m'] / ref['baseline']['diameter_m']
    large = ref['combined_large']['diameter_m'] / ref['baseline']['diameter_m']
    plan = []
    for water, factor in LUNAR_CASES:
        plan.append((f'moon_q{water:g}_f{factor:g}_s2', dict(water=water, flux_factor=factor)))
    for water, factor in ((3., 1.), (3., .01)):
        base = f'moon_q{water:g}_f{factor:g}'
        plan += [(base + '_speed_low', dict(water=water, flux_factor=factor, band_speed=band['speed']['low']['factor'])),
                 (base + '_speed_high', dict(water=water, flux_factor=factor, band_speed=band['speed']['high']['factor'])),
                 (base + '_size_small', dict(water=water, flux_factor=factor, size_factor=small)),
                 (base + '_size_large', dict(water=water, flux_factor=factor, size_factor=large))]
    for stick in (1, 3, 4):
        plan.append((f'moon_q3_f1_s{stick}', dict(water=3., flux_factor=1., stick=stick)))
    plan.append(('moon_q3_f1_s2_uniform', dict(water=3., flux_factor=1., source_fraction=1.)))
    return plan


def lunar_matrix(part=None, parts=2):
    """Run the lunar plan, or every parts-th entry starting at part, for parallel halves."""
    plan = lunar_plan()
    if part is not None:
        plan = plan[part::parts]
    rows = []
    for name, options in plan:
        cfg = lunar_config(**options)
        r = run(name, 'moon', cfg)
        if r['returncode'] != 0:
            raise RuntimeError('McSnow failed: ' + name)
        rows.append(dict(name=name, options=options, wall_s=r['wall_s'], profile=column_profile(name),
                         **arriving_snow(name, layer_m=2000.)))
        label = 'all' if part is None else str(part)
        (CACHE / f'lunar_matrix_{label}.json').write_text(json.dumps(rows, indent=1) + '\n')
    return rows


def column_profile(name, bin_m=1000., thresholds_m=(.10, .14)):
    """Mass-flux-weighted sizes by height in the final snapshot, and the depth below the top where the
    90th-percentile span first exceeds each threshold (None if it never does)."""
    rows = []
    with (CACHE / 'runs' / name / 'mass2fr.dat').open() as f:
        header = [h.strip() for h in f.readline().split(',')]
        for line in f:
            if ',' in line:
                rows.append([float(x) for x in line.split(',')])
    a = np.asarray(rows)
    col = {h: i for i, h in enumerate(header)}
    top = json.loads((CACHE / 'runs' / name / 'run.json').read_text())['configuration']['grid_nml']['dom_top']
    z, d = a[:, col['sHeight']], a[:, col['dia']]
    w = a[:, col['sMult']] * a[:, col['mTot']] * a[:, col['vel']]
    bins = []
    for upper in np.arange(top, 0, -bin_m):
        s = (z < upper) & (z >= upper - bin_m)
        if s.sum() >= 3:
            bins.append(dict(depth_below_top_km=float((top - upper + bin_m / 2) / 1000), superparticles=int(s.sum()),
                             d50_mm=weighted_quantile(d[s], w[s], .5) * 1e3, d90_mm=weighted_quantile(d[s], w[s], .9) * 1e3))
    onset = {f'{t * 100:g}_cm': next((b['depth_below_top_km'] for b in bins if b['d90_mm'] > t * 1e3), None)
             for t in thresholds_m}
    return dict(bins=bins, onset_depth_km=onset, column_depth_km=top / 1000)


def assemble():
    """Compact product: Earth check, lunar sweep, superparticle refinement and build parity."""
    earth = json.loads((CACHE / 'earth_matrix.json').read_text())
    lunar = []
    for part in ('0', '1'):
        path = CACHE / f'lunar_matrix_{part}.json'
        if path.exists():
            lunar += json.loads(path.read_text())
    refinement = []
    for name, per_step in (('moon_q3_f1_s2', 2), ('moon_q3_f1_s2_sp10', 10), ('moon_q3_f1_s2_sp50', 50)):
        prof = column_profile(name)
        refinement.append(dict(name=name, superparticles_per_step=per_step, profile=prof,
                               **{k: v for k, v in arriving_snow(name, layer_m=2000.).items()
                                  if k.startswith('column_')}))
    parity = json.loads((CACHE / 'parity.json').read_text())
    builds = {w: verify_build(w) for w in ('earth', 'moon')}
    product = dict(
        schema='terluna.atmosphere.snow-mcsnow-aggregation/1',
        evidence=('McSnow aggregation-only columns: an empirical comparison of Earth snow from measured crystals with '
                  'observed snowflake size distributions, then numerical lunar columns from the frozen component\'s '
                  'crystals; no breakup, vapour growth or coupled storm.'),
        sticking=STICKING, earth_monomer=dict(EARTH_MONOMER, mass_kg=EARTH_MONOMER_MASS), lunar_monomer=LUNAR_MONOMER,
        observed_relations={k: dict(lambda_mm=a, exponent=b) for k, (a, b) in OBSERVED.items()},
        earth=earth, lunar=lunar, superparticle_refinement=refinement, author_parity=parity,
        builds={w: dict(executable_sha256=b['executable_sha256'], patch_sha256=b['patch_sha256'],
                        revision=b['revision']) for w, b in builds.items()},
        producer_sha256=sha(Path(__file__)), raw_directory=str(CACHE.relative_to(ROOT)))
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def author_parity(duration=3600., seed=20261021):
    """The authors' aggregation case on the unpatched rebuild and on the patched Earth build."""
    config = author.configuration(True, duration=duration)
    out = {}
    for label, exe in (('reference', CACHE / 'author_reference_build/MCSNOW'), ('patched', CACHE / 'build_earth/MCSNOW')):
        path = CACHE / 'parity' / label
        path.mkdir(parents=True, exist_ok=True)
        (path / 'input.nml').write_text(author.namelist(config))
        for item in json.loads(builder.MANIFEST.read_text())['files']:
            if item['path'].startswith('input/'):
                shutil.copyfile(SOURCE / item['path'], path / Path(item['path']).name)
        start = time.monotonic()
        with (path / 'run.log').open('w') as log:
            proc = subprocess.run([str(exe), str(seed)], cwd=path, stdout=log, stderr=subprocess.STDOUT)
        out[label] = dict(returncode=proc.returncode, wall_s=time.monotonic() - start,
                          files={p.name: sha(p) for p in sorted(path.glob('*.dat'))})
    shared = sorted(set(out['reference']['files']) & set(out['patched']['files']))
    differ = [n for n in shared if out['reference']['files'][n] != out['patched']['files'][n]]
    return dict(duration_s=duration, seed=seed, compared_files=len(shared), differing_files=differ,
                identical=not differ and out['reference']['returncode'] == out['patched']['returncode'] == 0,
                wall_s={k: v['wall_s'] for k, v in out.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'parity', 'earth', 'moon', 'assemble'))
    parser.add_argument('--part', type=int, default=None)
    args = parser.parse_args()
    if args.action == 'build':
        for world in ('earth', 'moon'):
            print(world, build(world)['executable_sha256'])
    elif args.action == 'parity':
        print(json.dumps(author_parity(), indent=1))
    elif args.action == 'moon':
        for r in lunar_matrix(args.part):
            print(f"{r['name']:30} onset depth (D90 > 10 / 14 cm): {r['profile']['onset_depth_km']}; column largest"
                  f" {r['column_largest_mm']:.0f} mm; wall {r['wall_s']:.0f} s", flush=True)
    elif args.action == 'assemble':
        product = assemble()
        print(len(product['earth']), 'Earth runs;', len(product['lunar']), 'lunar runs')
    elif args.action == 'earth':
        for r in earth_matrix():
            print(f"{r['name']:22} D0 melted {r['melted_d0_mm']:.2f} mm (GM {r['observed_melted_d0_mm']['gunn_marshall_1958']:.2f},"
                  f" SS {r['observed_melted_d0_mm']['sekhon_srivastava_1970']:.2f}); Dmax D0 {r['max_dimension_d0_mm']:.1f},"
                  f" D99 {r['max_dimension_d99_mm']:.1f}, largest {r['largest_mm']:.1f} mm; arriving {r['arriving_water_mm_h']:.2f} mm/h")


if __name__ == '__main__':
    main()
