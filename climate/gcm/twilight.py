"""The sunlight the GCM leaves out: the twilight beyond each column's horizon and the low Sun's slant paths.

    $PY -m climate.gcm.twilight optics      # the near-infrared line optics, cached
    $PY -m climate.gcm.twilight longwave    # the longwave line optics, cached
    $PY -m climate.gcm.twilight solve       # the spherical sky's multiple scattering, in batches on the drive
    $PY -m climate.gcm.twilight ledger      # ../results/gcm/twilight.json

$PY is a Python with numba, which the sky model's solver needs (terluna-research/venvs/optics).

ExoPlaSim's radiation (plasim/src/radmod.f90) lights a column only while its Sun stands above the horizon, gmu0 > 0,
and stretches the low Sun's absorber paths with Earth's horizon magnification, 35 / sqrt(1 + 1224 mu0^2); CM1's RRTMG
likewise lights only sunlit columns. The Moon's air is tall against its radius, so its twilight is long and bright,
its horizon air mass near 7 and the air above each square metre of ground a widening cone. This module sets the
energy the curved air takes, by solar elevation and height, against flat columns of the same air and against
PlaSim's own clear-sky arithmetic.

The design column's line-by-line optics from 360 nm to 5 um (the sky model's fine visible grid and the radiative-
convective model's near infrared) are grouped into channels that keep the absorbed energy (faithful_channels). The
channels below 1250 nm that are not opaque take the spherical sky model's multiple scattering (illumination/sky);
the rest take the direct beam through the same 100 spherical shells, with the diffuse and reflected light a
pseudo-spherical two-stream column adds to its own beam. The flat and pseudo-spherical columns of the same layers
are the radiative-convective model's two-stream; PlaSim's is plasim_shortwave.py with the design run's calibration.
Each is integrated over the sphere against the sunlight the solid Moon's disk intercepts, the GCM's footing. The
infrared leaving the top is traced through the same layers, flat and as shells (outgoing_infrared), since the
geometry that takes in more sunlight also lets out more infrared. The terms the GCM's and CM1's twilight patches take
split the geometry into the shells' widening above each square metre of ground (the cone) and the curved light and
escape, on both sides.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE.parent / 'results' / 'gcm'
OUT = RESULTS / 'twilight.json'
RUNS = Path('/media/projectspace/terluna-research/gcm-runs/twilight')
SKY_RUNS = RUNS / 'sky'
SCHEMA = 'terluna.climate.gcm-twilight/1'
WORLD = 'moon_1.2atm'
NIR_NU = (2000.0, 1e7 / 830.0)            # cm^-1: 5 um to 830 nm, inside the shortwave line grid
PROCESSES = 4
GCM_LEVELS = 10
PLASIM_SOLAR1 = 0.502127834309164          # PlaSim's own share of the design sunlight below 0.75 um (its MOST_DIAG)
DESIGN_RUN = 'A28_dim5_moon'               # whose progress record holds that share
SUN_CLASSES = (('below_horizon', -90.0, 0.0), ('low_sun_0_10_deg', 0.0, 10.0), ('day_above_10_deg', 10.0, 90.0))
FINE_STEP = 0.025                          # cm^-1, the sky model's finer visible line grid
SOLVED_NM = 1250.0                         # multiple scattering below this wavelength
SKY_EDGES_NM = np.r_[np.arange(360.0, 500.0, 10.0), np.arange(500.0, 830.0, 20.0), 830.0,
                     np.arange(870.0, SOLVED_NM, 40.0), SOLVED_NM]
NIR_EDGES_NM = np.r_[np.sort(1e7 / np.arange(2000.0 + 50.0, 1e7 / SOLVED_NM + 1e-6, 50.0)), 5000.0 + 1e-6]
TOLERANCE = 0.01                           # faithful_channels' budget, a share of each channel's sunlight
OPAQUE = 3.0                               # whole-column absorption optical depth from which the beam alone serves
SKY_QUALITY = 'standard'
SKY_ALBEDO = 0.1                           # the sky model's ground, as in its products
SKY_BATCH = 24
CHUNK = 200
DIRECT_ELEVATIONS = np.radians(np.r_[np.linspace(-40.0, -10.0, 31)[:-1], np.linspace(-10.0, 10.0, 81)[:-1],
                                     np.linspace(10.0, 90.0, 41)])
STANDARD = ROOT / 'research' / 'runs' / 'optical_comfort' / 'spherical' / f'{WORLD}_standard.npz'
CHECK_MU = (1.0, 0.5, 0.2, 0.05)
TERM_FROM_DEG = -40.0                      # the terms' deepest Sun: 30 degrees down the air takes 0.04 W/m2 of 1167


def log(msg, t0=[time.time()]):
    print(f'{time.time() - t0[0]:8.1f} s  {msg}', file=sys.stderr, flush=True)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def near_infrared_cache(col, cfg):
    """The drive file holding the near infrared's optics for a column and line configuration."""
    key = hashlib.sha256(json.dumps([WORLD, cfg.nu_min, cfg.nu_max, cfg.step, cfg.molecules,
                                     [list(map(float, col['z_m']))]], sort_keys=True).encode()).hexdigest()[:16]
    return RUNS / f'{WORLD}_near_infrared_{key}.npz'


def near_infrared_optics():
    """Absorption and Rayleigh optical depths of the design column (100 layers) from 830 nm to 5 um on the
    radiative-convective model's shortwave line grid, with the shielded sunlight's energy, cached on the drive."""
    from atmosphere.radiative_convective import optics, shortwave as sw
    from illumination.surface_light import model as surface
    col, _, _, _ = surface.build(surface.COLUMNS[WORLD])
    cfg = replace(optics.SHORTWAVE, processes=PROCESSES)
    path = near_infrared_cache(col, cfg)
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            return {k: data[k].copy() for k in data.files}
    log(f'line absorption on {len(col["z_m"]) - 1} layers, {cfg.nu_min:g}-{cfg.nu_max:g} cm-1')
    nu, tau_abs = optics.optical_depth(col, cfg)
    keep = (nu >= NIR_NU[0]) & (nu < NIR_NU[1])
    nu, tau_abs = nu[keep], tau_abs[:, keep]
    tau_sca = optics.rayleigh_optical_depth(col, nu)
    shield, _ = surface.sunlight('design')
    energy, _ = sw.solar_spectrum(nu, shield=shield)
    energy = energy * np.gradient(nu)
    out = dict(nu=nu, absorption=tau_abs.astype(np.float32), scattering=tau_sca.astype(np.float32), energy=energy,
               height=np.asarray(col['z_m'], float))
    RUNS.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, **out)
    tmp.replace(path)
    log(f'cached {path.name}: {nu.size} wavenumbers, {energy.sum():.1f} W/m2 of shielded sunlight')
    return out


def longwave_optics():
    """Absorption optical depths of the design column (100 layers) on the radiative-convective model's longwave line
    grid, cached on the drive."""
    from atmosphere.radiative_convective import optics
    from illumination.surface_light import model as surface
    col, _, _, _ = surface.build(surface.COLUMNS[WORLD])
    cfg = replace(optics.LONGWAVE, processes=PROCESSES)
    key = hashlib.sha256(json.dumps(['longwave', WORLD, cfg.nu_min, cfg.nu_max, cfg.step, cfg.molecules,
                                     [list(map(float, col['z_m']))]], sort_keys=True).encode()).hexdigest()[:16]
    path = RUNS / f'{WORLD}_longwave_{key}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            return {k: data[k].copy() for k in data.files}
    log(f'longwave line absorption on {len(col["z_m"]) - 1} layers, {cfg.nu_min:g}-{cfg.nu_max:g} cm-1')
    nu, tau = optics.optical_depth(col, cfg)
    out = dict(nu=nu, absorption=tau.astype(np.float32), height=np.asarray(col['z_m'], float))
    RUNS.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, **out)
    tmp.replace(path)
    log(f'cached {path.name}: {nu.size} wavenumbers')
    return out


LW_NODES = 8                               # Gauss nodes in each panel of rays that reach the ground
LW_SHELL_NODES = 4                         # and in each shell's panel of the rays above it
LW_PANELS_KM = (0.0, 0.1, 1.0, 5.0, 20.0, 60.0, 150.0)   # the ground rays' panels by depth below the ground's limb
LW_CHUNK = 20000


def outgoing_infrared(nu, tau, edges, t_layers, t_surface):
    """The infrared leaving the top of a column of isothermal layers (optical depths bottom-up) over a black ground,
    per square metre of ground: through flat layers, and through the same layers as spherical shells. The shells'
    rays are taken by their distance b from the Moon's centre, (2 pi / R^2) * integral of I(b) b db per square metre
    of ground, in panels that close in on the ground's limb from inside (rays that reach the ground) and, one panel
    to a shell, between the shells' edges outside it (rays that cross the air twice, whose paths kink where their
    tangent point crosses an edge). Each layer's share is kept, so the shells' widening (a shell emits from (r/R)^2
    of the flat layer's volume) and the curved escape split as the sunlight's cone and curved light do. W/m2."""
    from atmosphere.radiative_convective import longwave as lw
    R, top = edges[0], edges[-1]
    dz = np.diff(edges)
    nl = len(dz)
    step = nu[1] - nu[0]
    x, w = np.polynomial.legendre.leggauss(LW_NODES)
    unit = lambda lo, hi: (lo + (hi - lo) * 0.5 * (x + 1), (hi - lo) * 0.5 * w)
    rays = []                                  # (geometry, weight, [(layer, path)], reaches the ground)
    for m, wt in zip(*unit(0.0, 1.0)):
        rays.append(('flat', 2.0 * math.pi * wt * m, [(l, dz[l] / m) for l in range(nl - 1, -1, -1)], True))
    panels = np.array(LW_PANELS_KM) * 1e3
    inside = np.r_[panels, R][::-1]                              # depths R - b below the limb, deepest first
    xs, ws = np.polynomial.legendre.leggauss(LW_SHELL_NODES)
    spans = [(R - a, R - c, x, w) for a, c in zip(inside[:-1], inside[1:])] + \
            [(lo, hi, xs, ws) for lo, hi in zip(edges[:-1], edges[1:])]
    for lo, hi, xn, wn in spans:
        for b, wt in zip(lo + (hi - lo) * 0.5 * (xn + 1), (hi - lo) * 0.5 * wn):
            along = lambda l: (math.sqrt(max(edges[l + 1] ** 2 - b * b, 0.0))
                               - math.sqrt(max(max(edges[l], b) ** 2 - b * b, 0.0)))
            near = [(l, along(l)) for l in range(nl - 1, -1, -1) if edges[l + 1] > b]
            ground = b < R
            rays.append(('sphere', 2.0 * math.pi * wt * b / R ** 2, near if ground else near + near[::-1], ground))
    out = {g: dict(layers=np.zeros(nl), ground=0.0) for g in ('flat', 'sphere')}
    for c0 in range(0, nu.size, LW_CHUNK):
        sl = slice(c0, min(c0 + LW_CHUNK, nu.size))
        k = tau[:, sl] / dz[:, None]
        B = lw.planck(nu[sl][None, :], np.asarray(t_layers)[:, None])
        Bs = lw.planck(nu[sl], t_surface)
        edge_w = np.full(nu[sl].size, step)
        if c0 == 0:
            edge_w[0] *= 0.5
        if sl.stop == nu.size:
            edge_w[-1] *= 0.5
        for geometry, weight, segments, ground in rays:
            seen = np.ones(nu[sl].size)
            for l, ds in segments:
                e = np.exp(-k[l] * ds)
                out[geometry]['layers'][l] += weight * float(((B[l] * (1.0 - e)) * seen) @ edge_w)
                seen *= e
            if ground:
                out[geometry]['ground'] += weight * float((Bs * seen) @ edge_w)
    widening = ((0.5 * (edges[1:] + edges[:-1])) / R) ** 2
    sphere = out['sphere']['layers']                            # per square metre of ground
    shallow = sphere / widening                                 # the same escape from a flat layer's volume
    flat = out['flat']['layers']
    total = lambda layers, g: float(layers.sum() + g)
    return dict(flat=total(flat, out['flat']['ground']), sphere=total(sphere, out['sphere']['ground']),
                sphere_shallow=total(shallow, out['sphere']['ground']),
                ground=dict(flat=out['flat']['ground'], sphere=out['sphere']['ground']),
                layers=dict(flat=flat, sphere=sphere, sphere_shallow=shallow))


def fine_inputs():
    """Line-by-line optics of the design column's 100 layers from 360 nm to 5 um, by rising wavelength: the sky
    model's fine visible grid (illumination/sky/solved_optics, with its trace gases and colour weights) and the near
    infrared's; with the files they come from."""
    from atmosphere.radiative_convective import optics
    from illumination.sky import solved_optics
    from illumination.surface_light import model as surface
    fine, meta = solved_optics.fine_optics(WORLD, STANDARD.parent, step_cm1=FINE_STEP)
    nir = near_infrared_optics()
    order = np.argsort(1e7 / nir['nu'])
    col, _, _, _ = surface.build(surface.COLUMNS[WORLD])
    files = (STANDARD.parent / f'{WORLD}_fine_{meta["key"][:16]}.npz',
             near_infrared_cache(col, replace(optics.SHORTWAVE, processes=PROCESSES)))
    sources = {str(f): digest(f) for f in files}
    return dict(height=fine['height'], wavelength=np.r_[fine['wavelength'], (1e7 / nir['nu'])[order]],
                energy=np.r_[fine['energy'], nir['energy'][order]],
                absorption=np.c_[fine['absorption'], nir['absorption'][:, order]],
                scattering=np.c_[fine['scattering'], nir['scattering'][:, order]],
                xyz=np.r_[fine['xyz'], np.zeros((order.size, 3))]), sources


FAITHFUL_PATHS = (0.03, 0.1, 0.3, 1.0, 3.0, 13.0)    # slant columns, in whole vertical columns


def faithful_channels(absorption, scattering, energy, wavelength, edges_nm, tolerance, smooth_below_nm=0.0,
                      xyz=None):
    """Channels that keep the sunlight's absorbed energy. Within each band (edges in nm) the wavenumbers are sorted
    by whole-column absorption and bisected by energy until each group, carrying its layers' energy-weighted optical
    depths, absorbs a direct beam along slant columns of 0.03 to 13 whole columns within tolerance (a share of the
    group's sunlight) of what its members absorb. An averaged optical depth always absorbs more than its members,
    so the bias of coarse groups runs one way. Bands below smooth_below_nm keep one channel."""
    col = absorption.sum(0)
    m = np.asarray(FAITHFUL_PATHS)
    groups = []
    for lo, hi in zip(edges_nm[:-1], edges_nm[1:]):
        ix = np.flatnonzero((wavelength >= lo) & (wavelength < hi))
        if ix.size == 0:
            continue
        ix = ix[np.argsort(col[ix], kind='stable')]
        if hi <= smooth_below_nm:
            groups.append((ix, lo, hi))
            continue
        stack = [ix]
        while stack:
            g = stack.pop()
            e = energy[g]
            total = e.sum()
            mean = col[g] @ e / total
            exact = (1.0 - np.exp(-m[:, None] * col[g][None, :])) @ e
            if g.size == 1 or np.max(np.abs(total * (1.0 - np.exp(-m * mean)) - exact)) <= tolerance * total:
                groups.append((g, lo, hi))
            else:
                k = int(np.clip(np.searchsorted(np.cumsum(e), 0.5 * total), 1, g.size - 1))
                stack += [g[:k], g[k:]]
    groups.sort(key=lambda x: (x[1], col[x[0]].mean()))
    rows = []
    for g, lo, hi in groups:
        w = energy[g] / energy[g].sum()
        rows.append((energy[g].sum(), absorption[:, g] @ w, scattering[:, g] @ w, float(wavelength[g] @ w), (lo, hi),
                     np.zeros(3) if xyz is None else xyz[g].sum(0)))
    return dict(energy=np.array([r[0] for r in rows]), absorption=np.stack([r[1] for r in rows], 1),
                scattering=np.stack([r[2] for r in rows], 1), wavelength=np.array([r[3] for r in rows]),
                band=np.array([r[4] for r in rows]), xyz=np.array([r[5] for r in rows]))


def channels(fine):
    """The ledger's channels over 360 nm-5 um, each marked for the multiple-scattering solve (below SOLVED_NM and not
    opaque) or the direct beam."""
    parts = [faithful_channels(fine['absorption'], fine['scattering'], fine['energy'], fine['wavelength'], edges,
                               TOLERANCE, smooth, fine['xyz'])
             for edges, smooth in ((SKY_EDGES_NM, 500.0), (NIR_EDGES_NM, 0.0))]
    ch = {k: np.concatenate([p[k] for p in parts], axis=1 if k in ('absorption', 'scattering') else 0) for k in parts[0]}
    ch['solved'] = (ch['wavelength'] < SOLVED_NM) & (ch['absorption'].sum(0) < OPAQUE)
    ch['height'] = fine['height']
    return ch


def sky_batches(ch, solve=False):
    """The solved channels' multiple scattering by the spherical sky model's solver (its standard quality and
    ground), in batches kept on the drive so an interruption loses one; with solve, the missing ones are solved."""
    from illumination.sky import solved_transport as st
    from illumination.surface_light import model as surface
    pick = np.flatnonzero(ch['solved'])
    key = hashlib.sha256(b''.join(np.ascontiguousarray(ch[k][..., pick]).tobytes()
                                  for k in ('absorption', 'scattering', 'energy'))
                         + json.dumps([SKY_QUALITY, SKY_ALBEDO, digest(st.__file__)]).encode()).hexdigest()[:16]
    starts = range(0, pick.size, SKY_BATCH)
    paths = [SKY_RUNS / f'sky_{key}_{b:02d}.npz' for b in range(len(starts))]
    for b, (start, path) in enumerate(zip(starts, paths)):
        if path.exists() or not solve:
            continue
        sel = pick[start:start + SKY_BATCH]
        log(f'sky batch {b + 1} of {len(paths)}: {sel.size} channels, {ch["wavelength"][sel].min():.0f}-'
            f'{ch["wavelength"][sel].max():.0f} nm')
        sub = dict(height=ch['height'], absorption=ch['absorption'][:, sel], scattering=ch['scattering'][:, sel],
                   energy=ch['energy'][sel], xyz=ch['xyz'][sel], wavelength=ch['wavelength'][sel], band=ch['band'][sel])
        SKY_RUNS.mkdir(parents=True, exist_ok=True)
        st.solve(sub, surface.COLUMNS[WORLD].planet.radius_m, SKY_QUALITY, albedo=SKY_ALBEDO, path=path)
    missing = [p.name for p in paths if not p.exists()]
    if missing and not solve:
        raise SystemExit(f'{len(missing)} of {len(paths)} sky batches are not solved yet (task solve): {missing[0]}...')
    return pick, paths


def bracket(grid, x):
    i = int(np.clip(np.searchsorted(grid, x) - 1, 0, len(grid) - 2))
    return i, (x - grid[i]) / (grid[i + 1] - grid[i])


def column_quadrature(edges, nodes):
    """Gauss points in radius splitting at the layers' and the field's nodes, with their layers."""
    gauss, weights = np.polynomial.legendre.leggauss(3)
    splits = np.unique(np.r_[edges, nodes])
    radii = ((splits[1:, None] + splits[:-1, None]) / 2 + np.diff(splits)[:, None] * gauss / 2).ravel()
    w = (np.diff(splits)[:, None] * weights / 2).ravel()
    layer = np.clip(np.searchsorted(edges, radii) - 1, 0, len(edges) - 2)
    return radii, w, layer


def field(path, albedo, elevations):
    """What a solved batch's light leaves, by solar elevation at a ground point: absorbed by the ground and in each
    layer of the cone above it, per unit area of ground, and in the same layers of a flat column under the same light
    (without the cone's widening, (r/R)^2), per channel per unit incident flux. The diffuse field is the retained
    series (moments less the last order, as the sky product evaluates it) and the beam the solver's own table, each
    taken linearly to the elevations."""
    from illumination.sky.solved_transport import solar_table
    with np.load(path, allow_pickle=False) as z:
        s = {k: z[k] for k in ('r', 'a', 'sr', 'sa', 'beam', 'moments', 'last_order', 'edges', 'scattering',
                               'absorption')}
    r, a, edges = s['r'], s['a'], s['edges']
    m = at_elevations(a, np.moveaxis(s['moments'] - s['last_order'], 1, 0), elevations)    # (elev, r, ch, 7)
    direct = solar_table(np.array([r[0]]), elevations, edges, s['scattering'] + s['absorption'])[0]
    ground = (1.0 - albedo) * (m[:, 0, :, 5] + direct[:, :, 1])
    beam = at_elevations(s['sa'], np.moveaxis(s['beam'][..., 0], 1, 0), elevations)          # (elev, sr, ch)
    radii, w, layers = column_quadrature(edges, r)
    air = np.zeros((len(edges) - 1, elevations.size, s['absorption'].shape[1]))
    flat = np.zeros_like(air)
    for rad, wt, lay in zip(radii, w, layers):
        ir, u = bracket(r, rad)
        bi, bu = bracket(s['sr'], rad)
        j = (1 - u) * m[:, ir, :, 0] + u * m[:, ir + 1, :, 0]
        f = (1 - bu) * beam[:, bi] + bu * beam[:, bi + 1]
        taken = wt * s['absorption'][lay][None, :] * (j + f)
        air[lay] += (rad / r[0]) ** 2 * taken
        flat[lay] += taken
    return ground, air, flat


def direct_beams(tau_abs, tau_sca, edges, albedo, elevations):
    """The direct beam through the spherical shells (layer optical depths, bottom-up): absorbed in each layer of the
    cone, in the same layers of a flat column under the same beam, and by the ground, per unit ground area, by solar
    elevation, per channel per unit incident flux."""
    from illumination.sky.solved_transport import solar_table
    dz = np.diff(edges)
    kabs, kext = tau_abs / dz[:, None], (tau_abs + tau_sca) / dz[:, None]
    radii, w, layers = column_quadrature(edges, edges)
    n = tau_abs.shape[1]
    air = np.zeros((len(dz), elevations.size, n)); ground = np.zeros((elevations.size, n))
    flat = np.zeros_like(air)
    for c0 in range(0, n, CHUNK):
        sl = slice(c0, min(c0 + CHUNK, n))
        t = solar_table(np.r_[edges[0], radii], elevations, edges, np.ascontiguousarray(kext[:, sl]))
        ground[:, sl] = (1.0 - albedo) * t[0, :, :, 1]
        for i, (rad, wt, lay) in enumerate(zip(radii, w, layers)):
            taken = wt * kabs[lay, sl][None, :] * t[i + 1, :, :, 0]
            air[lay, :, sl] += (rad / edges[0]) ** 2 * taken
            flat[lay, :, sl] += taken
    return air, ground, flat


def two_stream(tau_abs, tau_sca, albedo, elevation, edges=None):
    """The radiative-convective model's two-stream columns (layer optical depths bottom-up) for the same layers, flat
    or with its pseudo-spherical direct beam when edges are given: absorbed in each layer (bottom-up) and by the
    ground per unit ground area, by solar elevation, per channel per unit incident flux; nothing below the horizon."""
    from atmosphere.radiative_convective import shortwave as sw
    ta, ts = tau_abs[::-1], tau_sca[::-1]
    g = np.zeros_like(ta)
    nl, nc = ta.shape
    air = np.zeros((nl, elevation.size, nc)); ground = np.zeros((elevation.size, nc))
    radius = None if edges is None else edges[::-1]
    for i, el in enumerate(elevation):
        mu0 = math.sin(el)
        if mu0 <= 1e-4:
            continue
        direct, diffuse, up = sw.column_fluxes(ta, ts, g, albedo, mu0, radius)
        net = direct + diffuse - up                              # levels top-down
        air[:, i] = (net[:-1] - net[1:])[::-1]
        ground[i] = net[-1]
    return air, ground


def direct_only(tau_abs, tau_sca, albedo, elevation, edges):
    """two_stream's pseudo-spherical direct beam alone: what it leaves in each layer (bottom-up) and on the ground,
    per unit ground area, by solar elevation, per channel per unit incident flux."""
    from atmosphere.radiative_convective import shortwave as sw
    ext = (tau_abs + tau_sca)[::-1]
    share = np.where(ext > 0, tau_abs[::-1] / np.maximum(ext, 1e-300), 0.0)
    zero = np.zeros_like(ext)
    air = np.zeros((ext.shape[0], elevation.size, ext.shape[1])); ground = np.zeros((elevation.size, ext.shape[1]))
    for i, el in enumerate(elevation):
        mu0 = math.sin(el)
        if mu0 <= 1e-4:
            continue
        direct, _, _ = sw.column_fluxes(ext, zero, zero, 0.0, mu0, edges[::-1])
        air[:, i] = ((direct[:-1] - direct[1:]) * share)[::-1]
        ground[i] = (1.0 - albedo) * direct[-1]
    return air, ground


def column_absorption(tau_abs, tau_sca, energy, albedo, mu0, edges, chunk=30000):
    """The pseudo-spherical two-stream's absorption in each layer (bottom-up) and by the ground, W/m2 at the beam for
    the channels' energy, the Sun at mu0; in chunks of channels."""
    from atmosphere.radiative_convective import shortwave as sw
    air, ground = np.zeros(tau_abs.shape[0]), 0.0
    for s in range(0, tau_abs.shape[1], chunk):
        sl = slice(s, min(s + chunk, tau_abs.shape[1]))
        ta = np.ascontiguousarray(tau_abs[::-1, sl])
        d, f, u = sw.column_fluxes(ta, np.ascontiguousarray(tau_sca[::-1, sl]), np.zeros_like(ta), albedo, mu0,
                                   edges[::-1])
        net = (d + f - u) @ energy[sl]
        air += (net[:-1] - net[1:])[::-1]
        ground += net[-1]
    return air, ground


def plasim_sigma(levels=GCM_LEVELS):
    """PlaSim's default vertical grid (ExoPlaSim's vtype 0), full levels top first."""
    s = np.arange(1, levels + 1) / levels
    half = 0.75 * s + 1.75 * s ** 3 - 1.5 * s ** 4
    return 0.5 * (np.r_[0.0, half[:-1]] + half)


def water_fraction(x, molar_mass):
    """Water vapour's share of the mass of moist air holding mole fraction x of it, the dry air's molar mass given."""
    from atmosphere.radiative_convective import thermodynamics as th
    mw = th.MOLAR_MASS['H2O']
    return x * mw / (x * mw + molar_mass * (1 - x))


def plasim_column(col, molar_mass, albedo):
    """The solved column as the design run's PlaSim holds it: the GCM's levels, surface pressure, gravity and
    radiation settings, with the solved column's temperature and water at each level's share of its pressure."""
    from climate.gcm import exoplasim_run as er, plasim_shortwave as ps
    sigma = plasim_sigma()
    lp = np.log(col['layer_p_pa'][::-1])
    at = lambda v: np.interp(np.log(sigma * col['surface_pressure_pa']), lp, np.asarray(v)[::-1])
    planet = json.loads((er.PRODUCTS / 'moon_gcm_configuration.json').read_text())['planet']
    return dict(sigma=sigma, ps=er.EXPERIMENTS['A28']['pressure_pa'], t=at(col['layer_t_k']),
                q=water_fraction(at(col['layer_x_h2o']), molar_mass), albedo=albedo,
                rcoeff=ps.RCOEFF_PER_RAYSCALE * er.MODEL['rayleigh_scale'], solar1=PLASIM_SOLAR1,
                gravity=planet['gravity_m_s2'])


def plasim(column, elevation):
    """PlaSim's clear sky on the column with the design run's water absorption: absorbed in each of its layers (top
    first) and by the ground, and reflected, per unit ground area per unit incident flux, by solar elevation;
    nothing while the Sun is below the horizon (its dawn, gmu0 > 0)."""
    from climate.gcm import exoplasim_run as er, plasim_shortwave as ps
    layers = np.zeros((len(column['sigma']), elevation.size)); ground = np.zeros(elevation.size)
    reflected = np.zeros(elevation.size)
    for i, el in enumerate(elevation):
        mu = math.sin(el)
        if mu <= 0.0:
            continue
        net = ps.net_flux(column, mu, 'scale', er.MODEL['sw_water_absorption'])
        layers[:, i], ground[i], reflected[i] = mu * (net[:-1] - net[1:]), mu * net[-1], mu * (1.0 - net[0])
    return layers, ground, reflected


def onto_gcm_layers(air, p_edges):
    """Absorption in the shells (bottom-up, first axis) shared among the GCM's layers (top first) by their overlap
    in sigma, each level's share of the column's surface pressure."""
    from climate.gcm import plasim_shortwave as ps
    s = np.asarray(p_edges, float) / p_edges[0]
    half = np.r_[0.0, ps.half_levels(plasim_sigma())]
    lo, hi = s[1:], s[:-1]
    out = np.zeros((half.size - 1,) + air.shape[1:])
    for l in range(half.size - 1):
        overlap = np.clip(np.minimum(hi, half[l + 1]) - np.maximum(lo, half[l]), 0.0, None) / (hi - lo)
        out[l] = np.tensordot(overlap, air, axes=(0, 0))
    return out


def air_mass(edges, density, elevation):
    """Slant column over vertical column from the ground towards a Sun at each elevation, through spherical shells
    of constant density (per metre of height, bottom-up)."""
    b = edges[0] * np.cos(np.atleast_1d(elevation))[:, None]
    path = np.diff(np.sqrt(np.maximum(edges[None, :] ** 2 - b ** 2, 0.0)), axis=1)
    return path @ density / (density @ np.diff(edges))


MAGNIFICATION_ELEVATIONS = np.radians(np.r_[np.linspace(0.0, 10.0, 101), np.linspace(10.5, 90.0, 160)])


def magnification(edges, water):
    """The Moon's slant factor for its water vapour (per metre of height, bottom-up) by solar elevation, and the
    rational form PlaSim's patch takes for it, (mu + a) / (mu^2 + b mu + c): exactly 1 overhead, c = a / m(0) at the
    horizon, a and b fitted in the logarithm."""
    from scipy.optimize import least_squares
    el = MAGNIFICATION_ELEVATIONS
    mu = np.sin(el)
    m = air_mass(edges, water, el)
    form = lambda p: (mu + p[0]) / (mu * mu + p[1] * mu + p[0] / m[0])
    fit = least_squares(lambda p: np.log(form(p)) - np.log(m), [0.18, 0.17], bounds=([1e-4, -0.5], [5.0, 5.0]))
    a, b = (float(x) for x in fit.x)
    error = form(fit.x) / m - 1.0
    return dict(form='(mu + a) / (mu^2 + b mu + c)', a=a, b=b, c=a / float(m[0]), horizon=float(m[0]),
                max_error=float(np.max(np.abs(error))), elevation_deg=list(np.degrees(el[::10])), water=list(m[::10]))


def sphere_share(x, elevation):
    """The share of the solid disk's intercepted power a per-ground-area quantity x(elevation) (per unit incident
    flux) makes over the sphere: 2 * integral of x d(sin a)."""
    return 2.0 * np.trapezoid(x, np.sin(elevation), axis=0)


def at_elevations(grid, x, points):
    """x (elevation first) taken linearly to points inside the grid."""
    i = np.clip(np.searchsorted(grid, points) - 1, 0, len(grid) - 2)
    u = ((points - grid[i]) / (grid[i + 1] - grid[i])).reshape((-1,) + (1,) * (x.ndim - 1))
    return (1 - u) * x[i] + u * x[i + 1]


def onto(grid, src, x):
    """x (elevation first, on the elevations src) taken to grid, nothing outside src's range."""
    inside = ((grid >= src[0]) & (grid <= src[-1])).reshape((-1,) + (1,) * (x.ndim - 1))
    return np.where(inside, at_elevations(src, x, np.clip(grid, src[0], src[-1])), 0.0)


def class_share(x, grid, lo_deg, hi_deg):
    """sphere_share over the solar elevations from lo to hi degrees only."""
    lo, hi = max(math.radians(lo_deg), grid[0]), min(math.radians(hi_deg), grid[-1])
    pts = np.unique(np.r_[lo, grid[(grid > lo) & (grid < hi)], hi])
    return sphere_share(at_elevations(grid, x, pts), pts)


def global_means(ground, layers, grid):
    """Global-mean W/m2 of absorption given per unit ground area (W/m2 at the beam's flux) by solar elevation: the
    ground and the air, each whole and by class of solar elevation, and the air by layer."""
    air = layers.sum(0)
    out = dict(ground=sphere_share(ground, grid) / 4.0, air=sphere_share(air, grid) / 4.0)
    out['by_sun'] = {name: dict(ground=class_share(ground, grid, lo, hi) / 4.0, air=class_share(air, grid, lo, hi) / 4.0)
                     for name, lo, hi in SUN_CLASSES}
    out['air_by_layer'] = list(sphere_share(layers.T, grid) / 4.0)
    return out


def rounded(x, digits=6):
    if isinstance(x, dict):
        return {k: rounded(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [rounded(v, digits) for v in (x.tolist() if isinstance(x, np.ndarray) else x)]
    if isinstance(x, (float, np.floating)):
        return float(f'{float(x):.{digits}g}')
    return x


def reduction_check(fine, ch, edges, p_edges):
    """The channels' pseudo-spherical two-stream absorption against the fine grid's, and the sky product's own
    production channels' (360-830 nm), with the Sun at CHECK_MU: W/m2 at the beam, the air by GCM layer."""
    with np.load(STANDARD, allow_pickle=False) as z:
        dz = np.diff(z['edges'])
        production = (z['absorption'] * dz[:, None], z['scattering'] * dz[:, None], z['energy'], z['wavelength'])
    rows = []
    for lo, hi in ((360.0, 830.0), (830.0, 5000.1)):
        f = (fine['wavelength'] >= lo) & (fine['wavelength'] < hi)
        c = (ch['wavelength'] >= lo) & (ch['wavelength'] < hi)
        for mu in CHECK_MU:
            row = dict(band_nm=[lo, min(hi, 5000.0)], mu0=mu)
            for name, (ta, ts, e) in dict(fine=(fine['absorption'][:, f], fine['scattering'][:, f], fine['energy'][f]),
                                          channels=(ch['absorption'][:, c], ch['scattering'][:, c], ch['energy'][c])).items():
                air, ground = column_absorption(ta, ts, e, SKY_ALBEDO, mu, edges)
                row[name] = dict(air=float(air.sum()) / mu, ground=float(ground) / mu,
                                 air_by_layer=list(onto_gcm_layers(air[:, None], p_edges)[:, 0] / mu))
            if hi < 1000.0:
                air, ground = column_absorption(*production[:3], SKY_ALBEDO, mu, edges)
                row['production'] = dict(air=float(air.sum()) / mu, ground=float(ground) / mu)
            rows.append(row)
    return rows


def ledger() -> dict:
    """The energy the curved air takes by solar elevation, ground and height, against flat columns of the same air
    and PlaSim's clear sky on the design run's column, integrated over the sphere."""
    from atmosphere.radiative_convective import optics, shortwave as sw
    from illumination.sky import solved_optics, solved_transport
    from illumination.surface_light import model as surface
    from climate.gcm import exoplasim_run as er, plasim_shortwave as ps
    col, air_spec, _, _ = surface.build(surface.COLUMNS[WORLD])
    R = surface.COLUMNS[WORLD].planet.radius_m
    fine, sources = fine_inputs()
    ch = channels(fine)
    edges, p_edges, albedo = R + ch['height'], col['p_pa'], SKY_ALBEDO
    pick, paths = sky_batches(ch)
    sky_a, meta = None, []
    for path in paths:
        with np.load(path, allow_pickle=False) as z:
            sky_a = z['a'] if sky_a is None else sky_a
            meta.append(json.loads(str(z['meta'])))
    # the two grids' shared angles differ in their last bits; one of each, so no interval is empty
    grid = np.radians(np.unique(np.round(np.degrees(np.r_[sky_a, DIRECT_ELEVATIONS]), 6)))
    beam = np.flatnonzero(~ch['solved'])
    e = ch['energy']
    log(f'{e.size} channels, {pick.size} solved in {len(paths)} batches, {beam.size} by their beam')
    sphere_g, sphere_a = np.zeros((grid.size, e.size)), np.zeros((len(edges) - 1, grid.size, e.size))
    shallow_a = np.zeros_like(sphere_a)
    start = 0
    for path in paths:
        g, a, fl = field(path, albedo, grid)
        sel = pick[start:start + g.shape[1]]
        sphere_g[:, sel], sphere_a[:, :, sel], shallow_a[:, :, sel] = g, a, fl
        start += g.shape[1]
    log('the beam channels: the curved beam and the two-stream')
    ta, ts = ch['absorption'], ch['scattering']
    ba, bg, bf = direct_beams(ta[:, beam], ts[:, beam], edges, albedo, grid)
    pa, pg = two_stream(ta, ts, albedo, grid, edges)
    da, dg = direct_only(ta[:, beam], ts[:, beam], albedo, grid, edges)
    # The beam channels' sphere: the curved direct beam, with the diffuse and reflected light the pseudo-spherical
    # two-stream adds to its own direct beam.
    sphere_g[:, beam] = bg + pg[:, beam] - dg
    sphere_a[:, :, beam] = ba + pa[:, :, beam] - da
    shallow_a[:, :, beam] = bf + pa[:, :, beam] - da
    beams = dict(sphere=(bg @ e[beam], ba @ e[beam]), pseudo_spherical=(dg @ e[beam], da @ e[beam]))
    del ba, bf, da
    fa, fg = two_stream(ta, ts, albedo, grid)
    per_channel = dict(sphere=(sphere_g, sphere_a), sphere_shallow=(sphere_g, shallow_a), pseudo_spherical=(pg, pa),
                       flat=(fg, fa))
    covered = float(e.sum())
    bands = dict(visible=ch['wavelength'] < 830.0, near_infrared=ch['wavelength'] >= 830.0)
    by_band = {name: {k: global_means(g[:, m] @ e[m], onto_gcm_layers(a[:, :, m] @ e[m], p_edges), grid)
                      for k, (g, a) in per_channel.items()} for name, m in bands.items()}
    by_band['beam_channels_direct_beam'] = {k: global_means(g, onto_gcm_layers(a, p_edges), grid)
                                            for k, (g, a) in beams.items()}
    tables = {k: (g @ e, onto_gcm_layers(a @ e, p_edges)) for k, (g, a) in per_channel.items()}
    shell_delta = shallow_a @ e - pa @ e                        # (shells bottom-up, elevations), W/m2 at the beam
    shell_cone = sphere_a @ e - shallow_a @ e
    del per_channel, sphere_a, shallow_a, pa, fa
    log('PlaSim on the design run\'s column')
    column = plasim_column(col, air_spec.molar_mass, albedo)
    p_layers, p_ground, p_reflected = plasim(column, grid)
    tables['plasim'] = (covered * p_ground, covered * p_layers)
    whole = {k: global_means(g, a, grid) for k, (g, a) in tables.items()}

    def difference(x, y):
        out = {key: whole[x][key] - whole[y][key] for key in ('ground', 'air')}
        out['by_sun'] = {c: {k: whole[x]['by_sun'][c][k] - whole[y]['by_sun'][c][k] for k in ('ground', 'air')}
                         for c, _, _ in SUN_CLASSES}
        out['air_by_layer'] = list(np.subtract(whole[x]['air_by_layer'], whole[y]['air_by_layer']))
        return out
    log('the reduction against the fine grid')
    check = reduction_check(fine, ch, edges, p_edges)
    # Column masses per square metre of ground and the low Sun's slant paths.
    mass = col['layer_mass_kg_m2']
    water = mass * water_fraction(col['layer_x_h2o'], air_spec.molar_mass)
    mid = 0.5 * (edges[1:] + edges[:-1])
    dz = np.diff(edges)
    elev = np.radians([0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 60.0, 90.0])
    mu = np.sin(elev)
    half = np.r_[0.0, ps.half_levels(plasim_sigma())]
    heights = np.interp(-np.log(np.clip(half * col['surface_pressure_pa'], 1e-3, None)), -np.log(p_edges), edges - R)
    solved_e = float(e[pick].sum())
    term = grid >= np.radians(TERM_FROM_DEG)
    term_g = (tables['sphere_shallow'][0] - tables['pseudo_spherical'][0])[term] / covered
    term_q = (tables['sphere_shallow'][1] - tables['pseudo_spherical'][1])[:, term] / covered
    lunar = magnification(edges, water / dz)
    log('the outgoing infrared through flat layers and spherical shells')
    lwo = longwave_optics()
    lw = outgoing_infrared(lwo['nu'], lwo['absorption'].astype(float), edges, col['layer_t_k'], float(col['t_k'][0]))
    lw_cone = lw['layers']['sphere'] - lw['layers']['sphere_shallow']
    lw_curved = lw['layers']['sphere_shallow'] - lw['layers']['flat']
    weights = onto_gcm_layers(np.diag(mass), p_edges)                 # each shell's mass in each GCM layer
    layer_t = (weights @ col['layer_t_k']) / weights.sum(1)
    out = dict(
        schema=SCHEMA, world=WORLD,
        evidence=(f'A runnable calculation of clear, cloudless, aerosol-free air over ground of albedo {albedo}, without '
                  'refraction and without empirical check. The design column\'s line-by-line optics from 360 nm to 5 um '
                  f'(the sky model\'s {FINE_STEP} cm-1 visible grid, the radiative-convective model\'s near infrared) '
                  f'are grouped into {e.size} channels that keep the absorbed energy within {TOLERANCE:.0%} on slant '
                  f'columns (reduction_check). The {pick.size} channels below {SOLVED_NM:.0f} nm with whole-column '
                  f'absorption below {OPAQUE:g} take the spherical sky model\'s multiple scattering at its {SKY_QUALITY} '
                  'quality; the others take the direct beam through the same spherical shells with the diffuse and '
                  'reflected light of a pseudo-spherical two-stream column, so their twilight holds no scattered '
                  'light. The flat and pseudo-spherical columns are the radiative-convective model\'s delta two-stream '
                  'on the same 100 layers; PlaSim\'s is a copy of its clear-sky arithmetic on the solved column at the '
                  'design run\'s levels and settings.'),
        reading_rule=('Global means are W/m2 over the whole sphere from the sunlight inside 360 nm-5 um, '
                      f'{covered:.1f} W/m2 at the beam, scaling with the run\'s sunlight; PlaSim\'s shares are applied to '
                      'the same sunlight. The sphere is the reference; sphere_shallow is its light on the flat column\'s '
                      'own layers, without the shells\' widening above each square metre of ground. '
                      'sphere_minus_pseudo_spherical is what the curved geometry adds to columns of the same air and '
                      'optics, by class of solar elevation at the ground point and by GCM layer (top first), the sum of '
                      'cone_w_m2 (the widening) and curved_light_w_m2 (the twilight, the side light and the curved '
                      'scattered light); pseudo_spherical_minus_plasim is how PlaSim\'s calibrated arithmetic divides the '
                      'same sunlight differently; sphere_minus_plasim is the whole clear-sky difference. by_elevation '
                      'gives the absorption per square metre of ground per unit of sunlight at the beam (the whole '
                      'covered spectrum), from which the sphere integral is 2 * integral of x d(sin a).'),
        sunlight=dict(covered_w_m2=covered, visible_w_m2=float(e[bands['visible']].sum()),
                      near_infrared_w_m2=float(e[bands['near_infrared']].sum()), band_nm=[360.0, 5000.0],
                      solved_w_m2=solved_e, plasim_share_below_750_nm=PLASIM_SOLAR1),
        channels=dict(count=int(e.size), solved=int(pick.size), tolerance=TOLERANCE, opaque_column_depth=OPAQUE,
                      solved_below_nm=SOLVED_NM, sky_quality=SKY_QUALITY, batches=len(paths),
                      orders_converged=all(m['order_converged'] for m in meta),
                      solve_seconds=float(sum(m['seconds'] for m in meta))),
        reduction_check=check,
        column=dict(
            mass_kg_m2=dict(sphere=float(mass @ (mid / R) ** 2), flat=float(mass.sum()),
                            gcm=float(column['ps'] / column['gravity'])),
            mass_mean_height_km=float(mass @ (mid - R) / mass.sum() / 1e3),
            water_kg_m2=float(water.sum()), surface_pressure_pa=float(col['surface_pressure_pa']),
            gcm_surface_pressure_pa=float(column['ps']), surface_albedo=albedo,
            slant=dict(elevation_deg=list(np.degrees(elev)),
                       air=list(air_mass(edges, mass / dz, elev)), water=list(air_mass(edges, water / dz, elev)),
                       plasim=list(35.0 / np.sqrt(1.0 + 1224.0 * mu * mu)), flat=list(1.0 / np.maximum(mu, 1e-300)))),
        gcm=dict(run=DESIGN_RUN, sigma=list(plasim_sigma()), half_level_height_km=list(heights / 1e3),
                 rayleigh_scale=er.MODEL['rayleigh_scale'], sw_water_absorption=er.MODEL['sw_water_absorption'],
                 rcoeff=column['rcoeff'], gravity_m_s2=column['gravity']),
        global_mean_w_m2=dict(by_band=by_band, whole_spectrum=whole),
        sphere_minus_pseudo_spherical_w_m2=difference('sphere', 'pseudo_spherical'),
        cone_w_m2=difference('sphere', 'sphere_shallow'),
        curved_light_w_m2=difference('sphere_shallow', 'pseudo_spherical'),
        pseudo_spherical_minus_plasim_w_m2=difference('pseudo_spherical', 'plasim'),
        sphere_minus_plasim_w_m2=difference('sphere', 'plasim'),
        gcm_term=dict(
            reading_rule=('What the curved light adds to a flat column of the same air (curved_light_w_m2): the '
                          'sphere\'s light on the flat column\'s own layers, less the pseudo-spherical column\'s, per '
                          'square metre of ground per unit of sunlight at the beam, by the solar elevation at the ground '
                          'point, to the ground and to each GCM layer (top first). Nothing outside the elevations given. '
                          'Added to the GCM\'s own shortwave it keeps PlaSim\'s calibrated column and supplies the '
                          'twilight, the low Sun\'s side light and the curved scattered light. The cone (cone_w_m2), the '
                          'more sunlight the widening shells above each square metre take, is left out: those shells '
                          'emit their infrared from the same wider area, which a flat model\'s longwave leaves out too.'),
            elevation_deg=list(np.degrees(grid[term])), ground=term_g, layers=term_q,
            global_mean_w_m2=dict(ground=covered / 4.0 * float(sphere_share(term_g, grid[term])),
                                  air=covered / 4.0 * float(sphere_share(term_q.sum(0), grid[term])))),
        cm1_term=dict(
            reading_rule=('gcm_term on the 100 shells of the design column instead of the GCM\'s layers, for a model '
                          'on heights: what the curved light adds to each shell of a flat column (bottom-up, between the '
                          'heights given) per square metre of ground per unit of sunlight at the beam, by solar '
                          'elevation; the ground as in gcm_term.'),
            height_edges_km=list((edges - R) / 1e3), elevation_deg=list(np.degrees(grid[term])), ground=term_g,
            shells=shell_delta[:, term] / covered),
        cm1_term_cone=dict(
            reading_rule=('The cone the twilight terms leave out, by shell (bottom-up) and solar elevation, as in '
                          'cm1_term: the more sunlight the widening shells above each square metre of ground take.'),
            shells=shell_cone[:, term] / covered),
        gcm_term_cone=dict(
            reading_rule='The cone the GCM term leaves out, by GCM layer (top first) and solar elevation, as in gcm_term.',
            layers=(tables['sphere'][1] - tables['sphere_shallow'][1])[:, term] / covered),
        longwave=dict(
            reading_rule=('The infrared leaving the top of the design column (clear sky, its temperatures), W/m2 per '
                          'square metre of ground: through flat layers, and through the same layers as spherical '
                          'shells, which take in more sunlight (the cone and the curved light) and lose more infrared '
                          'by the same geometry: the shells\' widening (cone) and the escape below the horizontal past '
                          'the Moon\'s limb (curved). By GCM layer (top first) and by shell (bottom-up) each with its '
                          'design temperature, so a model may scale them with its own; the ground\'s share apart.'),
            olr_w_m2=dict(flat=lw['flat'], sphere=lw['sphere'], sphere_shallow=lw['sphere_shallow']),
            cone_w_m2=float(lw_cone.sum()), curved_w_m2=float(lw['sphere_shallow'] - lw['flat']),
            ground_curved_w_m2=lw['ground']['sphere'] - lw['ground']['flat'], ground_temperature_k=float(col['t_k'][0]),
            gcm_layers=dict(temperature_k=layer_t, cone=onto_gcm_layers(lw_cone[:, None], p_edges)[:, 0],
                            curved=onto_gcm_layers(lw_curved[:, None], p_edges)[:, 0]),
            shells=dict(temperature_k=col['layer_t_k'], cone=lw_cone, curved=lw_curved),
            geometry_net_w_m2=dict(
                cone=(float(sphere_share((tables['sphere'][1] - tables['sphere_shallow'][1]).sum(0), grid)) / 4.0
                      - float(lw_cone.sum())),
                curved=(float(sphere_share(term_g, grid[term]) + sphere_share(term_q.sum(0), grid[term])) * covered / 4.0
                        - float(lw['sphere_shallow'] - lw['flat'])))),
        magnification=lunar,
        by_elevation=dict(elevation_deg=list(np.degrees(grid)),
                          **{k: dict(ground=g / covered, layers=a / covered) for k, (g, a) in tables.items()},
                          plasim_reflected=p_reflected),
        inputs={**sources, str(STANDARD.relative_to(ROOT)): digest(STANDARD), **{str(p): digest(p) for p in paths},
                **{str(Path(m.__file__).resolve().relative_to(ROOT)): digest(m.__file__)
                   for m in (solved_optics, solved_transport, surface, sw, optics)}},
        producer=dict(domain='climate', files={f'gcm/{Path(m.__file__).name}': digest(m.__file__)
                                               for m in (sys.modules[__name__], ps, er)}))
    return rounded(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('task', choices=('optics', 'longwave', 'solve', 'ledger'))
    args = ap.parse_args(argv)
    if args.task == 'optics':
        near_infrared_optics()
        return 0
    if args.task == 'longwave':
        longwave_optics()
        return 0
    if args.task == 'solve':
        fine, _ = fine_inputs()
        ch = channels(fine)
        del fine
        log(f'{ch["energy"].size} channels, {int(ch["solved"].sum())} to solve')
        sky_batches(ch, solve=True)
        return 0
    out = ledger()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix('.tmp')
    tmp.write_text(json.dumps(out, indent=1) + '\n')
    tmp.replace(OUT)
    log(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
