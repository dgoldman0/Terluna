"""Visible molecular optics from the solved columns, with checked spectral reduction.

The fine grid retains each layer's H2O/O2 line absorption above 500 nm. Within
each wavelength band we sort whole-column absorption and average layer optical
depths in solar-energy quantiles. Every reduced channel carries its own energy
and CIE XYZ weights. Direct slant beams provide a fine-grid comparison before
these channels are used for multiple scattering.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from atmosphere.middle_atmosphere import photolysis as ph
from atmosphere.radiative_convective import optics, shortwave as sw
from illumination.cloud_light.model import MolecularColumn
from illumination.sky.colour_matching import cmf
from illumination.surface_light import model as surface

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = 'terluna.illumination.solved-sky-optics/1'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs(column):
    """Hash the producer, atmospheric products and restored spectroscopy bytes."""
    files = [
        'illumination/sky/solved_optics.py', 'illumination/sky/colour_matching.py',
        'illumination/sky/atmospheres.py', 'illumination/cloud_light/model.py',
        'illumination/surface_light/model.py', 'shared/constants.json', 'shared/constants.py',
        'shared/external_inputs.py',
        'atmosphere/radiative_convective/optics.py', 'atmosphere/radiative_convective/spectroscopy.py',
        'atmosphere/radiative_convective/shortwave.py', 'atmosphere/radiative_convective/thermodynamics.py',
        'atmosphere/radiative_convective/climate.py', 'atmosphere/radiative_convective/longwave.py',
        'atmosphere/radiative_convective/ck.py', 'atmosphere/radiative_convective/fetch_inputs.py',
        'atmosphere/middle_atmosphere/equilibrium.py', 'atmosphere/middle_atmosphere/photolysis.py',
        'atmosphere/middle_atmosphere/chemistry.py', 'atmosphere/middle_atmosphere/fetch_inputs.py',
        'atmosphere/middle_atmosphere/results/middle_atmosphere.json',
        'protection/transmission/results/uv_transmission.json',
    ]
    for name, _ in column.stored:
        files += [f'atmosphere/middle_atmosphere/results/{folder}/{name}.{ext}'
                  for folder, ext in (('profiles', 'csv'), ('cases', 'json'))]
    external = {}
    for folder in ('atmosphere/radiative_convective', 'atmosphere/middle_atmosphere'):
        manifest = ROOT / folder / 'inputs.json'
        files.append(str(manifest.relative_to(ROOT)))
        data = json.loads(manifest.read_text())
        for record in data['files']:
            p = manifest.parent / data['input_directory'] / record['name']
            if p.is_file():
                sha = digest(p)
                if sha != record['sha256'] or p.stat().st_size != record['bytes']:
                    raise ValueError(f'Spectral input checksum differs: {p}')
                external[str(p.relative_to(ROOT))] = sha
    return dict(files={f: digest(ROOT / f) for f in files}, external_inputs=external,
                software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


def fine_optics(world, directory, step_cm1=0.05, processes=2):
    """Layer optical depths and incident energy on the 1-nm / line-by-line grids.

    All 100 solved-column layers are retained. Above 500 nm the line model
    supplies H2O, O2, water continuum and collision-induced absorption; the
    photolysis tables add each solved trace gas once.
    """
    column = surface.COLUMNS[world]
    provenance = inputs(column)
    key = hashlib.sha256(json.dumps([world, step_cm1, provenance], sort_keys=True).encode()).hexdigest()
    path = Path(directory) / f'{world}_fine_{key[:16]}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            if json.loads(str(data['meta']))['key'] != key:
                raise ValueError('Fine-optics cache fingerprint differs')
            return {k: data[k].copy() for k in data.files if k != 'meta'}, json.loads(str(data['meta']))
    col, air, densities, profiles = surface.build(column)
    light_name = 'design' if world.startswith('moon') else 'unfiltered'
    shield, shield_hash = surface.sunlight(light_name)
    sunlight = ph.Sunlight(shield)
    abs_uv, sca_uv, _ = sunlight.optics(col, densities)
    uv = (ph.CENTRES_NM >= 360) & (ph.CENTRES_NM < 500)
    cfg = replace(surface.LBL, nu_min=1e7 / 830.0, step=step_cm1, processes=processes)
    print(f'{world}: computing line absorption on {len(col["z_m"])-1} layers', flush=True)
    nu, abs_lbl = optics.optical_depth(col, cfg)
    wave = 1e7 / nu
    for species, density in densities.items():
        sigma = sunlight.xs.sigma(species, col['layer_t_k'])
        for layer in range(len(abs_lbl)):
            abs_lbl[layer] += density[layer] * col['layer_dz_cm'][layer] * np.interp(
                wave, ph.CENTRES_NM, sigma[layer], left=0, right=0)
    sca_lbl = optics.rayleigh_optical_depth(col, nu)
    energy, _ = sw.solar_spectrum(nu, shield=shield)
    dw = np.gradient(nu)
    dw[[0, -1]] *= 0.5
    energy *= dw
    arrays = dict(
        height=col['z_m'], wavelength=np.r_[ph.CENTRES_NM[uv], wave[::-1]],
        absorption=np.c_[abs_uv[:, uv], abs_lbl[:, ::-1]],
        scattering=np.c_[sca_uv[:, uv], sca_lbl[:, ::-1]],
        energy=np.r_[sunlight.energy[uv], energy[::-1]],
    )
    arrays['xyz'] = 683.0 * arrays['energy'][:, None] * cmf(arrays['wavelength'])
    meta = dict(schema=SCHEMA, key=key, world=world, sunlight=light_name,
                radius_m=column.planet.radius_m, profiles=profiles, shield_sha256=shield_hash,
                column=surface.column_facts(col, air, densities), producer=provenance,
                line_spacing_cm1=step_cm1, layers=len(col['z_m'])-1,
                solar_source='TSIS-1 HSRS, admitted stride-100 table; normalized by the atmosphere model',
                evidence='Solved one-dimensional molecular column and prescribed shield spectrum; '
                         'piecewise homogeneous spherical layers; scalar elastic Rayleigh scattering.')
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, **arrays, meta=json.dumps(meta))
    tmp.replace(path)
    return arrays, meta


def reduce_spectrum(fine, band_nm=10, quantiles=4):
    """Solar-weighted column-rank bins with a common ordering for every layer.

    Smooth bands below 500 nm use one channel each. Quantile membership is
    split at exact energy boundaries, preserving incident energy and XYZ sums.
    A single ordering across layers is an approximation checked on slant rays.
    """
    if band_nm <= 0 or quantiles < 1:
        raise ValueError('Positive band width and quantile count required')
    wave, energy = fine['wavelength'], fine['energy']
    grouped = {k: [] for k in ('wavelength', 'absorption', 'scattering', 'energy', 'xyz', 'band')}
    for lo in np.arange(360., 830., band_nm):
        hi = min(lo + band_nm, 830.)
        ix = np.flatnonzero((wave >= lo - 1e-9) & ((wave < hi) if hi < 830 else (wave <= hi + 1e-8)))
        if len(ix) == 0:
            continue
        ix = ix[np.argsort(fine['absorption'][:, ix].sum(axis=0), kind='stable')]
        weights = energy[ix]
        cum = np.r_[0., np.cumsum(weights)]
        count = 1 if hi <= 500 else quantiles
        for left, right in zip(np.linspace(0, cum[-1], count + 1)[:-1],
                               np.linspace(0, cum[-1], count + 1)[1:]):
            w = np.maximum(0, np.minimum(cum[1:], right) - np.maximum(cum[:-1], left))
            good = w > 0
            j, w = ix[good], w[good]
            total = w.sum()
            if total <= 0:
                continue
            grouped['wavelength'].append(float(w @ wave[j] / total))
            grouped['absorption'].append(fine['absorption'][:, j] @ w / total)
            grouped['scattering'].append(fine['scattering'][:, j] @ w / total)
            grouped['energy'].append(total)
            grouped['xyz'].append((fine['xyz'][j] * (w / energy[j])[:, None]).sum(axis=0))
            grouped['band'].append([lo, hi])
    result = {k: np.asarray(v) for k, v in grouped.items()}
    result['absorption'] = result['absorption'].T
    result['scattering'] = result['scattering'].T
    result['height'] = fine['height'].copy()
    if not np.allclose(result['xyz'].sum(axis=0), fine['xyz'].sum(axis=0), rtol=1e-10):
        raise ValueError('Spectral reduction lost source colour weights')
    return result


def beam_audit(fine, reduced, radius):
    """Compare exact layered slant beams before and after spectral reduction."""
    shape = MolecularColumn(radius, fine['height'], np.zeros((len(fine['height'])-1, 1)), np.array([550.]))
    rows = []
    for height in (0., 10000., 30000., 80000.):
        if height > fine['height'][-1]:
            continue
        for sun in (-10., -5., -2., 0., 1., 3., 10., 30., 90.):
            lengths, blocked = shape.path(height, np.sin(np.radians(sun)))
            if blocked:
                continue
            factors = lengths / np.diff(fine['height'])
            exact = np.exp(-factors @ (fine['absorption'] + fine['scattering'])) @ fine['xyz']
            approx = np.exp(-factors @ (reduced['absorption'] + reduced['scattering'])) @ reduced['xyz']
            if exact[1] > 1e-6:
                rows.append(dict(height_m=height, sun_deg=sun, reference_lux=float(exact[1]),
                                 reduced_lux=float(approx[1]), relative_lux_error=float(approx[1]/exact[1]-1),
                                 xy_error=float(np.max(abs(approx[:2]/approx.sum()-exact[:2]/exact.sum())))))
    return dict(samples=rows, maximum_absolute_relative_lux_error=max(abs(r['relative_lux_error']) for r in rows),
                maximum_absolute_xy_error=max(r['xy_error'] for r in rows))


def adaptive_spectrum(fine, radius, band_nm=10., tolerance=.002, max_channels=500):
    """Bisect energy quantiles until sampled slant-beam XYZ errors meet a budget.

    Training rays span several heights and slant columns. Summed absolute bin
    errors provide a conservative sampled budget. beam_audit uses different
    heights and angles as a separate check of the reduction.
    """
    shape = MolecularColumn(radius,fine['height'],np.zeros((len(fine['height'])-1,1)),np.array([550.]))
    factors = []
    for h in (0.,5000.,20000.,50000.,min(100000.,fine['height'][-1])):
        if h>fine['height'][-1]:
            continue
        for angle in (-12.,-7.,-3.,.3,2.,7.,20.,60.):
            lengths,blocked = shape.path(h,np.sin(np.radians(angle)))
            if not blocked:
                factors.append(lengths/np.diff(fine['height']))
    factors = np.array(factors)
    optical = fine['absorption']+fine['scattering']
    transmission = np.exp(-factors@optical)
    reference = transmission@fine['xyz']
    # Photopic relative error controls brightness. X and Z use total tristimulus
    # magnitude, keeping the colour budget useful for strongly reddened beams.
    normalization = np.repeat(np.maximum(reference.sum(axis=1),1e-4)[:,None],3,axis=1)
    normalization[:,1] = np.maximum(reference[:,1],1e-4)
    def group(ix,band):
        weights = fine['energy'][ix]
        total = weights.sum()
        ta = fine['absorption'][:,ix]@weights/total
        ts = fine['scattering'][:,ix]@weights/total
        xyz = fine['xyz'][ix].sum(axis=0)
        exact = transmission[:,ix]@fine['xyz'][ix]
        estimate = np.exp(-factors@(ta+ts))[:,None]*xyz
        error = abs(exact-estimate)/normalization
        return dict(indices=ix,band=band,energy=total,absorption=ta,scattering=ts,xyz=xyz,
                    wavelength=float(weights@fine['wavelength'][ix]/total),error=error)
    groups = []
    for lo in np.arange(360.,830.,band_nm):
        hi = min(lo+band_nm,830.)
        wave = fine['wavelength']
        ix = np.flatnonzero((wave>=lo-1e-9)&((wave<hi) if hi<830 else (wave<=hi+1e-8)))
        ix = ix[np.argsort(optical[:,ix].sum(axis=0),kind='stable')]
        if len(ix):
            groups.append(group(ix,[lo,hi]))
    budget = sum(g['error'] for g in groups)
    while budget.max()>tolerance and len(groups)<max_channels:
        worst = np.unravel_index(budget.argmax(),budget.shape)
        scores = [g['error'][worst] if len(g['indices'])>1 else -1 for g in groups]
        index = int(np.argmax(scores))
        old = groups.pop(index)
        ix = old['indices']
        cut = np.searchsorted(np.cumsum(fine['energy'][ix]),old['energy']/2)+1
        cut = min(len(ix)-1,max(1,cut))
        new = [group(ix[:cut],old['band']),group(ix[cut:],old['band'])]
        budget += new[0]['error']+new[1]['error']-old['error']
        groups.extend(new)
    groups.sort(key=lambda g:(g['band'][0],g['absorption'].sum()))
    result = {k:np.asarray([g[k] for g in groups]) for k in
              ('band','energy','xyz','wavelength','scattering','absorption')}
    result['absorption'] = result['absorption'].T
    result['scattering'] = result['scattering'].T
    result['height'] = fine['height'].copy()
    if budget.max()>tolerance:
        raise ValueError(f'Spectral channel cap reached before tolerance: {budget.max()}')
    audit = beam_audit(fine,result,radius)
    audit.update(training_rays=len(factors),maximum_training_xyz_error_bound=float(budget.max()),
                 target=tolerance,channels=len(groups),maximum_band_width_nm=band_nm)
    return result,audit
