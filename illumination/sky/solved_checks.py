"""Reproduce numerical checks and complete residual scattering orders.

python -m illumination.sky.solved_checks --grid --line-grid --monte-carlo --finish
After exporting atlases: python -m illumination.sky.solved_checks --angular --summarize
Each flag can be run separately against the ignored spherical run directory.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from illumination.sky.solved_optics import fine_optics,reduce_spectrum,digest
from illumination.sky.solved_reference import check,spectral_flux_check
from illumination.sky.solved_sky import DEFAULT_DIRECTORY,WORLDS,SUNS,dump,load,photometry,horizontal_integral
from illumination.sky.solved_transport import build_paths,energy_audit,evaluate,quadrature,solve
from illumination.sky.solver import transport
from illumination.cloud_light.model import MolecularColumn


def line_grid_checks(directory):
    result = {}
    for world in WORLDS:
        fine,meta = fine_optics(world,directory)
        finer,_ = fine_optics(world,directory,step_cm1=.025)
        geometry = MolecularColumn(meta['radius_m'],fine['height'],
                                   np.zeros((len(fine['height'])-1,1)),np.array([550.]))
        rows = []
        for h in (0.,10000.,30000.,80000.):
            for sun in (-10.,-5.,-2.,0.,1.,3.,10.,30.,90.):
                lengths,blocked = geometry.path(h,np.sin(np.radians(sun)))
                if blocked:
                    continue
                factors = lengths/np.diff(fine['height'])
                a = np.exp(-factors@(fine['scattering']+fine['absorption']))@fine['xyz'][:,1]
                b = np.exp(-factors@(finer['scattering']+finer['absorption']))@finer['xyz'][:,1]
                if b>1e-6:
                    rows.append(dict(height_m=h,sun_deg=sun,coarse_lux=float(a),fine_lux=float(b),relative_difference=float(a/b-1)))
        result[world] = dict(coarse_spacing_cm1=.05,fine_spacing_cm1=.025,samples=rows,
                             maximum_absolute_relative_difference=max(abs(r['relative_difference']) for r in rows))
    dump(directory/'line_grid_checks.json',result)


def angular_checks(directory):
    result = {}
    suns = np.array([90.,10.,1.,0.,-2.,-6.,-12.])
    elevations = 90*np.linspace(0,1,81)**2
    azimuths = np.linspace(0,180,129)
    for world in WORLDS:
        solution,settings = load(directory/f'{world}_standard.npz')
        radiance = evaluate(solution,settings,suns,elevations,azimuths)
        fine = horizontal_integral(radiance@solution['xyz'][:,1],elevations,azimuths)
        with np.load(directory/f'{world}_atlas.npz',allow_pickle=False) as atlas:
            indices = [np.flatnonzero(atlas['suns']==a)[0] for a in suns]
            coarse = horizontal_integral(atlas['xyz'][indices,:,:,1],atlas['elevations'],atlas['azimuths'])
        result[world] = dict(coarse_shape=[41,65],fine_shape=[81,129],samples=[dict(
            sun_deg=float(a),coarse_diffuse_lux=float(c),fine_diffuse_lux=float(f),relative_difference=float(c/f-1))
            for a,c,f in zip(suns,coarse,fine,strict=True)],
            maximum_absolute_relative_difference=float(abs(coarse/fine-1).max()))
        print(world,'angular integral refinement',result[world]['maximum_absolute_relative_difference'],flush=True)
        dump(directory/'angular_readout_checks.json',result)


def grid_checks(directory):
    fine,meta = fine_optics('moon_1.2atm',directory)
    optical = reduce_spectrum(fine,10.,2)
    indices = [4,23,43]
    for k in ('scattering','absorption'):
        optical[k] = optical[k][:,indices]
    for k in ('xyz','energy','wavelength','band'):
        optical[k] = optical[k][indices]
    for quality in ('draft','standard','refined'):
        solution,settings = solve(optical,meta['radius_m'],quality,path=directory/f'grid_{quality}.npz')
        audit = energy_audit(solution,settings)
        audit.update(spot_suns=[90,30,10,0,-3,-6,-12,-18],spot_elevations=[0,10,45,90],spot_azimuths=[0,90,180],
                     wavelength_nm=solution['wavelength'].tolist())
        audit['radiance'] = evaluate(solution,settings,audit['spot_suns'],audit['spot_elevations'],audit['spot_azimuths']).tolist()
        dump(directory/f'grid_{quality}_audit.json',audit)


def monte_carlo_checks(directory):
    solution,settings = load(directory/'grid_standard.npz')
    rows = []
    for j,sun,elev,az,flux in [(0,90,90,90,False),(1,90,90,90,False),(1,10,10,0,False),
                              (1,-6,10,0,False),(2,-6,10,0,False),(1,90,90,90,True),
                              (1,10,90,90,True),(1,-6,90,90,True)]:
        row = check(solution,settings,j,sun,elev,az,photons=120000,flux=flux)
        value = (np.interp(np.radians(sun),solution['a'],solution['moments'][0,:,j,5]) if flux else
                 evaluate(solution,settings,[sun],[elev],[az])[0,0,0,j])
        row.update(deterministic=float(value),relative_difference=float(value/row['mean']-1),
                   difference_standard_errors=float((value-row['mean'])/row['standard_error']))
        rows.append(row)
    dump(directory/'monte_carlo.json',dict(schema='terluna.illumination.solved-sky-monte-carlo/1',
         source='grid_standard.npz',finite_sun_in_deterministic=True,point_sun_in_reference=True,spots=rows))
    rows = []
    for world in WORLDS:
        fine,meta = fine_optics(world,directory)
        suns = (90.,10.,-6.,-12.,0.,-2.) if world=='earth_control' else (90.,10.,-6.,-12.,1.)
        for sun in suns:
            extra = world=='earth_control' and sun in (0.,-2.)
            row = spectral_flux_check(fine,meta['radius_m'],sun,photons=480000 if extra else 240000,
                                      seed=118 if extra else 117)
            row['world'] = world
            rows.append(row)
    dump(directory/'spectral_monte_carlo.json',dict(schema='terluna.illumination.solved-sky-spectral-monte-carlo/1',
         source='Fine 1-nm/line-by-line optics; photopic wavelength sampling on every path; point Sun',spots=rows))


def finish_orders(path,maximum_order=160):
    """Continue each unconverged channel from its saved final scattering order."""
    solution,settings = load(path)
    if settings['order_converged']:
        return
    start = time.monotonic()
    r,a,sr,sa = [solution[k] for k in ('r','a','sr','sa')]
    roi = a>=np.deg2rad(-24)
    total,last = solution['moments'],solution['last_order']
    tolerance = settings['increment_tolerance']
    per_channel = np.max(last[0,roi,:,5]/np.maximum(total[0,roi,:,5],1e-12),axis=0)
    per_maximum = np.max(last[:,:,:,0],axis=(0,1))/np.maximum(np.max(total[:,:,:,0],axis=(0,1)),1e-30)
    active = np.flatnonzero((per_channel>=tolerance)|(per_maximum>=tolerance))
    q = quadrature(r,r[0],settings['angular_order_per_interval'],settings['azimuth_samples'])
    orders = np.asarray(settings['channel_orders'])
    continuation = []
    for channel in active:
        paths = build_paths(r,q[0],solution['edges'],solution['scattering'][:,channel:channel+1],
                            solution['absorption'][:,channel:channel+1],solution['shells'])
        beam = np.ascontiguousarray(solution['beam'][:,:,channel:channel+1,:])
        previous = np.ascontiguousarray(last[:,:,channel:channel+1,:])
        for order in range(orders[channel]+1,maximum_order+1):
            current = transport(previous,beam,r,a,sr,sa,*q,*paths,settings['albedo'],False)
            total[:,:,channel:channel+1,:] += current
            last[:,:,channel:channel+1,:] = current
            orders[channel] = order
            rel = float(np.max(current[0,roi,0,5]/np.maximum(total[0,roi,channel,5],1e-12)))
            largest = float(current[:,:,:,0].max()/max(total[:,:,channel,0].max(),1e-30))
            previous = current
            if rel<tolerance and largest<tolerance:
                break
        if rel>=tolerance or largest>=tolerance:
            raise ValueError(f'Scattering order cap reached for channel {channel}')
        continuation.append(dict(channel=int(channel),final_order=int(order),ground_increment_fraction=rel,
                                 maximum_increment_fraction=largest))
    settings.update(channel_orders=orders.tolist(),order_converged=True,continued_orders=continuation,
                    continuation_seconds=time.monotonic()-start,continuation_runner='illumination/sky/solved_checks.py')
    temporary = path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**solution,meta=json.dumps(settings))
    temporary.replace(path)
    print(path.name,'completed orders',continuation,flush=True)


def summarize(directory):
    spectral = json.loads((directory/'spectral_monte_carlo.json').read_text())
    rows,grids = [],{}
    for world in WORLDS:
        s,m = load(directory/f'{world}_standard.npz')
        coarse,_ = load(directory/f'{world}_draft.npz')
        _,_,flux = photometry(s,SUNS)
        _,_,coarse_flux = photometry(coarse,SUNS)
        difference = coarse_flux[:,1]/np.maximum(flux[:,1],1e-20)-1
        with np.load(directory/f'{world}_atlas.npz',allow_pickle=False) as atlas:
            if not np.array_equal(atlas['suns'],SUNS):
                raise ValueError('Atlas solar grid differs')
            angular_flux = horizontal_integral(atlas['xyz'][:,:,:,1],atlas['elevations'],atlas['azimuths'])
        grids[world] = dict(maximum_daylight_diffuse_relative_difference=float(abs(difference[SUNS>=0]).max()),
                           maximum_twilight_through_minus18_relative_difference=float(abs(difference[(SUNS<0)&(SUNS>=-18)]).max()))
        for row in spectral['spots']:
            if row['world']!=world:
                continue
            i = np.flatnonzero(SUNS==row['sun_deg'])[0]
            value = float(angular_flux[i])
            mean,error = row['mean_diffuse_lux'],row['standard_error_lux']
            relative_error = error/mean if mean>0 else None
            rows.append(dict(**row,deterministic_diffuse_lux=value,
                             moment_grid_diffuse_lux=float(flux[i,1]),
                             moment_grid_relative_difference=float(flux[i,1]/mean-1) if mean>0 else None,
                             relative_difference=value/mean-1 if mean>0 else None,
                             difference_standard_errors=(value-mean)/error if error>0 else None,
                             relative_standard_error=relative_error,
                             percent_comparison_resolved=relative_error is not None and relative_error<=.02,
                             sampling_status='sampled positive contributions' if mean>0 else
                             'Zero sampled contributions; additional importance sampling is needed to bound this faint flux'))
    standard = json.loads((directory/'grid_standard_audit.json').read_text())
    refined = json.loads((directory/'grid_refined_audit.json').read_text())
    change = np.asarray(standard['radiance'])/np.maximum(refined['radiance'],1e-30)-1
    grid = dict(wavelength_nm=refined['wavelength_nm'],sun_deg=refined['spot_suns'],
                maximum_radiance_relative_difference_by_sun=np.max(abs(change),axis=(1,2,3)).tolist(),
                standard_energy_residual=standard['residual'],refined_energy_residual=refined['residual'])
    files = [directory/f'{world}_{kind}.npz' for world in WORLDS for kind in ('draft','standard','production_optics','atlas')]
    files += [directory/name for name in ('grid_standard.npz','grid_refined.npz','spectral_monte_carlo.json','monte_carlo.json')]
    result = dict(schema='terluna.illumination.solved-sky-numerical-checks/1',
                  inputs={p.name:digest(p) for p in files},spectral_monte_carlo=rows,
                  broadband_grid_changes=grids,three_channel_refinement=grid,
                  maximum_resolved_broadband_monte_carlo_relative_difference=max(
                      abs(r['relative_difference']) for r in rows if r['percent_comparison_resolved']),
                  percent_comparison_relative_standard_error_limit=.02,
                  all_monte_carlo_paths_completed=all(r['truncated_paths']==0 for r in rows))
    dump(directory/'numerical_checks.json',result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=DEFAULT_DIRECTORY)
    for flag in ('grid','line-grid','angular','monte-carlo','finish','summarize'):
        parser.add_argument('--'+flag,action='store_true')
    args = parser.parse_args()
    if args.grid:
        grid_checks(args.directory)
    if args.line_grid:
        line_grid_checks(args.directory)
    if args.angular:
        angular_checks(args.directory)
    if args.monte_carlo:
        monte_carlo_checks(args.directory)
    if args.finish:
        for world in WORLDS:
            for quality in ('draft','standard'):
                path = args.directory/f'{world}_{quality}.npz'
                finish_orders(path)
                s,m = load(path)
                dump(args.directory/f'{world}_{quality}_energy.json',energy_audit(s,m))
    if args.summarize:
        print(json.dumps(summarize(args.directory),indent=2))


if __name__=='__main__':
    main()
