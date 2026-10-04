"""Compute and export a shielded spherical sky on the solved atmospheric columns.

OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=3 python -m illumination.sky.solved_sky --compute
python -m illumination.sky.solved_checks --finish
python -m illumination.sky.solved_sky --export

Bulk fields live under the ignored research-run directory. The compact JSON
product records those archives, their schemas and hashes, and numerical checks.
"""
from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np
import numba

from illumination.sky.solved_optics import ROOT, adaptive_spectrum, digest, fine_optics
from illumination.sky.solved_transport import energy_audit, evaluate, solar_table, solve
from shared.provenance import constants_used

SCHEMA = 'terluna.illumination.solved-spherical-sky/1'
ARCHIVE_SCHEMA = 'terluna.illumination.solved-spherical-sky-atlas/1'
DEFAULT_DIRECTORY = ROOT / 'research/runs/optical_comfort/spherical'
PRODUCT = ROOT / 'illumination/sky/results/solved_sky.json'
WORLDS = ('moon_1.2atm','earth_control')
SUNS = np.array([-30.,-24.,-20.,-18.,-16.,-14.,-12.,-10.,-8.,-6.,-5.,-4.,-3.,-2.,-1.,0.,
                 .5,1.,2.,3.,5.,7.,10.,15.,20.,30.,45.,60.,75.,90.])
MODEL_FILES = ('illumination/sky/solved_optics.py','illumination/sky/solved_transport.py',
               'illumination/sky/solved_reference.py','illumination/sky/solved_sky.py',
               'illumination/sky/solved_checks.py',
               'illumination/sky/solver.py','illumination/sky/reference_mc.py',
               'illumination/sky/colour_matching.py')


def load(path):
    with np.load(path,allow_pickle=False) as data:
        return {k:data[k].copy() for k in data.files if k!='meta'},json.loads(str(data['meta']))


def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def compute(directory,worlds=WORLDS):
    directory = Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    for world in worlds:
        fine,metadata = fine_optics(world,directory)
        reduced,audit = adaptive_spectrum(fine,metadata['radius_m'],tolerance=.01)
        np.savez_compressed(directory/f'{world}_production_optics.npz',**reduced,meta=json.dumps(metadata))
        dump(directory/f'{world}_production_beam.json',audit)
        del fine
        for quality in ('draft','standard'):
            solution,settings = solve(reduced,metadata['radius_m'],quality,path=directory/f'{world}_{quality}.npz')
            dump(directory/f'{world}_{quality}_energy.json',energy_audit(solution,settings))
            del solution


def photometry(solution,suns):
    angles = np.radians(suns)
    direct = solar_table(solution['r'][:1],angles,solution['edges'],
                         solution['scattering']+solution['absorption'])[0]
    diffuse = np.array([np.interp(angles,solution['a'],solution['moments'][0,:,l,5])
                        for l in range(len(solution['energy']))]).T
    return direct[:,:,0]@solution['xyz'],direct[:,:,1]@solution['xyz'],diffuse@solution['xyz']


def horizontal_integral(luminance,elevations,azimuths):
    """Exact horizontal integral of bilinear interpolation in sin(elevation), azimuth."""
    mu,az = np.sin(np.radians(elevations)),np.radians(azimuths)
    ring = np.sum((luminance[:,:,:-1]+luminance[:,:,1:])*np.diff(az),axis=2)
    a,b = mu[:-1],mu[1:]
    return np.sum((b-a)/6*((2*a+b)*ring[:,:-1]+(a+2*b)*ring[:,1:]),axis=1)


def export(directory,worlds=WORLDS,reuse_atlases=False):
    directory = Path(directory)
    previous = json.loads(PRODUCT.read_text()) if reuse_atlases else None
    if previous is not None and previous['schema']!=SCHEMA:
        raise ValueError('Unsupported cached sky product schema')
    results = {}
    for world in worlds:
        solution,settings = load(directory/f'{world}_standard.npz')
        if not settings['order_converged']:
            raise ValueError('Complete the residual scattering orders with solved_checks --finish')
        _,optical = load(directory/f'{world}_production_optics.npz')
        coarse,coarse_settings = load(directory/f'{world}_draft.npz')
        dni,direct,moment_diffuse = photometry(solution,SUNS)
        _,_,coarse_diffuse = photometry(coarse,SUNS)
        del coarse
        elevations = 90*np.linspace(0,1,41)**2
        azimuths = np.linspace(0,180,65)
        archive = directory/f'{world}_atlas.npz'
        if reuse_atlases:
            if digest(archive)!=previous['worlds'][world]['atlas']['sha256']:
                raise ValueError('Cached atlas hash differs from its admitted product')
            for file in ('illumination/sky/solved_transport.py','illumination/sky/solver.py','illumination/sky/colour_matching.py'):
                if previous['producer']['files'][file]!=digest(ROOT/file):
                    raise ValueError('Ray evaluator changed; export fresh atlases')
            with np.load(archive,allow_pickle=False) as saved:
                head = json.loads(str(saved['meta']))
                if head['source_solution_sha256']!=digest(directory/f'{world}_standard.npz'):
                    raise ValueError('Source solution differs from cached atlas')
                if any(not np.array_equal(saved[key],value) for key,value in
                       (('suns',SUNS),('elevations',elevations),('azimuths',azimuths))):
                    raise ValueError('Cached atlas angular grid differs')
                xyz = saved['xyz'].copy()
        else:
            spectrum = evaluate(solution,settings,SUNS,elevations,azimuths)
            xyz = spectrum@solution['xyz']
            del spectrum
        diffuse = np.array([horizontal_integral(xyz[:,:,:,c],elevations,azimuths) for c in range(3)]).T
        header = dict(schema=ARCHIVE_SCHEMA,world=world,albedo=settings['albedo'],
                      angular_interpolation='bilinear in sin(elevation) and symmetric solar-relative azimuth',
                      radiance_units='CIE XYZ photopic scale; Y in cd m^-2',
                      flux_units='CIE XYZ photopic scale; Y in lux',
                      diffuse_flux_rule='Exact horizontal integral of the stored angular interpolant',
                      source_solution_sha256=digest(directory/f'{world}_standard.npz'))
        np.savez_compressed(archive,suns=SUNS,elevations=elevations,azimuths=azimuths,xyz=xyz,
                            direct_normal_xyz=dni,direct_horizontal_xyz=direct,diffuse_horizontal_xyz=diffuse,
                            moment_diffuse_horizontal_xyz=moment_diffuse,
                            meta=json.dumps(header))
        rows = []
        for i,sun in enumerate(SUNS):
            total = direct[i]+diffuse[i]
            rows.append(dict(sun_deg=float(sun),direct_normal_lux=float(dni[i,1]),
                             direct_horizontal_lux=float(direct[i,1]),diffuse_horizontal_lux=float(diffuse[i,1]),
                             total_horizontal_lux=float(total[1]),
                             diffuse_xy=(diffuse[i,:2]/diffuse[i].sum()).tolist() if diffuse[i].sum()>0 else None,
                             zenith_luminance_cd_m2=float(xyz[i,-1,0,1]),
                             sunsetward_horizon_luminance_cd_m2=float(xyz[i,0,0,1]),
                             moment_diffuse_horizontal_lux=float(moment_diffuse[i,1]),
                             coarse_grid_diffuse_relative_difference=float(coarse_diffuse[i,1]/max(moment_diffuse[i,1],1e-30)-1),
                             atlas_over_moment_flux_relative_difference=float(diffuse[i,1]/max(moment_diffuse[i,1],1e-30)-1)))
        results[world] = dict(optical_column=optical,settings=settings,coarse_grid_settings=coarse_settings,
                              spectral_check=json.loads((directory/f'{world}_production_beam.json').read_text()),
                              energy=json.loads((directory/f'{world}_standard_energy.json').read_text()),
                              atlas=dict(path=str(archive.relative_to(ROOT)) if archive.is_relative_to(ROOT) else str(archive),
                                         sha256=digest(archive),schema=ARCHIVE_SCHEMA),samples=rows)
        print(world,'noon',rows[-1],'horizon',rows[15],flush=True)
    checks = {}
    for name in ('grid_draft_audit.json','grid_standard_audit.json','grid_refined_audit.json','monte_carlo.json',
                 'spectral_monte_carlo.json','line_grid_checks.json','angular_readout_checks.json','numerical_checks.json'):
        if (directory/name).is_file():
            checks[name] = json.loads((directory/name).read_text())
    product = dict(schema=SCHEMA,producer=dict(lane='illumination',runner='illumination/sky/solved_sky.py',
                                              files={f:digest(ROOT/f) for f in MODEL_FILES},
                                              constants=constants_used(MODEL_FILES)),
                   software=dict(python=platform.python_version(),numpy=np.__version__,numba=numba.__version__),
                   evidence='Scalar clear-sky multiple scattering on the current solved molecular columns. '
                            'The Moon uses titania-stack transmission times 0.95; Earth uses unfiltered sunlight. '
                            'Ground albedo is 0.1. Numerical checks describe this conditional atmosphere.',
                   reading_rule='Use each world’s own absolute sky radiance and flux. Sun angles are geometric. '
                                'Atlas Y is luminance; preserve signed XYZ-to-RGB components if converting colour. '
                                'Surface diffuse flux is integrated from the angular atlas; the internal moment-grid '
                                'flux is retained separately. The global energy ledger uses the moment-grid quadrature. '
                                'Interpolation error grows in twilight: compare stored solar samples first.',
                   assumptions=['Horizontally uniform solved one-dimensional atmosphere',
                                'Piecewise homogeneous optical layers, straight rays, uniform finite solar disk',
                                'Scalar Rayleigh phase, Lambertian ground, 360–830 nm photopic calculation',
                                'Clouds, aerosols, polarization, refraction and terrain require additional models',
                                'Solar-weighted spectral reduction tested separately from spatial convergence'],
                   worlds=results,additional_checks=checks)
    PRODUCT.parent.mkdir(parents=True,exist_ok=True)
    dump(PRODUCT,product)
    return product


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=DEFAULT_DIRECTORY)
    parser.add_argument('--compute',action='store_true')
    parser.add_argument('--export',action='store_true')
    parser.add_argument('--reuse-atlases',action='store_true',help='reuse existing radiance only after hash and source checks')
    parser.add_argument('--world',choices=WORLDS)
    args = parser.parse_args()
    worlds = (args.world,) if args.world else WORLDS
    if args.compute:
        compute(args.directory,worlds)
    if args.export:
        export(args.directory,worlds,args.reuse_atlases)
    if not (args.compute or args.export):
        parser.error('Choose --compute, --export, or both')


if __name__ == '__main__':
    main()
