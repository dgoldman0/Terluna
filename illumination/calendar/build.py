"""Reproduce the clear-sky calendar's compact spectral transfer product.

OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=2 python -m illumination.calendar.build --compute --export

The heavy point-source solution is ignored research data. Only its surface
readout, source weights and verification records are published. Existing solved
products remain immutable. Missing external inputs fail explicitly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from illumination.earthlight import model as earthlight
from illumination.sky.solved_optics import fine_optics, adaptive_spectrum, digest
from illumination.calendar.transport import solve, point_table
from illumination.surface_light.model import sunlight
from shared.constants import AU, EARTH_MOON_DISTANCE, EARTH_RADIUS, MOON_RADIUS, SUN_RADIUS
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/'research/runs/light_calendar'
PRODUCT = Path(__file__).parent/'results/transfer.json'
SCHEMA = 'terluna.illumination.calendar-transfer/1'


def numbers(values):
    """Ten significant digits retain faint flux without a fixed decimal floor."""
    a=np.asarray(values)
    if a.ndim==0:return float(f'{float(a):.10g}')
    return [numbers(v) for v in a]


def check_optics_provenance(meta):
    for group in ('files','external_inputs'):
        for name,expected in meta['producer'][group].items():
            path=ROOT/name
            if (path.exists() or group=='files') and digest(path)!=expected:
                raise ValueError(f'Cached optics producer/input differs: {name}')


def transport_identity(optical_file):
    return dict(optics_sha256=digest(optical_file),
                producers={p:digest(ROOT/p) for p in ['illumination/calendar/transport.py','illumination/sky/solved_transport.py',
                                                    'illumination/sky/solver.py']},
                settings=dict(quality='standard',albedo=.1,max_orders=240,tolerance=1e-6,
                              source_radius_rad=0.,convergence_min_elevation_deg=-90.))


def compute(directory=DIRECTORY):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    optical_file=directory/'optics.npz'
    if optical_file.exists():
        with np.load(optical_file,allow_pickle=False) as z:
            reduced={k:z[k].copy() for k in z.files if k!='meta'}
            meta=json.loads(str(z['meta']))
            check_optics_provenance(meta)
    else:
        fine,meta=fine_optics('moon_1.2atm',directory,processes=2)
        reduced,audit=adaptive_spectrum(fine,meta['radius_m'],tolerance=.01)
        shield,_=sunlight('design')
        transmission=shield(fine['wavelength'])
        if np.any(transmission<=0):
            raise ValueError('Cannot recover unfiltered source where shield transmission is zero')
        bands=np.unique(reduced['band'],axis=0)
        reduced['unshielded_bands']=bands
        reduced['unshielded_energy']=np.array([
            np.sum((fine['energy']/transmission)[(fine['wavelength']>=lo-1e-9)&
                   ((fine['wavelength']<hi) if hi<830 else (fine['wavelength']<=hi+1e-8))])
            for lo,hi in bands])
        meta['spectral_reduction']=audit
        meta['calendar_optics_producer_sha256']=digest(Path(__file__))
        np.savez_compressed(optical_file,**reduced,meta=json.dumps(meta))
        del fine
    solution_path=directory/'point_standard.npz'
    identity=transport_identity(optical_file)
    identity_file=directory/'point_inputs.json'
    if solution_path.exists():
        if not identity_file.exists() or json.loads(identity_file.read_text())!=identity:
            raise ValueError('Point-source cache identity differs; use a fresh --directory')
        with np.load(solution_path,allow_pickle=False) as z:
            if not json.loads(str(z['meta']))['order_converged']:
                raise ValueError('Cached point-source transport did not converge')
        print('Reusing completed point-source solution',flush=True)
        return
    # Tighten the retained orders because the calculator reaches the antisolar point.
    solution,settings=solve(reduced,meta['radius_m'],'standard',max_orders=240,tolerance=1e-6,path=solution_path)
    if not settings['order_converged']:
        raise RuntimeError('Point-source transport has unconverged channels')
    identity_file.write_text(json.dumps(identity,indent=2)+'\n')


def export(directory=DIRECTORY, product=PRODUCT):
    directory=Path(directory)
    if json.loads((directory/'point_inputs.json').read_text())!=transport_identity(directory/'optics.npz'):
        raise ValueError('Point-source producer/input identity differs')
    with np.load(directory/'point_standard.npz',allow_pickle=False) as z:
        settings=json.loads(str(z['meta']))
        if not settings['order_converged'] or settings['solar_disk_radius_deg']!=0:
            raise ValueError('A converged point-source transport is required')
        angles=np.degrees(z['a']);diffuse=z['moments'][0,:,:,5].copy()
        last=z['last_order'][0,:,:,5].copy();xyz=z['xyz'].copy();energy=z['energy'].copy();channel_bands=z['band'].copy()
        # Fine direct-beam sampling resolves grazing rays before disk integration.
        beam_angles=np.unique(np.r_[-90.,0.,np.arange(.01,1.001,.01),np.arange(1.05,10.01,.05),np.arange(10.25,90.01,.25)])
        direct=point_table(z['r'][:1],np.radians(beam_angles),z['edges'],z['scattering']+z['absorption'])[0,:,:,1]
    with np.load(directory/'optics.npz',allow_pickle=False) as z:
        bands=z['unshielded_bands'];unfiltered=z['unshielded_energy']
        optical_meta=json.loads(str(z['meta']))
        check_optics_provenance(optical_meta)
    shielded=np.array([energy[(channel_bands==b).all(1)].sum() for b in bands])
    index=np.array([int(np.flatnonzero((bands==b).all(1))[0]) for b in channel_bands])
    phases=np.arange(0.,180.001,2.)
    earth_irradiance=np.array([earthlight.calibrated_irradiance(bands,unfiltered,a,EARTH_MOON_DISTANCE) for a in phases])
    weights=(earth_irradiance/shielded)[:,index]
    photopic=weights*xyz[:,1]
    solar_diffuse=diffuse@xyz[:,1];solar_direct=direct@xyz[:,1]
    earth_diffuse=photopic@diffuse.T;earth_direct=photopic@direct.T
    # The browser receives absolute illuminance, preserving zero-valued direct beams.
    files=['illumination/calendar/build.py','illumination/calendar/geometry.py','illumination/earthlight/model.py',
           'illumination/earthlight/inputs.json','illumination/sky/solved_transport.py','illumination/sky/solver.py',
           'geography/lunar_ephemeris.py','illumination/calendar/disk.py','illumination/calendar/evaluator.mjs',
           'illumination/sky/solved_optics.py','illumination/calendar/transport.py']
    data=dict(schema=SCHEMA,units='lux on an unobstructed horizontal surface',
              evidence='Scalar spherical multiple scattering on the solved 1.2-atm column; finite sources integrated over a point-source response.',
              reading_rule='Linear interpolation in native source-elevation nodes and phase. Disk integration follows evaluator.mjs. No stellar/airglow floor is added.',
              producer=dict(files={p:digest(ROOT/p) for p in files},constants=constants_used(files)),
              inputs={p.name:digest(p) for p in [directory/'point_standard.npz',directory/'optics.npz']},
              optical_producer=optical_meta['producer'],spectral_reduction=optical_meta['spectral_reduction'],
              transport_settings=settings,
              constants=dict(au_m=AU,earth_radius_m=EARTH_RADIUS,moon_radius_m=MOON_RADIUS,sun_radius_m=SUN_RADIUS,
                             earth_reference_distance_m=EARTH_MOON_DISTANCE),
              atmosphere=dict(scenario='moon_1.2atm',scenario_pressure_atm=1.2,
                              column=optical_meta['column'],ground_albedo=.1,shield='Titania stack × 0.95',
                              geometry='Spherical, horizontally uniform molecular column; straight rays'),
              source_model=dict(earth_spectrum='Glenar et al. 2019; visual normalization Robinson et al. 2025',
                                phase_bounds='Spectral fit derived below 60 degrees and extrapolated; spectral shape held beyond 144 degrees. Visual observations span 5 to 144 degrees; a continuous Lambert-shaped tail reaches new Earth.',
                                disk='Lambert incidence pattern normalized to measured disk-integrated phase brightness; spatial spectrum uniform',
                                finite_distance='Topocentric directions, angular radii and distances; the atmosphere response to each disk element uses parallel rays',
                                variation='Mean Earth spectrum; continents, seasons and cloud changes are not predicted'),
              serialization_significant_digits=10,
              diffuse_elevation_deg=numbers(angles),direct_elevation_deg=numbers(beam_angles),phase_deg=numbers(phases),
              solar=dict(direct_lux=numbers(solar_direct),diffuse_lux=numbers(solar_diffuse)),
              earth=dict(direct_lux=numbers(earth_direct),diffuse_lux=numbers(earth_diffuse),above_air_lux=numbers(photopic.sum(1)),
                         spectrum=dict(band_edges_nm=numbers(bands),solar_irradiance_w_m2=numbers(unfiltered),
                                       earth_irradiance_w_m2_by_phase=numbers(earth_irradiance),
                                       apparent_disk_reflectance_by_phase=numbers(earth_irradiance/unfiltered/(EARTH_RADIUS/EARTH_MOON_DISTANCE)**2),
                                       reading_rule='Band-integrated flux at 1 AU solar distance and the reference Earth/Moon distance. Apparent disk reflectance is Earth irradiance divided by incident solar band irradiance and (Earth radius/reference distance)^2; phase zero is the geometric-albedo spectrum.')),
              checks=dict(last_order_fraction_max=float(np.max((last@xyz[:,1])/np.maximum(solar_diffuse,1e-300))),
                          scattering_orders=max(settings['channel_orders'])),
              sources=[dict(title='Glenar et al., Earthshine as an illumination source at the Moon',url='https://doi.org/10.1016/j.icarus.2018.12.025'),
                       dict(title='Earthshine spectrum table, CC BY 4.0',url='https://doi.org/10.17632/xfjm6nmh3m.2'),
                       dict(title='Earthshine data licence, CC BY 4.0',url='https://creativecommons.org/licenses/by/4.0/'),
                       dict(title='Robinson et al. (2025), Earth phase curve',url='https://arxiv.org/abs/2507.22258')])
    product=Path(product);product.parent.mkdir(parents=True,exist_ok=True)
    product.write_text(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')
    print('Exported',product,product.stat().st_size,'bytes',data['checks'],flush=True)
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compute',action='store_true');p.add_argument('--export',action='store_true')
    p.add_argument('--directory',type=Path,default=DIRECTORY)
    args=p.parse_args()
    if args.compute:compute(args.directory)
    if args.export:export(args.directory)
    if not(args.compute or args.export):p.error('Choose --compute or --export')


if __name__=='__main__':main()
