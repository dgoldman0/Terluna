"""Clear evening light, darkness thresholds and elevated-cloud geometry.

python -m illumination.cloud_light.evening_context
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

from illumination.sky.solved_sky import load, horizontal_integral, digest, photometry
from illumination.sky.solved_transport import evaluate, SUN_RADIUS_RAD
from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT/'research/runs/optical_comfort/evening'
SCHEMA = 'terluna.illumination.evening-context/1'
SUNS = np.r_[np.arange(-60.,-24.,2.),np.arange(-24.,1.,1.),2.,5.,10.]


def last_sunlit_depression(height_m, minimum_elevation_deg=0., observer_height_m=1.6):
    """Most favourable sunsetward geometry, including the Sun's last limb."""
    e = np.radians(minimum_elevation_deg)
    observer = MOON_RADIUS+observer_height_m
    cloud = MOON_RADIUS+height_m
    separation = np.arccos(np.clip(observer/cloud*np.cos(e),-1,1))-e
    return np.degrees(separation+np.arccos(MOON_RADIUS/cloud)+SUN_RADIUS_RAD)


def threshold_depression(suns,illuminance,level):
    if level<illuminance.min() or level>illuminance.max():
        return None
    return float(-np.interp(np.log(level),np.log(illuminance),suns))


def main():
    RUNS.mkdir(parents=True,exist_ok=True)
    source = ROOT/'illumination/sky/results/solved_sky.json'
    p = json.loads(source.read_text())
    folder = ROOT/'research/runs/optical_comfort/spherical'
    expected = p['additional_checks']['numerical_checks.json']['inputs']
    flux,arcs,settings = {},{},{}
    elevations = 90*np.linspace(0,1,41)**2
    azimuths = np.linspace(0,180,65)
    files = {}
    for quality in ('draft','standard'):
        path = folder/f'moon_1.2atm_{quality}.npz'
        if digest(path)!=expected[path.name]:
            raise ValueError('Molecular field hash differs')
        s,m = load(path)
        spectrum = evaluate(s,m,SUNS,elevations,azimuths)
        xyz = spectrum@s['xyz']
        values = np.array([horizontal_integral(xyz[:,:,:,c],elevations,azimuths) for c in range(3)]).T
        _,direct,_ = photometry(s,SUNS)
        flux[quality] = values[:,1]+direct[:,1]
        arcs[quality] = xyz
        settings[quality] = m
        files[str(path.relative_to(ROOT))] = digest(path)
        print(quality, 'deepest flux',flux[quality][0],flush=True)
    out = RUNS/'clear_evening.npz'
    meta = dict(schema=SCHEMA, source_sha256=digest(source), files=files,
                units=dict(xyz='photopic CIE XYZ; Y radiance in cd m^-2',surface_illuminance='lux'),
                reading_rule='Bilinear angular interpolation in sin(elevation), azimuth; linear solar-angle radiance. '
                             'Solar illumination only. Molecular air and ground albedo 0.1 match solved_sky.json. '
                             'Deep-twilight values extend the same scattering series; refinement is recorded separately.')
    np.savez_compressed(out,metadata=json.dumps(meta),suns=SUNS,elevations=elevations,azimuths=azimuths,
                        xyz=arcs['standard'],total_horizontal_lux=flux['standard'],
                        coarse_total_horizontal_lux=flux['draft'])
    thresholds = [dict(level_lux=v,sun_depression_deg=threshold_depression(SUNS,flux['standard'],v),
                       coarse_sun_depression_deg=threshold_depression(SUNS,flux['draft'],v)) for v in (100.,10.,1.,.1)]
    bounds = []
    for height in (20.,24.,36.,40.,60.,80.,100.):
        for elevation in (0.,5.,10.):
            dep = last_sunlit_depression(height*1000,elevation)
            light = float(np.exp(np.interp(-dep,SUNS,np.log(flux['standard']))))
            bounds.append(dict(cloud_height_km=height,minimum_view_elevation_deg=elevation,
                               last_geometric_sun_depression_deg=dep,clear_ground_lux=light))
    rows = []
    for lat in (0.,15.,30.,45.,60.,70.,80.):
        row = dict(latitude_deg=lat,darkest_sun_elevation_deg=lat-90,
                   midnight_lux=(float(np.exp(np.interp(lat-90,SUNS,np.log(flux['standard'])))) if lat>=30 else None))
        for t in thresholds:
            d = t['sun_depression_deg']
            if d is None or d>90-lat:
                hours = None
            else:
                angle = np.degrees(np.arcsin(np.sin(np.radians(d))/np.cos(np.radians(lat))))
                hours = float(angle*SYNODIC_MONTH_DAYS*24/360)
            row[f'hours_after_sunset_to_{t["level_lux"]:g}_lux'] = hours
        rows.append(row)
    result = dict(schema=SCHEMA,producer=dict(files={str(Path(__file__).relative_to(ROOT)):digest(__file__)},inputs=files),
                  source=dict(path=str(source.relative_to(ROOT)),sha256=digest(source)),
                  archive=dict(path=str(out.relative_to(ROOT)),sha256=digest(out)),
                  assumptions=meta['reading_rule'],sun_deg=SUNS.tolist(),surface_lux=flux['standard'].tolist(),
                  coarse_relative_difference=(flux['draft']/flux['standard']-1).tolist(),
                  thresholds=thresholds,geometry_bounds=bounds,latitude_comparison=rows,
                  evidence='Clear solar twilight on the solved mean column. Geometric bounds allow a cloud at any '
                           'sunsetward distance and exclude extinction from the illumination criterion. '
                           'Cloud radiance, terrain and earthlight are separate calculations.')
    target=ROOT/'illumination/cloud_light/results/evening_context.json'
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(thresholds,indent=2),flush=True)


if __name__=='__main__':
    main()
