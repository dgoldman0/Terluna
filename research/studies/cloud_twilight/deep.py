"""Search long sunsetward sightlines for the latest illuminated cloud candidates."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from climate.crm.cloud_columns import sha256
from climate.crm.cloud_scene import unit_vectors
from shared.constants import MOON_RADIUS,SYNODIC_MONTH_DAYS

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
RUNS=ROOT/'research/runs/optical_comfort/evening'


def main():
    path=RUNS/'regional_inputs.json';manifest=json.loads(path.read_text())
    sky_path=ROOT/'illumination/sky/results/solved_sky.json';sky=json.loads(sky_path.read_text())
    model_path=ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz'
    if sha256(model_path)!=sky['additional_checks']['numerical_checks.json']['inputs'][model_path.name]:
        raise ValueError('Molecular archive differs')
    context_path=ROOT/'illumination/cloud_light/results/evening_context.json'
    context=json.loads(context_path.read_text())
    with np.load(model_path) as m:
        light=RegularGridInterpolator((m['sr']-MOON_RADIUS,np.degrees(m['sa'])),m['beam'][:,:,:,0]@m['xyz'][:,1])
    candidates=[];counts=[]
    offsets=(-96,-80,-64,-48,-32,32,48,64,80,96)
    for entry in manifest['products'][:3]:
        archive=ROOT/entry['path']
        if sha256(archive)!=entry['sha256']:
            raise ValueError('Cloud archive differs')
        with np.load(archive) as f:
            d={k:f[k] for k in f.files if k!='metadata'};meta=json.loads(str(f['metadata']))
        lat,lon=d['latitude_deg'],d['longitude_deg'];points=unit_vectors(lat,lon)
        alt=d['sun_elevation_deg'];h=d['hour_angle_deg'];n=len(lat)
        for offset in offsets:
            target=(np.arange(n)+offset)%n;top=d['cloud_top_m'][:,target]
            separation=np.sum(points*points[target],axis=1)
            cr=MOON_RADIUS+top;obs=MOON_RADIUS+1.6
            distance=np.sqrt(np.maximum(1,cr*cr+obs*obs-2*cr*obs*separation))
            mu=(cr*separation-obs)/distance
            elevation=np.degrees(np.arcsin(np.clip(mu,-1,1)))
            incident=light(np.stack((top.ravel(),alt[:,target].ravel()),axis=-1)).reshape(top.shape)
            ground=np.exp(np.interp(alt,context['sun_deg'],np.log(context['surface_lux'])))
            mass=d['cloud_path_above_20000_kg_m2'][:,target]
            possible=(h>=90)&(alt<=-25)&(alt>=-45)&(top>=20000)&(elevation>=.5)&(elevation<=15)
            possible &= d['land'][None,:].astype(bool)&(d['low_cloud']==0)&(d['qc_path_kg_m2']<.005)
            score=incident*np.maximum(mu,0)*(1-np.exp(-100*mass))/np.sqrt(np.maximum(ground,1))
            for lo,hi in [(25,30),(30,35),(35,40),(40,45)]:
                window=possible&(-alt>=lo)&(-alt<hi)
                counts.append(dict(case=entry['case'],offset_columns=offset,sun_depression_deg=[lo,hi],
                                   geometrically_eligible_pairs=int(window.sum()),
                                   pairs_above_0p1lux_incident=int(np.sum(window&(incident>=.1))),
                                   pairs_above_1lux_incident=int(np.sum(window&(incident>=1))),
                                   pairs_above_10lux_incident=int(np.sum(window&(incident>=10)))))
                indices=np.argpartition(np.where(window&(incident>=.1),score,0).ravel(),-4)[-4:]
                for flat in indices:
                    t,i=np.unravel_index(flat,score.shape)
                    if not window[t,i] or incident[t,i]<.1:
                        continue
                    j=int(target[i])
                    candidates.append(dict(case=entry['case'],snapshot=int(d['snapshot'][t]),observer_column=int(i),cloud_column=j,
                        observer_latitude_deg=float(lat[i]),observer_longitude_deg=float(lon[i]),observer_land=True,
                        cloud_latitude_deg=float(lat[j]),cloud_longitude_deg=float(lon[j]),time_day=float(d['time_day'][t]),
                        sun_elevation_deg=float(alt[t,i]),hours_after_sunset=float((h[t,i]-90)*SYNODIC_MONTH_DAYS*24/360),
                        cloud_top_km=float(top[t,i]/1000),view_elevation_deg=float(elevation[t,i]),range_km=float(distance[t,i]/1000),
                        cloud_direct_normal_lux=float(incident[t,i]),clear_ground_lux=float(ground[t,i]),
                        upper_cloud_water_path_kg_m2=float(mass[t,i]),selection_score=float(score[t,i]),
                        source_sha256=next(r['sha256'] for r in meta['raw_inputs'] if r['snapshot']==int(d['snapshot'][t]))))
    candidates.sort(key=lambda c:c['selection_score'],reverse=True)
    chosen=[]
    for lo,hi in [(25,30),(30,35),(35,40),(40,45)]:
        pool=[c for c in candidates if lo<=-c['sun_elevation_deg']<hi]
        if pool:
            chosen.append(dict(selection=f'long-range source candidate, Sun {lo}–{hi} degrees down',**pool[0]))
    product=dict(schema='terluna.research.deep-evening-candidates/1',
        producer=dict(file=str(Path(__file__).relative_to(ROOT)),sha256=sha256(__file__),
                      inputs={str(p.relative_to(ROOT)):sha256(p) for p in (path,sky_path,context_path)}),
        evidence='Expanded geometric and incident-beam search on the same saved second-cycle tracks. '
                 'Sun depression 25–45 degrees, elevations 0.5–15 degrees, target offsets to 96 columns (about 577 km). '
                 'Pair counts share observers and clouds and serve as search diagnostics. Foreground extinction and cloud radiance are further calculations.',
        selected_scenes=chosen,candidates=candidates[:80],pair_counts=counts)
    (HERE/'results/deep_evenings.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(json.dumps(chosen,indent=2))


if __name__=='__main__':
    main()
