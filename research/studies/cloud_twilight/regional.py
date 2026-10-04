"""Regional cloud occurrence and reproducible scene selection for lunar evenings."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
RUNS=ROOT/'research/runs/optical_comfort/evening'
SCHEMA='terluna.research.regional-cloud-evenings/1'


def digest(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit_vectors(lat,lon):
    a,b=np.radians(lat),np.radians(lon)
    return np.stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)),axis=-1)


def describe(top,mask):
    values=top[mask]/1000
    clouds=values[values>0]
    return dict(samples=int(mask.sum()),cloud_fraction=float(np.mean(values>0)) if len(values) else None,
                cloud_top_km_quantiles=np.percentile(clouds,[50,90,99,100]).tolist() if len(clouds) else None,
                fraction_above_40km=float(np.mean(values>=40)) if len(values) else None)


def build():
    from shared.constants import MOON_RADIUS,SYNODIC_MONTH_DAYS
    manifest=RUNS/'regional_inputs.json'
    inputs=json.loads(manifest.read_text())
    context_path=ROOT/'illumination/cloud_light/results/evening_context.json'
    context=json.loads(context_path.read_text())
    sky_path=ROOT/'illumination/sky/results/solved_sky.json'
    sky=json.loads(sky_path.read_text())
    model_path=ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz'
    if digest(model_path)!=sky['additional_checks']['numerical_checks.json']['inputs'][model_path.name]:
        raise ValueError('Sky moment archive differs')
    with np.load(model_path) as m:
        beam=m['beam'][:,:,:,0]@m['xyz'][:,1]
        beam_lookup=RegularGridInterpolator((m['sr']-MOON_RADIUS,np.degrees(m['sa'])),beam,bounds_error=True)
    summaries,regions,candidates=[],[],[]
    for record in inputs['products']:
        path=ROOT/record['path']
        if digest(path)!=record['sha256']:
            raise ValueError('Regional cloud archive differs')
        with np.load(path,allow_pickle=False) as raw:
            d={k:raw[k] for k in raw.files if k!='metadata'}
            meta=json.loads(str(raw['metadata']))
        lat,lon=d['latitude_deg'],d['longitude_deg']
        hour,alt,top=d['hour_angle_deg'],d['sun_elevation_deg'],d['cloud_top_m']
        evening=(hour>=90)&(hour<135)
        per_band=[]
        for lo,hi in ((0,15),(15,30),(30,45),(45,60),(60,75)):
            m=evening&(abs(lat[None,:])>=lo)&(abs(lat[None,:])<hi)
            if m.any():
                per_band.append(dict(abs_latitude_deg=[lo,hi],**describe(top,m)))
        by_time=[]
        for lo,hi in ((0,12),(12,24),(24,48),(48,72),(72,120)):
            time=(hour-90)*SYNODIC_MONTH_DAYS*24/360
            m=(time>=lo)&(time<hi)&(hour>=90)
            by_time.append(dict(hours_after_sunset=[lo,hi],**describe(top,m)))
        summaries.append(dict(case=record['case'],span_days=record['span_days'],
                              columns=record['columns'],latitude_bands=per_band,time_bins=by_time,
                              overall=describe(top,evening),geometry=meta['geometry']))
        if meta.get('path') is None:
            continue
        # Regional bins retain only sampled track columns, with a recorded denominator.
        latkey=np.floor((lat+90)/10).astype(int)
        lonkey=np.floor((lon%360)/15).astype(int)
        keys=latkey*24+lonkey
        for key in np.unique(keys):
            columns=keys==key
            mask=evening&columns[None,:]
            regional=dict(case=record['case'],latitude_deg=float(lat[columns].mean()),
                          longitude_deg=float(lon[columns].mean()),columns=int(columns.sum()),
                          land_fraction=float(d['land'][columns].mean()),**describe(top,mask))
            regions.append(regional)
        # Fast geometric screen of observers 48, 96 or 192 km along either side of a cloud.
        points=unit_vectors(lat,lon)
        n=len(lat)
        for offset in (-32,-16,-8,8,16,32):
            target=(np.arange(n)+offset)%n
            height=top[:,target]
            cossep=np.sum(points*points[target],axis=1)
            cloud_r=MOON_RADIUS+height
            obs_r=MOON_RADIUS+1.6
            distance=np.sqrt(np.maximum(1,cloud_r**2+obs_r**2-2*cloud_r*obs_r*cossep))
            mu=(cloud_r*cossep-obs_r)/distance
            elevation=np.degrees(np.arcsin(np.clip(mu,-1,1)))
            cloud_sun=alt[:,target]
            light=beam_lookup(np.stack((height.ravel(),cloud_sun.ravel()),axis=-1)).reshape(height.shape)
            clear_ground=np.exp(np.interp(alt,context['sun_deg'],np.log(context['surface_lux'])))
            water=d['cloud_path_above_20000_kg_m2'][:,target]
            eligible=(hour>=90)&(hour<175)&(alt<0)&(alt>-40)&(height>=20000)&(elevation>=5)&(light>=10)
            eligible &= (d['low_cloud']==0)&(d['qc_path_kg_m2']<.005)
            # A selection score: incident light, projected view and available upper cloud mass.
            # It supplies candidate scenes; cloud radiance is computed separately.
            score=light*np.maximum(mu,0)*(1-np.exp(-100*water))/np.sqrt(np.maximum(clear_ground,1))
            score=np.where(eligible,score,0)
            for low,high,cap in ((0,10,150000),(10,20,150000),(20,30,150000),(30,40,150000),(5,15,40000)):
                window=(-alt>=low)&(-alt<high)&(height<=cap)
                selected=np.argpartition((score*window).ravel(),-8)[-8:]
                for flat in selected:
                    t,i=np.unravel_index(flat,score.shape)
                    if score[t,i]<=0 or not window[t,i]:
                        continue
                    j=int(target[i])
                    candidates.append(dict(case=record['case'],snapshot=int(d['snapshot'][t]),time_day=float(d['time_day'][t]),
                         observer_column=int(i),cloud_column=j,observer_latitude_deg=float(lat[i]),observer_longitude_deg=float(lon[i]),
                         cloud_latitude_deg=float(lat[j]),cloud_longitude_deg=float(lon[j]),observer_land=bool(d['land'][i]),
                         sun_elevation_deg=float(alt[t,i]),hours_after_sunset=float((hour[t,i]-90)*SYNODIC_MONTH_DAYS*24/360),
                         cloud_top_km=float(height[t,i]/1000),view_elevation_deg=float(elevation[t,i]),range_km=float(distance[t,i]/1000),
                         cloud_direct_normal_lux=float(light[t,i]),clear_ground_lux=float(clear_ground[t,i]),
                         upper_cloud_water_path_kg_m2=float(water[t,i]),selection_score=float(score[t,i]),
                         source_sha256=next(r['sha256'] for r in meta['raw_inputs'] if r['snapshot']==int(d['snapshot'][t]))))
    candidates.sort(key=lambda r:r['selection_score'],reverse=True)
    chosen=[]
    # Strong early and deep examples, plus a different regional ring where available.
    for lo,hi in ((0,10),(10,20),(20,30),(30,40)):
        pool=[r for r in candidates if lo<=-r['sun_elevation_deg']<hi and r['observer_land']]
        if pool:
            chosen.append(dict(selection=f'largest land-observer source score at Sun depression {lo}–{hi}°',**pool[0]))
    pool=[r for r in candidates if 5<=-r['sun_elevation_deg']<15 and r['cloud_top_km']<=40 and r['observer_land']]
    if pool:
        chosen.append(dict(selection='strong 20–40 km cloud example in established twilight',**pool[0]))
    for case in ('ring_70_45e','ring_70_135e'):
        pool=[r for r in candidates if r['case']==case and r['observer_land']]
        if pool:
            chosen.append(dict(selection=f'largest source score on {case}',**pool[0]))
    sources={str(p.relative_to(ROOT)):digest(p) for p in (manifest,context_path,sky_path)}
    product=dict(schema=SCHEMA,producer=dict(file=str(Path(__file__).relative_to(ROOT)),sha256=digest(__file__),inputs=sources),
        evidence='One saved post-spin-up lunar cycle on three corrected flat great-circle rings, plus a corrected '
                 'terrain box. Regional occurrence describes the sampled tracks. Scene selection uses geometric '
                 'exposure and incident illumination; full cloud appearance is a separate transport product.',
        reading_rule='Cloud counts use qc+qi >= 1e-5 kg/kg; each sampled column has equal weight. '
                     'Regional bins span 10 degrees latitude and 15 degrees longitude and contain only track samples. '
                     'Times are hours after equinox geometric solar-centre sunset. Sea-level clear light uses a common '
                     'solved molecular column. Ranking scores select examples from one cycle; repeatability requires further cycles.',
        climate_inputs=inputs['products'],case_summaries=summaries,regional_bins=regions,
        selected_scenes=chosen,candidates=candidates[:120],darkness_context=context,
        remaining_inputs=['Three-dimensional cloud structure beyond the highland box','Regional molecular optical profiles',
                          'Ice habits and measured phase functions','Terrain horizons and earthlight for individual sites'])
    output=HERE/'results/regional_evenings.json'
    output.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(json.dumps(chosen,indent=2),flush=True)
    return product


if __name__=='__main__':
    build()
