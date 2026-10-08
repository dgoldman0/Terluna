"""Export ephemeris coefficients, constants, sites, map and Python golden geometry."""
import json
from pathlib import Path
import numpy as np
from geography import lunar_ephemeris as ephemeris
from illumination.calendar.geometry import geometry
import hashlib
from shared.constants import AU, EARTH_RADIUS, MOON_RADIUS, SUN_RADIUS
from shared.provenance import constants_used

ROOT=Path(__file__).resolve().parents[2]
DEST=Path(__file__).parent/'results'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build():
    source=ROOT/'research/studies/sea_appearance/results/lighting_calendar.json'
    calendar=json.loads(source.read_text())
    sites=[dict(id=str(i),name=s['short'],description=s['name'],longitude=s['lon_lat'][0],latitude=s['lon_lat'][1])
           for i,s in enumerate(calendar['sites'])]
    sites += [dict(id='north',name='Northern high latitudes',description='A clear-sky view at 80° north',longitude=0,latitude=80),
              dict(id='farside',name='Far-side equator',description='The deep-night reference on the far side',longitude=180,latitude=0)]
    with np.load(ROOT/'research/studies/sea_appearance/results/earth_over_seas.npz',allow_pickle=False) as z:
        mask=np.zeros((180,360),dtype=bool);mask[z['rows'],z['cols']]=True
        map_grid=dict(latitude_centres=z['lat_deg'].tolist(),longitude_centres=z['lon_deg'].tolist())
    spans=[]
    for i,row in enumerate(mask):
        changes=np.diff(np.r_[0,row.astype(int),0]);starts=np.flatnonzero(changes==1);ends=np.flatnonzero(changes==-1)
        spans += [[int(i),int(a),int(b)] for a,b in zip(starts,ends)]
    checks=[]
    for date in [2451544.5,2461041.5,2465476.3704272127,2488069.125,2547339.5,2634531.4791666665]:
        for lon,lat in [(0,0),(93.4,2.34),(-73.875,25.125),(180,-33.7),(0,89.9)]:
            values=geometry(date,lon,lat)
            checks.append(dict(jd_tt=date,longitude=lon,latitude=lat,values={k:float(v) for k,v in values.items()}))
    validation=ROOT/'geography/calendar_ephemeris_check.json'
    files=[Path(__file__),ROOT/'geography/lunar_ephemeris.py',Path(__file__).with_name('geometry.py')]
    data=dict(schema='terluna.illumination.calendar-astronomy/1',
              producer=dict(files={str(p.relative_to(ROOT)):digest(p) for p in files},constants=constants_used(files)),
              inputs={str(p.relative_to(ROOT)):digest(p) for p in [source,validation,
                      ROOT/'research/studies/sea_appearance/results/earth_over_seas.npz']},
              evidence='Portable coefficients of the date ephemeris with sampled JPL comparisons; scenario coasts and water mask from the sea-appearance products.',
              reading_rule='Evaluate geometry with evaluator.mjs, using TT and degrees east/north. Read calendar_range as inclusive years; map_grid holds the water-mask cell centres.',
              time_scale='TT for the orbital model; UTC uses historical leap seconds through 2017 and a fixed 69.184-second TT−UTC assumption thereafter',
              calendar_range=[2000,2500],longitude_distance=ephemeris.LONGITUDE_DISTANCE.tolist(),
              latitude=ephemeris.LATITUDE.tolist(),equator_to_ecliptic=ephemeris.EQUATOR_TO_ECLIPTIC,
              au_m=AU,earth_radius_m=EARTH_RADIUS,moon_radius_m=MOON_RADIUS,sun_radius_m=SUN_RADIUS,
              sites=sites,water_spans=spans,map_grid=map_grid,golden_geometry=checks,
              ephemeris_validation=json.loads(validation.read_text()) if validation.exists() else None)
    DEST.mkdir(parents=True,exist_ok=True)
    (DEST/'astronomy.json').write_text(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')
    print('Exported portable astronomy',len(checks),'geometry cases')
    return data


if __name__=='__main__':build()
