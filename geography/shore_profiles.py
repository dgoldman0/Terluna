"""Extract explicit cross-shore rock profiles from the pinned LOLA terrain."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from shared.constants import MOON_RADIUS

ROOT = Path(__file__).resolve().parents[1]
TERRAIN = ROOT/"research/runs/waves/shore/terrain.npz"
STATIONS = ROOT/"climate/waves/results/shore.json"
OUTPUT = ROOT/"research/runs/waves/shore_history/profiles.json"


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def extract(terrain_path=TERRAIN, stations_path=STATIONS):
    with np.load(terrain_path,allow_pickle=False) as data:
        meta=json.loads(data['metadata'].item())
        if meta['schema']!='terluna.geography.coastal-grid/1':
            raise ValueError('Unsupported shore terrain')
        x,y=data['longitude_deg'],data['latitude_deg']
        depth=meta['sea_level_m']-data['height_m']
    stations=json.loads(stations_path.read_text())
    if stations['schema']!='terluna.climate.shore-exposure/1':
        raise ValueError('Unsupported station selection product')
    selected=next(c for c in stations['cases'] if c['name']=='h852_s1_d36')['sites']
    interpolate=RegularGridInterpolator((y,x),depth,bounds_error=True)
    profiles={}
    for name in ('west_face','east_face','southern_shore'):
        station=selected[name]
        xy=np.array([station['longitude_deg'],station['latitude_deg']])
        i,j=np.argmin(abs(x-xy[0])),np.argmin(abs(y-xy[1]))
        np.testing.assert_allclose([x[i],y[j]],xy,atol=1e-9,rtol=0)
        metres=MOON_RADIUS*np.pi/180*np.array([np.cos(np.deg2rad(xy[1])),1.])
        gradient=np.array([(depth[j,i+1]-depth[j,i-1])/((x[i+1]-x[i-1])*metres[0]),
                           (depth[j+1,i]-depth[j-1,i])/((y[j+1]-y[j-1])*metres[1])])
        normal=-gradient/np.linalg.norm(gradient)
        distance=np.arange(0.,1000.01,2.)
        coordinates=xy+distance[:,None]*normal/metres
        values=interpolate(coordinates[:,::-1])
        high=np.flatnonzero(values<=-15.)
        if not len(high) or values[0]<=0:
            raise ValueError('Profile must connect the water station to land 15 m above sea level')
        stop=high[0]+1
        distance,values,coordinates=distance[:stop],values[:stop],coordinates[:stop]
        shore=np.flatnonzero(values<=0)[0]
        shore_distance=float(np.interp(0.,-values[shore-1:shore+1],distance[shore-1:shore+1]))
        profiles[name]=dict(station_longitude_latitude=xy.tolist(),normal_east_north=normal.tolist(),
            station_depth_m=float(values[0]),shore_distance_m=shore_distance,
            mean_submerged_slope=float(values[0]/shore_distance),distance_m=distance.tolist(),
            depth_m=values.tolist(),longitude_latitude=coordinates.tolist())
    return dict(schema='terluna.geography.shore-profiles/1',
        producer=dict(path='geography/shore_profiles.py',sha256=digest(Path(__file__))),
        terrain_sha256=digest(terrain_path),station_product_sha256=digest(stations_path),
        evidence='Cross-shore samples of the 118 m altimetry-derived flooded-rock raster.',
        reading_rule='Bilinear sampling at 2 m supplies a continuous numerical bed. Geographic information retains the native raster spacing. Sediment profiles remain explicit alternative geometries.',
        units=dict(distance='m from each water station toward land',depth='m below fixed sea level',normal='unit east/north vector'),
        sea_level_m=meta['sea_level_m'],native_spacing_m=MOON_RADIUS*np.pi/180/256,
        profiles=profiles)


def numerical_profile(record,station,*,dx=.5,kind='rock',buffer_m=100.,slope_denominator=50):
    """Add an explicit constant-depth wavemaker reach and dry land."""
    if record['schema']!='terluna.geography.shore-profiles/1' or dx<=0 or buffer_m<0:
        raise ValueError('Expected a shore profile with positive numerical spacing')
    p=record['profiles'][station]
    depth=p['station_depth_m']
    if kind=='rock':
        length=buffer_m+p['distance_m'][-1]
    elif kind=='slope' and slope_denominator>0:
        length=buffer_m+(depth+15)*slope_denominator
    else:
        raise ValueError('Choose rock or a positive explicit slope')
    x=np.arange(int(np.ceil(length/dx))+1)*dx
    local=np.maximum(x-buffer_m,0)
    bottom=(np.interp(local,p['distance_m'],p['depth_m']) if kind=='rock'
            else depth-local/slope_denominator)
    return x,bottom,dict(kind=kind,dx_m=dx,constant_depth_buffer_m=buffer_m,
                        slope_denominator=slope_denominator if kind=='slope' else None,
                        station=station,station_depth_m=depth,
                        still_water_shore_m=buffer_m+(p['shore_distance_m'] if kind=='rock' else depth*slope_denominator))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(extract(),indent=2)+'\n')
    print(args.output)
