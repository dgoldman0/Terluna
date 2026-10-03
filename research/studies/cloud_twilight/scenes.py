"""Run selected regional cloud scenes and fixed-observer evening histories.

python -m research.studies.cloud_twilight.scenes --mode pilot
"""
from __future__ import annotations
import argparse,json,fcntl
from pathlib import Path
import numpy as np
from climate.crm.cloud_scene import export
from illumination.cloud_light.render import render
from climate.crm.cloud_columns import sha256

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
RUNS=ROOT/'research/runs/optical_comfort/evening'
SCHEMA='terluna.research.cloud-evening-scenes/1'


def scene_key(pick,photons,shape,delta,width,ice_g,height):
    return (f"{pick['case']}_{pick['observer_column']}_{pick['snapshot']+delta}_"
            f"w{width:g}_g{ice_g:g}_h{height:g}_{shape[0]}x{shape[1]}_{photons}")


def one(pick,photons,shape,delta=0,width=200.,ice_g=.8,height=1.6):
    key=scene_key(pick,photons,shape,delta,width,ice_g,height)
    folder=RUNS/'scenes';folder.mkdir(parents=True,exist_ok=True)
    with (folder/f'{key}.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _one(pick,photons,shape,delta,width,ice_g,height)


def _one(pick,photons,shape,delta=0,width=200.,ice_g=.8,height=1.6):
    name=pick['case'];n=pick['snapshot']+delta;i=pick['observer_column']
    archive=RUNS/f'{name}_columns.npz'
    with np.load(archive) as a:
        header=json.loads(str(a['metadata']))
    records=[r for r in header['raw_inputs'] if r['snapshot']==n]
    if not records:
        raise ValueError('Chosen frame lies outside admitted second-cycle snapshots')
    key=scene_key(pick,photons,shape,delta,width,ice_g,height)
    section=RUNS/'sections'/f'{name}_{i}_{n}.npz'
    admission=export(name,n,i,records[0]['sha256'],section)
    offset=(pick['cloud_column']-i+header['grid']['nx']//2)%header['grid']['nx']-header['grid']['nx']//2
    az=0. if offset>0 else 180.
    completed=RUNS/'scenes'/f'{key}.json'
    if completed.exists():
        saved=json.loads(completed.read_text())
        valid=(sha256(ROOT/saved['archive']['path'])==saved['archive']['sha256'] and
               saved['source']['cloud_sha256']==admission['sha256'] and
               saved['camera_azimuth_in_ring_basis_deg']==az and
               saved['source']['molecular_product_sha256']==sha256(ROOT/'illumination/sky/results/solved_sky.json') and
               saved['source']['molecular_fields_sha256']==sha256(ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz') and
               all(sha256(ROOT/file)==value for file,value in saved['producer'].items()) and
               saved['sampled_truncated_paths']==0)
        if not valid:
            raise ValueError(f'Cached scene differs from its admitted inputs: {completed}. Use a fresh run directory or explicitly regenerate this scene.')
        saved.update(selection=pick['selection'],initial_pick=pick,delta_hours=delta*3,section=admission)
        print('Admitted existing scene',key,flush=True)
        return saved
    print('Rendering',key,flush=True)
    result=render(section,RUNS/'scenes'/f'{key}.npz',az,photons=photons,shape=shape,width_km=width,
                  ice_g=ice_g,observer_height_m=height,surface_photons=max(8192,photons*64))
    result.update(selection=pick['selection'],initial_pick=pick,delta_hours=delta*3,section=admission)
    (RUNS/'scenes'/f'{key}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Complete',key,'seconds',result['seconds'],'surface lux',result['surface_lux'],
          'recess lux',result['recess_lux'],'pixel SE',result['median_pixel_relative_standard_error'],flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['pilot','selected','history','sensitivity'],default='pilot')
    parser.add_argument('--photons',type=int,default=1024)
    args=parser.parse_args()
    product=json.loads((HERE/'results/regional_evenings.json').read_text())
    picks=product['selected_scenes']
    if args.mode=='pilot':
        result=[one(picks[3],64,(12,24))]
    elif args.mode=='selected':
        result=[one(p,args.photons,(16,32)) for p in picks]
    elif args.mode=='history':
        p=picks[1]
        result=[one(p,args.photons,(16,32),delta=t) for t in (-4,0,4,8,12,16,20)]
    else:
        p=picks[3]
        result=[one(p,args.photons,(12,24),width=w,ice_g=g,height=h)
                for w,g,h in [(100,.8,1.6),(400,.8,1.6),(200,.7,1.6),(200,.9,1.6),(200,.8,6000.)]]
    target=HERE/'results'/f'evening_scenes_{args.mode}.json'
    target.write_text(json.dumps(dict(schema=SCHEMA,mode=args.mode,
        producer=dict(file=str(Path(__file__).relative_to(ROOT)),sha256=sha256(__file__),
                      regional_input_sha256=sha256(HERE/'results/regional_evenings.json')),
        evidence='Computed optical scenarios on saved evolving cloud fields. Cross-ring width and particle phase '
                 'functions are sensitivities; the sampled regional occurrence and the rendered scenes have separate denominators.',
        scenes=result),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
