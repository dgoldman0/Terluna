"""Fixed-state drag alternatives before any coupled growth refit.

Geometry alternatives are diagnostics. The native source remains byte pinned;
the Python translation is checked against the compiled native evaluation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np

from atmosphere.snow import mcsnow_column_density as column
from atmosphere.snow import mcsnow_joint_acceptance as acceptance
from atmosphere.snow import mcsnow_drag_trace as old_trace

native, model, diagnosis = column.native, column.model, column.diagnosis
HERE=Path(__file__).resolve().parent
CACHE=column.CACHE.parent/'joint_acceptance'
DRIVER=HERE/'mcsnow_joint_drag_driver.f90'
PRODUCT=HERE/'results/mcsnow_joint_drag.json'
FIELDS=('air_density_kg_m3','viscosity_Pa_s','gravity_m_s2','mass_kg','equatorial_m',
        'polar_m','volume_m3','aspect_ratio','diameter_m','area_m2','native_speed_m_s')
MODES=('native','paper_best','span_width','paper_best_span_width','prolate_only','cylinder_only','paper_cylinder_span','paper_prism_mean','paper_prism_min','paper_prism_max')


def reynolds(phi,mass,q,rho,eta,gravity,original_best=False):
    """Pinned Reynolds_Boehm; optional originally documented q denominator."""
    k=min(max(.82+.18*phi,.85),.37+.63/phi,1.33/(max(math.log(phi),0)+1.19))
    # Two literals in the native source are default real, promoted to wp.
    gamma=max(1.,min(1.98,3.76-8.41*phi+float(np.float32(9.18))*phi**2-float(np.float32(3.53))*phi**3))
    cd=max(.292*k*gamma,.492-.2/math.sqrt(phi))*max(1.,q*(1.46*q-.46))
    denominator=max(q**.25,q) if original_best else q**.25
    best=8*mass*gravity*rho/(math.pi*eta**2*max(phi,1.)*denominator)
    cdprime=cd*(1+1.6*(best/2.8e6)**2)/(1+(best/2.8e6)**2)
    beta=math.sqrt(1+cdprime/6/k*math.sqrt(best/cdprime))-1
    re0=6*k/cdprime*beta**2
    gamma_small=(4.5*k*k*max(phi,1.)-cd)/(4*cd)
    return re0*(1+2*beta*math.exp(-beta*gamma_small)/((2+beta)*(1+beta)))


def speed(s,mode='native'):
    if mode not in MODES:
        raise ValueError('Unknown drag diagnostic')
    if mode.startswith('paper_prism_'):
        projected=dict(s)
        if s['aspect_ratio']>1:
            f={'mean':3/math.pi,'min':math.sqrt(3)/2,'max':1.}[mode.rsplit('_',1)[1]]
            projected['area_m2']=f*s['equatorial_m']*s['polar_m']
        return speed(projected,'paper_best')
    phi,m,rho,eta,g=(s[k] for k in ('aspect_ratio','mass_kg','air_density_kg_m3','viscosity_Pa_s','gravity_m_s2'))
    diameter=s['diameter_m']; area=s['area_m2']
    paper=mode.startswith('paper_')
    def velocity(q,d):
        return reynolds(phi,m,q,rho,eta,g,paper)*eta/(d*rho)
    if phi<=1:
        return velocity(area/(math.pi/4*diameter**2),diameter)
    d=(4*s['volume_m3']/math.pi/phi)**(1/3)
    if 'span' in mode:
        d=diameter/phi
    cyl=velocity(4/math.pi,d)
    pro=velocity(area/(math.pi/4*diameter**2/phi),diameter/phi)
    if mode=='prolate_only':
        return pro
    if mode in ('cylinder_only','paper_cylinder_span'):
        return cyl
    weight=math.exp(-.3*(phi-1))
    return weight*pro+(1-weight)*cyl


def states():
    result=[]
    inventory=acceptance.baseline_inventory()
    groups={}
    for r in inventory['points']:
        if r['family']!='trajectory' or r['start']=='volume':
            continue
        key=(r['temperature_C'],r['start'],r['age_min'])
        groups.setdefault(key,{})[r['quantity']]=r
    references=set()
    for (temp,start,age),rows in sorted(groups.items()):
        if (temp,age) not in references:
            references.add((temp,age))
            vals=[rows[q]['target'] for q in acceptance.QUANTITIES]
            result.append(dict(temperature_C=temp,start='reference_axes',age_min=age,values=vals,target_speed=vals[3]))
        vals=[rows[q]['value'] for q in acceptance.QUANTITIES]
        result.append(dict(temperature_C=temp,start=start,age_min=age,values=vals,
                           target_speed=rows['fall_speed_m_s']['target']))
    return result


def build():
    path=CACHE/'drag_build'
    baseline=native.verify_build(model.BASELINE)
    if path.exists():
        return verify_build()
    path.mkdir(parents=True)
    driver=path/DRIVER.name;shutil.copyfile(DRIVER,driver)
    source=model.BASELINE/'source/src';compiler=baseline['compiler_path']
    sources=native.selected_sources(source,driver,compiler)
    objects=[];commands=[]
    with (path/'build.log').open('w') as log:
        for i,src in enumerate(sources):
            obj=path/f'{i:02d}_{src.stem}.o';objects.append(str(obj))
            cmd=[compiler,*native.FLAGS,'-I'+str(source),'-I'+str(path),'-J'+str(path),'-c',str(src),'-o',str(obj)]
            commands.append(cmd)
            subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler,*native.FLAGS,*objects,'-o',str(path/'trace')];commands.append(cmd)
        subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    record=dict(schema='terluna.atmosphere.mcsnow-joint-drag-build/1',
        baseline_record_sha256=native.sha(model.BASELINE/'build.json'),source_identity=baseline['source_identity'],
        compiled_sources={str(p):native.sha(p) for p in sources},compiler_path=compiler,
        compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],flags=native.FLAGS,
        commands=commands,physics_patches=[],
        files={p.name:native.sha(p) for p in (driver,path/'trace',path/'build.log')})
    model.save(path/'build.json',record)
    return record


def verify_build():
    path=CACHE/'drag_build';r=json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE/'build.json')!=r['baseline_record_sha256']:
        raise ValueError('Baseline changed')
    for src,digest in r['compiled_sources'].items():
        if native.sha(src)!=digest:raise ValueError('Compiled source changed: '+src)
    for name,digest in r['files'].items():
        if native.sha(path/name)!=digest:raise ValueError('Build artifact changed: '+name)
    if native.sha(DRIVER)!=r['files'][DRIVER.name]:raise ValueError('Driver changed')
    return r


def run():
    build();rows=states();path=CACHE/'drag_run';path.mkdir(exist_ok=True)
    data=''.join(' '.join(map(str,[i,r['temperature_C']+273.15,101000.,*r['values'][:3]]))+'\n' for i,r in enumerate(rows))
    if (path/'input.txt').exists() and (path/'input.txt').read_text()!=data:
        raise ValueError('Existing input changed')
    (path/'input.txt').write_text(data)
    if not (path/'run.json').exists():
        with (path/'output.txt').open('w') as out,(path/'stderr.log').open('w') as log:
            p=subprocess.run([str(CACHE/'drag_build/trace'),str(path/'input.txt')],cwd=path,env=native.ENV,stdout=out,stderr=log)
        model.save(path/'run.json',dict(returncode=p.returncode,build_record_sha256=native.sha(CACHE/'drag_build/build.json'),
            files={name:native.sha(path/name) for name in ('input.txt','output.txt','stderr.log')}))
    return assemble()


def assemble():
    build_record=verify_build();path=CACHE/'drag_run';record=json.loads((path/'run.json').read_text())
    if record['returncode'] or record['build_record_sha256']!=native.sha(CACHE/'drag_build/build.json'):
        raise ValueError('Failed or altered native run')
    for name,digest in record['files'].items():
        if native.sha(path/name)!=digest:raise ValueError('Raw trace changed')
    raw=np.loadtxt(path/'output.txt');rows=states()
    if raw.shape!=(len(rows),1+len(FIELDS)) or not np.array_equal(raw[:,0],np.arange(len(rows))):
        raise ValueError('Incomplete trace')
    helper,helper_identity=old_trace.author_helper()
    helper_identity=dict(helper_identity,scope='Both archived branches are evaluated. The prolate interpolation differs from the pinned native source; neither implementation is empirical validation.')
    for row,values in zip(rows,raw[:,1:],strict=True):
        s=dict(zip(FIELDS,values.tolist(),strict=True));row['native_state']=s
        row['speeds']={mode:speed(s,mode) for mode in MODES}
        row['translation_relative_error']=abs(row['speeds']['native']/s['native_speed_m_s']-1)
        # Independent archived helper; different interpolation, not observed data.
        row['author_helper_speed']=float(helper(s['mass_kg'],s['area_m2'],s['diameter_m'],s['aspect_ratio'],
                s['viscosity_Pa_s'],s['air_density_kg_m3']))
        row['relative_errors']={k:v/row['target_speed']-1 for k,v in row['speeds'].items()}
        if s['aspect_ratio']>1:
            capacitance=.58*s['equatorial_m']/2*(1+.95*s['aspect_ratio']**.75)
            v=s['mass_kg']*s['gravity_m_s2']/(6*math.pi*s['viscosity_Pa_s']*capacitance)
            row['hubbard_douglas_hexagon']=dict(capacitance_m=capacitance,speed_m_s=v,
                reynolds_max_dimension=s['air_density_kg_m3']*v*math.hypot(s['polar_m'],s['equatorial_m'])/s['viscosity_Pa_s'],
                scope='Random orientation, viscous limit. Finite-Re drag is omitted; diagnostic only.')
    translation=max(r['translation_relative_error'] for r in rows)
    if translation>1e-11:raise ValueError('Python/native translation disagreement: '+str(translation))
    summaries=[]
    for temp in (-5.3,-14.4):
        for group in ('reference_axes','sphere','axes'):
            selected=[r for r in rows if r['temperature_C']==temp and r['start']==group]
            for mode in MODES:
                errors=[abs(r['relative_errors'][mode]) for r in selected]
                summaries.append(dict(temperature_C=temp,start=group,mode=mode,states=len(errors),
                    maximum_absolute_relative_error=max(errors),above_ten_percent=sum(e>.1 for e in errors),
                    regressions=sum(abs(r['relative_errors'][mode])>abs(r['relative_errors']['native'])+acceptance.ROUNDING for r in selected)))
    return dict(schema='terluna.atmosphere.mcsnow-joint-drag/1',baseline_sha256=acceptance.BASELINE_SHA256,
        evidence='Fixed mass and dimensions. No growth integration; alternative drag expressions are diagnostic, not accepted model changes.',
        modes=dict(native='Unchanged McSnow prolate/cylinder blend.',paper_best='Original max(q**0.25,q) Best-number denominator.',
            span_width='Reported full equatorial span replaces volume-equivalent cylinder width.',
            paper_best_span_width='Both preceding alternatives.',prolate_only='Prolate spheroid only.',
            cylinder_only='Native volume-equivalent cylinder only.',paper_cylinder_span='Original Best-number expression and actual span, cylinder only.',
            paper_prism_mean='Original Best number and prolate branch projected side area (3/pi)*width*length, averaging hexagonal azimuth.',
            paper_prism_min='Same, with minimum hexagonal projected width sqrt(3)/2 times corner span.',
            paper_prism_max='Same, with maximum hexagonal projected width equal to corner span.'),
        translation_maximum_relative_error=translation,source_identity=build_record['source_identity'],
        build_sha256=native.sha(CACHE/'drag_build/build.json'),run_sha256=native.sha(path/'run.json'),
        author_helper_identity=helper_identity,producer_sha256=native.sha(__file__),rows=rows,summaries=summaries)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('run','assemble'))
    args=parser.parse_args()
    model.save(PRODUCT,run() if args.action=='run' else assemble())
