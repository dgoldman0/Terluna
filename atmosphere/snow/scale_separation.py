"""Test offset contacts between oblate ice particles in the pinned SCALE model.

Experimental extension of Shima et al. (2020), section 9.3.8. A prescribed
separation ratio changes the equatorial extent of two oblate particles; the
native density interpolation then uses that new extent. This is a sensitivity,
not a calibrated contact model or a replacement for the preserved baseline.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import difflib
import hashlib
import json
import math
import re
import shutil
import subprocess
import numpy as np
from atmosphere.snow import scale_size_diagnosis as size
from atmosphere.snow import scale_build

lunar = size.lunar
CACHE = lunar.CACHE
SOURCE = 'contrib/SDM/sdm_coalescence_cold.f90'
SOURCE_SHA = '22573a3b9241c3fd446ca2853c67674bc05671704b7da08c6079f1372120b487'
ROUTINE = 'sdm_outcome_aggregation_shima5_multi'
OUT = lunar.ROOT/'atmosphere/snow/results/scale_separation.json'
CASES = [
    dict(name='separation_zero_moon', separation=0.),
    dict(name='separation_half_moon', separation=.5),
    dict(name='separation_half_moon_dt15', separation=.5, dt=15., packets=1),
    dict(name='separation_half_earth_dt15', separation=.5, gravity=size.STANDARD_GRAVITY,
         depth_factor=size.DEPTH_RATIO, dt=15., packets=3),
    dict(name='separation_half_earth_dt5', separation=.5, gravity=size.STANDARD_GRAVITY,
         depth_factor=size.DEPTH_RATIO, dt=5., packets=1),
]
COUNTERS = '''  integer(kind=8), save :: separation_calls=0, separation_oblate=0, separation_multiple=0
'''


def wrapper(separation):
    if not math.isfinite(separation) or not 0 <= separation <= 1:
        raise ValueError('Separation ratio must lie in [0,1]')
    return f'''
  subroutine {ROUTINE}(sd_ncol,a1,a2,c1,c2,rho1,rho2)
    use scale_const, only: rhoi_mks => CONST_DICE
    use m_sdm_common, only: F_THRD,ONE_PI
    integer,intent(in)::sd_ncol
    real(RP),intent(in)::a1,c1,rho1
    real(RP),intent(inout)::a2,c2,rho2
    real(RP),parameter::separation={separation:.17e}_RP, rho_min=10.0_RP
    real(RP)::mass,vmin,rho_average,anew,cmin,cmax,cnew
    integer::event
    separation_calls=separation_calls+1
    if(sd_ncol>1)separation_multiple=separation_multiple+1
    if(separation==0.0_RP .or. a1<=c1 .or. a2<=c2)then
      call {ROUTINE}_native(sd_ncol,a1,a2,c1,c2,rho1,rho2)
      return
    endif
    ! Sequential binary encounters for a multiple-event SDM pair. Column runs
    ! also record how often this branch is needed; zero separation retains the
    ! original native batched operation exactly.
    do event=1,sd_ncol
      if(a1>c1 .and. a2>c2)then
        separation_oblate=separation_oblate+1
        mass=F_THRD*ONE_PI*(a1*a1*c1*rho1+a2*a2*c2*rho2)
        vmin=F_THRD*ONE_PI*(a1*a1*c1+a2*a2*c2)
        rho_average=mass/vmin
        ! D'=max(D1,D2,(1+s)(D1+D2)/2), with D=2a for oblate inputs.
        anew=max(a1,a2,(1.0_RP+separation)*(a1+a2)/2.0_RP)
        cmin=vmin/(F_THRD*ONE_PI*anew*anew)
        cmax=c1+c2
        ! Retain the native density interpolation, recomputing its two volume
        ! endpoints for the changed equatorial extent.
        cnew=(rhoi_mks-rho_min)/((rhoi_mks-rho_average)/cmin+(rho_average-rho_min)/cmax)
        a2=anew;c2=cnew;rho2=mass/(F_THRD*ONE_PI*a2*a2*c2)
      else
        call {ROUTINE}_native(1,a1,a2,c1,c2,rho1,rho2)
      endif
    enddo
  end subroutine {ROUTINE}
'''


def native_routine(text):
    start=text.index('  subroutine '+ROUTINE+'(')
    end=text.index('  end subroutine '+ROUTINE,start)+len('  end subroutine '+ROUTINE)
    return text[start:end]


def patch(text,separation):
    if hashlib.sha256(text.encode()).hexdigest()!=SOURCE_SHA:
        raise ValueError('Unidentified aggregation source')
    added=wrapper(separation)
    original=native_routine(text)
    result=text.replace(original,original.replace(ROUTINE,ROUTINE+'_native')+'\n'+added)
    if result.count('\ncontains\n')!=1:raise ValueError('Unexpected module structure')
    result=result.replace('\ncontains\n','\n'+COUNTERS+'\ncontains\n',1)
    anchor='    close(observer_unit)\n'
    if result.count(anchor)!=1:raise ValueError('Missing collision observer')
    result=result.replace(anchor,anchor+'''    if(observer_calls==1)then
      open(newunit=observer_unit,file='separation_observer.txt',status='new',action='write')
      write(observer_unit,'(A)') '# call cumulative_outcome_calls cumulative_oblate_events cumulative_multiple_pairs'
    else
      open(newunit=observer_unit,file='separation_observer.txt',status='old',position='append',action='write')
    endif
    write(observer_unit,'(4(I20,1X))') observer_calls,separation_calls,separation_oblate,separation_multiple
    close(observer_unit)
''')
    return result


def build_name(separation):
    if separation not in (0.,.5):raise ValueError('Declared builds are 0 and 0.5 only')
    return 'separation_zero' if separation==0 else 'separation_half'


def build(separation):
    target=CACHE/build_name(separation)
    if target.exists():
        record=json.loads((target/'manifest.json').read_text())
        if record['separation']!=separation or record['producer_sha256']!=lunar.base.digest(Path(__file__)):
            raise ValueError('Existing build has a different producer or separation')
        return record
    target.mkdir()
    work=target/'work'
    shutil.copytree(CACHE/'lunar/work',work,symlinks=True)
    source=work/SOURCE;original=source.read_text();source.write_text(patch(original,separation))
    diff=''.join(difflib.unified_diff(original.splitlines(True),source.read_text().splitlines(True),fromfile='a/'+SOURCE,tofile='b/'+SOURCE))
    (target/'aggregation.patch').write_text(diff)
    deps=CACHE/'deps/usr'
    result=scale_build.build(work,'-I'+str(deps/'include'),str(deps/'lib/x86_64-linux-gnu/libnetcdff.so.7')+' /usr/lib/x86_64-linux-gnu/libnetcdf.so.19 -Wl,-rpath,'+str(deps/'lib/x86_64-linux-gnu')+' -lm -lz',workers=1)
    record=dict(separation=separation,build=result,original_manifest=json.loads((CACHE/'lunar/manifest.json').read_text()),
        hook_sha256=lunar.base.digest(lunar.HOOK),producer_sha256=lunar.base.digest(Path(__file__)),
        source=SOURCE,original_sha256=SOURCE_SHA,patched_sha256=lunar.base.digest(source),patch_sha256=lunar.base.digest(target/'aggregation.patch'))
    lunar.boxes.save(target/'manifest.json',record)
    return record


def configure(name,separation,cache=CACHE,**kw):
    record=json.loads((Path(cache)/build_name(separation)/'manifest.json').read_text())
    if record['separation']!=separation:raise ValueError('Wrong separation binary')
    d=size.configure(name,cache=cache,**kw)
    s=json.loads((d/'settings.json').read_text())
    s.update(build=record,build_directory=build_name(separation),separation=separation,
        geometry_evidence='Experimental planar-contact extension; fixed separation ratio and native density interpolation. No contact-energy or observed snowfall validation.',
        geometry_producer_sha256=lunar.base.digest(Path(__file__)))
    lunar.boxes.save(d/'settings.json',s)
    lunar.boxes.save(d/'inputs.json',{p.name:lunar.base.digest(p) for p in d.iterdir() if p.is_file() and p.name!='inputs.json'})
    return d


def run(spec):
    d=CACHE/'runs'/spec['name']
    if not d.exists():d=configure(**spec)
    saved=json.loads((d/'settings.json').read_text())
    if saved['separation']!=spec['separation']:raise ValueError('Cached separation mismatch')
    if not (d/'analysis.json').exists():
        if not (d/'execution.json').exists():lunar.execute(d)
        lunar.analyze(d)
    return d.name


def pair_outcomes(rows,separation,directory):
    """Compile and call the actual patched/native outcome routines in isolation."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=CACHE/'lunar/work'/SOURCE
    text=source.read_text()
    if lunar.base.digest(source)!=SOURCE_SHA:raise ValueError('Unidentified pair reference')
    constants=CACHE/'lunar/work/scalelib/src/common/scale_const.F90'
    ice=re.search(r'CONST_DICE\s*=\s*([0-9.]+)_RP',constants.read_text()).group(1)
    code=f'''module scale_precision
 integer,parameter::RP=kind(1.0d0)
end module
module scale_const
 use scale_precision
 real(RP),parameter::CONST_DICE={ice}_RP
end module
module m_sdm_common
 use scale_precision
 real(RP),parameter::F_THRD=4.0_RP/3.0_RP,ONE_PI=acos(-1.0_RP)
end module
module pairs
 use scale_precision
 implicit none
{COUNTERS}
contains
{native_routine(text).replace(ROUTINE,ROUTINE+'_native')}
{wrapper(separation)}
end module
program compare
 use pairs
 implicit none
 integer::ios,n
 real(RP)::a1,a2,c1,c2,r1,r2,an,cn,rn
 do
  read(*,*,iostat=ios)n,a1,a2,c1,c2,r1,r2
  if(ios/=0)exit
  an=a2;cn=c2;rn=r2
  call {ROUTINE}_native(n,a1,an,c1,cn,r1,rn)
  call {ROUTINE}(n,a1,a2,c1,c2,r1,r2)
  write(*,'(6(ES25.16,1X))')an,cn,rn,a2,c2,r2
 enddo
end program
'''
    (directory/'pairs.f90').write_text(code)
    subprocess.run(['gfortran','-O2','-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow','pairs.f90','-o','pairs'],cwd=directory,check=True,capture_output=True,text=True)
    stdin=''.join(' '.join(format(float(v),'.17g') if i else str(int(v)) for i,v in enumerate(row))+'\n' for row in rows)
    result=subprocess.run([str((directory/'pairs').resolve())],input=stdin,check=True,capture_output=True,text=True)
    return np.array([[float(v) for v in line.split()] for line in result.stdout.splitlines()])


def assemble():
    cases={};missing=[]
    names=[s['name'] for s in CASES]+['size_moon_baseline','size_earth_shallow_dt15_p3','size_earth_shallow_dt5_p1']
    for name in names:
        d=CACHE/'runs'/name
        try:
            c=json.loads((d/'analysis.json').read_text())
            cases[name]=dict(group=size.comparison.native_group(c),population=size.population(d,c),
                arrivals={k:v for k,v in c['arrivals'].items() if k!='cdf'},settings=c['settings'],
                time_blocks=c['gathering_time_blocks'],final=c['final'],observer=c['observer'],
                water_relative_error=c['max_relative_water_residual'],analysis_sha256=lunar.base.digest(d/'analysis.json'),
                input_sha256=json.loads((d/'inputs.json').read_text()),execution=c['execution'])
            if (d/'separation_observer.txt').exists():
                last=np.atleast_2d(np.loadtxt(d/'separation_observer.txt',dtype=np.int64))[-1]
                cases[name]['geometry_counters']=dict(zip(['calls','outcome_calls','oblate_events','multiple_pairs'],map(int,last)))
        except (FileNotFoundError,ValueError) as e:missing.append(dict(name=name,error=str(e)))
    pairs=[('size_moon_baseline','separation_zero_moon','zero-separation full-column parity',True),
           ('size_moon_baseline','separation_half_moon','lunar geometry sensitivity',False),
           ('separation_half_moon','separation_half_moon_dt15','lunar timestep refinement',True),
           ('size_earth_shallow_dt15_p3','separation_half_earth_dt15','Earth geometry sensitivity at 15 s',False),
           ('size_earth_shallow_dt5_p1','separation_half_earth_dt5','Earth geometry sensitivity at 5 s',False),
           ('separation_half_earth_dt15','separation_half_earth_dt5','Earth timestep refinement',True)]
    comparisons=[]
    for left,right,label,screen in pairs:
        if left not in cases or right not in cases:continue
        c=size.comparison.comparison(cases[left]['group'],cases[right]['group'],label,screen)
        c.update(left=left,right=right,query_grid_mass_cdf_difference=c['maximum_mass_cdf_difference'])
        exact=size.exact_column_cdf(left,right)
        c.update(maximum_mass_cdf_difference=exact['distance'],exact_mass_cdf_difference=exact['distance'],
                 exact_mass_cdf_difference_at_m=exact['diameter_m'])
        if screen:
            c['passed']=c['passed'] and exact['distance']<=.05
            c['screen']='At most 5% flux/d50/d90/d99 change and 0.05 exact equal-column mass-CDF distance.'
        comparisons.append(c)
    result=dict(schema='terluna.atmosphere.scale-separation/1',
        evidence='Executed sensitivity to planar contact separation; neither a validated aggregate geometry nor an observed Earth/lunar snowfall prediction.',
        producer_sha256=lunar.base.digest(Path(__file__)),cases=cases,comparisons=comparisons,missing_or_failed=missing)
    lunar.boxes.save(OUT,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['build','run','assemble'])
    p.add_argument('--separation',type=float,default=.5)
    p.add_argument('--case',action='append');p.add_argument('--workers',type=int,choices=[1,2],default=1)
    args=p.parse_args()
    if args.action=='build':print(json.dumps(build(args.separation),indent=2))
    else:
        if args.action=='run':
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                jobs=[pool.submit(run,s) for s in CASES if args.case is None or s['name'] in args.case]
                for job in as_completed(jobs):print('complete',job.result(),flush=True)
        r=assemble();print('complete',len(r['cases']),'missing',len(r['missing_or_failed']))
