"""Compile pinned author aggregation components independently of Python physics.

The unchanged source routines supply geometry, capacitance, growth, speed,
interception and option1 sticking. The driver composes their phi=1 outer kernel;
it does not run the author's collision sampler or a coupled snow column.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import build_welss_reference as base

HERE=Path(__file__).resolve().parent

PROGRAM='''program aggregation_reference
use reference
implicit none
type(t_sp) :: a,b
type(t_atmo) :: atmo
type(t_CE_boehm) :: info1,info2
real(wp) :: rh,v1,v2,re1,re2,c1,c2,g1,g2,qdot,L,V_new,phi_tend,ec,stick,kernel
integer :: ios
do
 read(*,*,iostat=ios) atmo%T,atmo%p,rh,a%m_i,a%mm,b%m_i,b%mm
 if(ios/=0) exit
 atmo%psati=psat_ice(atmo%T);atmo%psatw=psat_water(atmo%T)
 atmo%pv=rh*atmo%psati
 atmo%rho=(atmo%p-atmo%pv)/(Rd*atmo%T)+atmo%pv/(Rv*atmo%T)
 atmo%eta=dynamic_viscosity(atmo%T)
 call m2d_particle(a);call m2d_particle(b)
 if(a%phi/=1._wp.or.b%phi/=1._wp) error stop 'excluded shape'
 v1=boehmVel(atmo,a,info1);v2=boehmVel(atmo,b,info2)
 re1=v1*a%d*atmo%rho/atmo%eta;re2=v2*b%d*atmo%rho/atmo%eta
 c1=capacitance(a,re1);c2=capacitance(b,re2)
 call vapourDiffusion(a,v1,1._wp,atmo,g1,qdot,L,V_new,phi_tend)
 call vapourDiffusion(b,v2,1._wp,atmo,g2,qdot,L,V_new,phi_tend)
 ec=Ecbohm(atmo,info1,info2);stick=stick_effi_tempOnly(atmo%T)
 kernel=ec*stick*pi*(.5_wp*a%d+.5_wp*b%d)**2*abs(v1-v2)
 write(*,'(15(es26.17e3,1x))') a%d,c1,g1,v1,re1,b%d,c2,g2,v2,re2, &
    info1%SIG,info2%SIG,ec,stick,kernel
enddo
end program
'''


def cases():
    pairs=[(1e-12,1e-11),(1e-11,1e-9),(1e-9,1e-8),(1e-8,1e-7),
           (1e-7,1e-5),(1e-8,1e-8),(1e-8,1e-8*(1+1e-8))]
    temperatures=[273.15,269.15,269.15-1e-7,267.15,267.15-1e-7,
                  264.15,264.15-1e-7,260.65,260.65-1e-7,
                  256.15,256.15-1e-7,253.15,253.15-1e-7]
    rows=[]
    for t in temperatures:
        for m1,m2 in pairs:
            for mm1,mm2 in ((1.,2.),(2.,1.)):
                rows.append(dict(temperature_k=t,pressure_pa=80000.,relative_humidity_ice=1.05,
                                 mass_left_kg=m1,monomers_left=mm1,mass_right_kg=m2,monomers_right=mm2))
    return rows


def source():
    text=base.reference_source().split('program run_reference')[0]
    text=text.replace('module reference\n','module mo_sp_nml\ninteger :: sp_stick_effi=1\nend module\nmodule reference\n',1)
    text=text.rsplit('end module',1)[0]+(HERE/'mcsnow_sticking.inc').read_text()+'\nend module\n'
    return text+PROGRAM


def generate():
    code=source();rows=cases()
    fields=('temperature_k','pressure_pa','relative_humidity_ice','mass_left_kg','monomers_left','mass_right_kg','monomers_right')
    stdin='\n'.join(' '.join(format(row[k],'.17g') for k in fields) for row in rows)+'\n'
    with tempfile.TemporaryDirectory(prefix='welss-aggregation-components-') as directory:
        work=Path(directory);(work/'reference.f90').write_text(code)
        command=['gfortran','-O0','-fno-fast-math','-ffree-line-length-none','-std=f2008','reference.f90','-o','reference']
        build=subprocess.run(command,cwd=work,capture_output=True,text=True)
        if build.returncode:raise RuntimeError(build.stderr)
        output=subprocess.run([str(work/'reference')],cwd=work,input=stdin,check=True,capture_output=True,text=True).stdout
    keys=('diameter_left_m','capacitance_left_m','deposition_left_kg_s','speed_left_m_s','Re_left',
          'diameter_right_m','capacitance_right_m','deposition_right_kg_s','speed_right_m_s','Re_right',
          'SIG_left','SIG_right','interception','sticking','kernel_m3_s')
    lines=output.splitlines()
    if len(lines)!=len(rows):raise ValueError('incorrect reference case count')
    for row,line in zip(rows,lines,strict=True):row['expected']=dict(zip(keys,map(float,line.split()),strict=True))
    names=(*base.FILES,'build_welss_reference.py','mcsnow_sticking.inc','build_welss_aggregation_reference.py')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    return dict(schema='terluna.atmosphere.welss-aggregation-components/1',source_revision='1e428d7e901cbc5b91f2040e575f736ec5963a74',
                source_files={str((HERE/name).relative_to(HERE.parent.parent.parent)):sha(HERE/name) for name in names},
                compiled_source_sha256=hashlib.sha256(code.encode()).hexdigest(),
                compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],compiler_flags=command[1:-3],
                evidence='Independent compiled author components and composed phi=1 outer kernel; source parity, not a population or empirical comparison.',
                scope=dict(habit=0,agggeo=1,phi=1,sticking_option=1,rough_SIG=.25,random_outer_speed=0.),cases=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.write_text(json.dumps(generate(),indent=2,allow_nan=False)+'\n')
