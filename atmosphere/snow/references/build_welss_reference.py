"""Compile pinned McSnow diagnostic-monomer microphysics independently.

    python atmosphere/snow/references/build_welss_reference.py --output PATH

No Python microphysics implementation is imported. The synthetic air states
exercise component formulas, not archived atmosphere/population identity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
FILES = ('mcsnow_welss_constants.inc', 'mcsnow_welss_functions.inc',
         'mcsnow_types.inc', 'mcsnow_functions.inc', 'welss_reference_sources.json',
         'mcsnow_bsd_3_clause.txt')

DRIVER = '''CONSTANT_SOURCE
module reference
use mo_constants
use mo_atmo_types
use geometry_constants
implicit none
integer :: sp_habit=0, agggeo_type=1, snowgeo_type=2
logical :: prolate_intp=.false.
type t_sp
  real(wp) :: m_w=0, m_r=0, m_f=0, m_i=0, phi=1, mm=1, d=0, p_area=0, v_i=0, v_r=0
  integer :: statusb=0
end type
type flags
  integer :: liquidOut=0, liquidSurf=1, liquidFillin=2
end type
type(flags) :: sp_stat
TYPE_SOURCE
contains
subroutine habit_predict(particle,mdot,dt,T,Dv,vt,ReSc,V_new,phi_tend)
  type(t_sp), intent(in) :: particle
  real(wp), intent(in) :: mdot,dt,T,Dv,vt,ReSc
  real(wp), intent(out) :: V_new,phi_tend
  error stop 'habit prediction outside diagnostic-monomer fixture'
end subroutine
BOEHM_SOURCE
FUNCTION_SOURCE
end module
program run_reference
use reference
implicit none
type(t_sp) :: particle
type(t_atmo) :: atmo
real(wp) :: relative_ice,v,re,C,Dv,Ka,L,mdot,qdot,V_new,phi_tend,ventv,venth,vcrit
integer :: ios
do
  read(*,*,iostat=ios) atmo%T,atmo%p,relative_ice,particle%m_i
  if(ios /= 0) exit
  atmo%psati=psat_ice(atmo%T); atmo%psatw=psat_water(atmo%T)
  atmo%pv=relative_ice*atmo%psati
  atmo%rho=(atmo%p-atmo%pv)/(Rd*atmo%T)+atmo%pv/(Rv*atmo%T)
  atmo%eta=dynamic_viscosity(atmo%T)
  call m2d_particle(particle)
  if(particle%mm/=1._wp .or. particle%phi/=1._wp .or. particle%statusb/=0) error stop 'scope'
  v=boehmVel(atmo,particle)
  re=v*particle%d*atmo%rho/atmo%eta
  C=capacitance(particle,re)
  Dv=diffusivity_vapour(atmo%T,atmo%p); Ka=thermalCond_dry(atmo%T)
  ventv=ventilation_coef_i(sqrt(re)*(atmo%eta/(atmo%rho*Dv))**z1_3,particle%phi)
  venth=ventilation_coef_i(sqrt(re)*(atmo%eta*c_pa/Ka)**z1_3,particle%phi)
  call vapourDiffusion(particle,v,1._wp,atmo,mdot,qdot,L,V_new,phi_tend)
  vcrit=pi_6*particle%d**3-particle%m_i*rhoii
  write(*,'(21(es26.17e3,1x))') atmo%pv,atmo%rho,atmo%eta,atmo%psati,atmo%psatw, &
    Dv,Ka,L,particle%d,particle%p_area,particle%p_area/(pi_4*particle%d**2),C,v,re,mdot, &
    ventv,venth,vcrit,C/particle%d,mth,Dth
enddo
end program
'''


def cases():
    # Values near the analytic compact-ice/power-law transition deliberately
    # probe the source's cancellation-sensitive capacitance branch.
    alpha = .0028 * 10**(2*2.1-3)
    diameter_transition = (919*math.pi/6/alpha)**(1/(2.1-3))
    mass_transition = 919*math.pi/6*diameter_transition**3
    masses = (1e-20, 1e-16, 1e-14, 1e-12, 1e-11,
              mass_transition*(1-1e-6), mass_transition*(1+1e-6),
              1e-10, 1e-9, 1e-7, 1e-5)
    states = ((273.15, 101325., 1.), (263.15, 85000., 1.05),
              (253.15, 70000., 1.10), (243.15, 58000., .95))
    return [dict(temperature_k=t, pressure_pa=p, relative_humidity_ice=rh, mass_kg=m)
            for t,p,rh in states for m in masses]


def reference_source():
    return (DRIVER.replace('CONSTANT_SOURCE', (HERE/FILES[0]).read_text())
            .replace('FUNCTION_SOURCE', (HERE/FILES[1]).read_text())
            .replace('TYPE_SOURCE', (HERE/FILES[2]).read_text())
            .replace('BOEHM_SOURCE', (HERE/FILES[3]).read_text()))


def generate():
    source = reference_source()
    inputs = cases()
    input_keys = ('temperature_k','pressure_pa','relative_humidity_ice','mass_kg')
    stdin = '\n'.join(' '.join(format(row[k],'.17g') for k in input_keys) for row in inputs)+'\n'
    with tempfile.TemporaryDirectory(prefix='welss-microphysics-fortran-') as directory:
        work = Path(directory)
        (work/'reference.f90').write_text(source)
        command = ['gfortran','-O0','-fno-fast-math','-ffree-line-length-none',
                   '-std=f2008','reference.f90','-o','reference']
        compiled = subprocess.run(command,cwd=work,capture_output=True,text=True)
        if compiled.returncode:
            raise RuntimeError(compiled.stderr)
        result = subprocess.run([str(work/'reference')],input=stdin,cwd=work,
                                check=True,capture_output=True,text=True)
    keys = ('vapour_pressure_pa','air_density_kg_m3','dynamic_viscosity_pa_s',
            'ice_saturation_pressure_pa','water_saturation_pressure_pa',
            'vapour_diffusivity_m2_s','thermal_conductivity_w_m_k','latent_heat_j_kg',
            'diameter_m','projected_area_m2','area_ratio','capacitance_m',
            'terminal_speed_m_s','Re','deposition_kg_s','vapour_ventilation_factor',
            'heat_ventilation_factor','v_crit_m3','capacitance_over_diameter',
            'source_transition_mass_kg','source_transition_diameter_m')
    lines = result.stdout.splitlines()
    if len(lines) != len(inputs):
        raise ValueError('Unexpected independent Fortran result count')
    for row,line in zip(inputs,lines):
        row['expected'] = dict(zip(keys,map(float,line.split()),strict=True))
        if not all(math.isfinite(x) for x in row['expected'].values()):
            raise ValueError('Nonfinite independent source output')
    sha = lambda raw: hashlib.sha256(raw).hexdigest()
    return dict(schema='terluna.atmosphere.welss-microphysics-reference/1',
                source_revision='1e428d7e901cbc5b91f2040e575f736ec5963a74',
                evidence='Independent scalar Fortran source parity at synthetic specified air states. '
                         'No archived-column identity, population agreement or empirical validation.',
                scope=dict(sp_habit=0,agggeo_type=1,snowgeo_type=2,mm=1,phi=1,
                           pure_ice=True,rough_surface=True,ventilation=True,
                           atmosphere_interpolation_tested=False),
                source_files={name:sha((HERE/name).read_bytes()) for name in FILES},
                generation_script_sha256=sha(Path(__file__).read_bytes()),
                compiled_source_sha256=sha(source.encode()),
                compiler=subprocess.run(['gfortran','--version'],check=True,capture_output=True,text=True).stdout.splitlines()[0],
                compiler_flags=command[1:-3],
                cases=inputs,
                reading_rule='Source functions and literals are retained unchanged, including unsuffixed '
                             'default-real constants and the v_crit>0 capacitance decision. Air density uses '
                             'the source ideal moist-gas expression at explicit T,p,pv; no grid interpolation '
                             'or cloud-water repartition is tested. The habit-prediction stub is a fail-fast '
                             'guard for an excluded branch and is never called. No supply cap is applied.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(generate(),indent=2,allow_nan=False)+'\n')
