"""Compile retained McSnow source slices and generate independent scalar fixtures.

Run from the repository root with --output PATH. Requires gfortran; compilation
and executable live in a TemporaryDirectory. The adapter is never imported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent

DRIVER = '''module reference
implicit none
integer, parameter :: wp=kind(1.d0)
real(wp), parameter :: pi=acos(-1._wp), pi_4=pi/4, pi_6=pi/6, z4_pi=4/pi
real(wp) :: grav
integer :: sp_habit=0
logical :: prolate_intp=.false.
type t_atmo
  real(wp) :: rho, eta
end type
type t_sp
  real(wp) :: m_w=0, m_r=0, m_f=0, m_i=0, phi=1, mm=2, d=0, p_area=0
  integer :: statusb=0
end type
type flags
  integer :: liquidOut=0, liquidSurf=1
end type
type(flags) :: sp_stat
TYPE_SOURCE
contains
subroutine spheroid_geometry(sp,r_a,r_c,V)
  type(t_sp), intent(in) :: sp
  real(wp), intent(out) :: r_a,r_c,V
  error stop 'prolate branch outside this driver scope'
end subroutine
FUNCTION_SOURCE
end module
program run_reference
use reference
implicit none
type(t_sp) :: one,two
type(t_CE_boehm) :: p1,p2
type(t_atmo) :: atmo
real(wp) :: d1,m1,q1,s1,d2,m2,q2,s2,v1,v2,ec,re1,re2
integer :: ios
do
  read(*,*,iostat=ios) d1,m1,q1,s1,d2,m2,q2,s2,grav,atmo%rho,atmo%eta
  if(ios /= 0) exit
  one%d=d1; one%m_i=m1; one%p_area=pi_4*d1**2*q1
  two%d=d2; two%m_i=m2; two%p_area=pi_4*d2**2*q2
  v1=boehmVel(atmo,one,p1); v2=boehmVel(atmo,two,p2)
  p1%sig=s1; p2%sig=s2
  ec=Ecbohm(atmo,p1,p2)
  re1=v1*d1*atmo%rho/atmo%eta; re2=v2*d2*atmo%rho/atmo%eta
  write(*,'(11(es26.17e3,1x))') v1,re1,p1%DELS,p1%CDI,v2,re2,p2%DELS,p2%CDI,ec,p1%alpha_e,p2%alpha_e
enddo
end program
'''


def cases():
    # These are specified air/gravity inputs, not observed terrestrial states.
    # The supplied mass law merely makes diverse numerical test particles.
    pairs = [(2e-5, 2e-5), (2e-5, 2.00002e-5), (2e-5, 2e-4),
             (2e-4, 2e-3), (2e-4, .02), (2e-4, .2), (.02, .2), (.2, .6)]
    out = []
    for g in (9.81, 1.62):
        for q in (.15, .35, 1.):
            for sigma in (0., .25):
                for index, (d1, d2) in enumerate(pairs):
                    # Use compact ice for the tiny tests and open-aggregate
                    # masses above 0.2 mm; this is fixture construction only.
                    def mass(d):
                        return min(800*3.141592653589793*d**3/6, .0185*d**1.9)
                    out.append(dict(D1=d1, m1=mass(d1), q1=q, sigma1=sigma,
                                    D2=d2, m2=mass(d2), q2=q, sigma2=sigma,
                                    g=g, rho=1.1 if index % 2 == 0 else .3, eta=1.7e-5))
    # Unequal area ratios and mixed internal widths exercise the pair mapping.
    for sigma1, sigma2 in ((0., .25), (.25, 0.)):
        out.append(dict(D1=.001,m1=1e-8,q1=.2,sigma1=sigma1,
                        D2=.003,m2=5e-8,q2=.8,sigma2=sigma2,g=1.62,rho=.5,eta=1.6e-5))
    out.extend(json.loads((HERE/'boehm_lunar_cases.json').read_text())['cases'])
    return out


def generate():
    paths = [HERE/'mcsnow_types.inc', HERE/'mcsnow_functions.inc']
    source = DRIVER.replace('TYPE_SOURCE', paths[0].read_text()).replace(
        'FUNCTION_SOURCE', paths[1].read_text())
    inputs = cases()
    input_keys = ('D1','m1','q1','sigma1','D2','m2','q2','sigma2','g','rho','eta')
    stdin = '\n'.join(' '.join(format(row[k], '.17g') for k in input_keys) for row in inputs)+'\n'
    with tempfile.TemporaryDirectory(prefix='snow-boehm-fortran-') as directory:
        work = Path(directory)
        (work/'reference.f90').write_text(source)
        command = ['gfortran','-O0','-fno-fast-math','-ffree-line-length-none',
                   '-std=f2008','reference.f90','-o','reference']
        subprocess.run(command, cwd=work, check=True, capture_output=True, text=True)
        result = subprocess.run([str(work/'reference')],input=stdin,cwd=work,
                                check=True,capture_output=True,text=True)
    expected_keys = ('v1','Re1','delta1','CDI1','v2','Re2','delta2','CDI2','Ec','alpha_e1','alpha_e2')
    outputs = result.stdout.splitlines()
    if len(outputs) != len(inputs):
        raise ValueError('unexpected Fortran output count')
    for row, line in zip(inputs, outputs):
        row['expected'] = dict(zip(expected_keys, map(float,line.split()), strict=True))
    sha = lambda raw: hashlib.sha256(raw).hexdigest()
    return dict(schema='terluna.atmosphere.snow-boehm-reference/1',
                evidence='Independent scalar Fortran source parity; no empirical or lunar validation.',
                source_revision='1e428d7e901cbc5b91f2040e575f736ec5963a74',
                compiler=subprocess.run(['gfortran','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],
                compiler_flags=command[1:-3],
                reference_source_sha256=sha(source.encode()),
                source_slices={p.name:sha(p.read_bytes()) for p in paths},
                generation_script_sha256=sha(Path(__file__).read_bytes()),
                additional_inputs_sha256={name:sha((HERE/name).read_bytes()) for name in
                                          ('boehm_lunar_cases.json','mcsnow_ground_graupel.json')},
                input_construction='Explicit synthetic cases; g=9.81 and 1.62 m/s2, rho=0.3–1.1 kg/m3. Masses specify numerical tests, not observations.',
                cases=inputs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(generate(),indent=2,allow_nan=False)+'\n')
