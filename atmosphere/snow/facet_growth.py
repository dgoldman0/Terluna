"""Stagnant compact-crystal growth with published facet attachment coefficients.

An explicit development control: constant dense ice, no branching, hollowing,
ventilation or drag. It does not replace the lunar population growth law.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from atmosphere.snow import facet_attachment as attachment
from atmosphere.snow import growth_humidity_validation as humidity
from atmosphere.snow import flakes
from shared.provenance import constants_used

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BUILD=ROOT/'research/runs/snow_scale/facet_attachment/double_stable_final'
RESULT=HERE/'results/facet_growth.json'
ICE_DENSITY=919. # Native McSnow comparison convention, not a new global constant.


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def effective_transport(state,a,c):
    """ZH14/ICTG axis-weighted resistances, f_v=f_T=1, native alpha_T=.96."""
    if not all(math.isfinite(x) and x>0 for x in (a,c)):raise ValueError('Positive finite semiaxes required')
    required=('capacitance_m','boundary_capacitance_m','diffusivity_m2_s','conductivity_W_m_K',
              'vapor_molecular_speed_m_s','alpha_a','alpha_c','thermal_jump_length_m')
    if any(not math.isfinite(state[k]) or state[k]<=0 for k in required):
        raise ValueError('Finite positive source transport and attachment values required')
    if any(state[k]>1 for k in ('alpha_a','alpha_c')):
        raise ValueError('Attachment efficiencies cannot exceed unity')
    if state['boundary_capacitance_m']<state['capacitance_m']:
        raise ValueError('Enlarging the particle cannot reduce capacitance')
    cap=state['capacitance_m'];ratio=cap/state['boundary_capacitance_m']
    dv=state['diffusivity_m2_s'];ka=state['conductivity_W_m_K'];speed=state['vapor_molecular_speed_m_s']
    dv_a=dv/(4*dv/(state['alpha_a']*speed)*cap/(a*c)+ratio)
    dv_c=dv/(4*dv/(state['alpha_c']*speed)*cap/(a*a)+ratio)
    lt=state['thermal_jump_length_m']
    ka_a=ka/(lt*cap/(a*c)+ratio);ka_c=ka/(lt*cap/(a*a)+ratio)
    return (2*dv_a+dv_c)/3,(2*ka_a+ka_c)/3


def mass_and_axis_rates(state,a,c,temperature,ice_supersaturation,density=ICE_DENSITY):
    """Facet dc/da=alpha_c/alpha_a, conserving dense spheroidal volume growth."""
    if not math.isfinite(temperature) or temperature<=0:raise ValueError('Positive finite temperature required')
    for key in ('vapor_gas_constant_J_kg_K','latent_sublimation_J_kg','ice_saturation_pressure_Pa'):
        if not math.isfinite(state[key]) or state[key]<=0:raise ValueError('Invalid thermodynamic source value: '+key)
    if not math.isfinite(density) or density<=0:raise ValueError('Positive deposition density required')
    if not math.isfinite(ice_supersaturation) or ice_supersaturation<0:raise ValueError('Growth-only control')
    dv,ka=effective_transport(state,a,c);rv=state['vapor_gas_constant_J_kg_K'];latent=state['latent_sublimation_J_kg']
    rhoeq=state['ice_saturation_pressure_Pa']/(rv*temperature)
    resistance=1/dv+rhoeq*latent/(ka*temperature)*(latent/(rv*temperature)-1)
    if not math.isfinite(resistance) or resistance<=0:raise ValueError('Invalid total transport resistance')
    deff=1/resistance
    dm=4*math.pi*state['capacitance_m']*rhoeq*ice_supersaturation*deff
    q=state['alpha_c']/state['alpha_a']
    da=dm/(density*(4*math.pi/3)*(2*a*c+a*a*q))
    return dict(mass_rate_kg_s=dm,equatorial_rate_m_s=da,polar_rate_m_s=q*da,
        effective_diffusivity_m2_s=dv,effective_conductivity_W_m_K=ka,
        coupled_diffusivity_m2_s=deff,facet_rate_ratio=q)


def integrate(library,datum,phi,m,rtol=1e-7,pressure=101000.):
    """Fixed observed forcing; unknown initial shape and source mechanisms are scenarios."""
    if phi not in (.1,.3,1.):raise ValueError('Undeclared seed shape')
    if m not in (1.,10.):raise ValueError('Undeclared attachment mechanism')
    if not math.isfinite(rtol) or rtol<=0:raise ValueError('Positive solver tolerance required')
    t=datum['nominal_temperature_C']+273.15
    e=flakes.ice_saturation_pressure(t)+datum['vapor_density_excess_g_m3']*1e-3*flakes.R_VAPOUR*t
    a0=5e-6;c0=phi*a0;m0=ICE_DENSITY*4*math.pi/3*a0*a0*c0
    reference=library.evaluate(.01,t,pressure,a0,c0,m)
    si=e/reference['ice_saturation_pressure_Pa']-1
    if si<=0:raise ValueError('Forcing is not supersaturated under source convention')
    # Log axes avoid inadmissible trial geometry. Integrated deposition is scaled
    # by initial mass; it is independent of the geometric mass diagnostic.
    def rhs(seconds,y):
        a=a0*math.exp(y[0]);c=c0*math.exp(y[1])
        state=library.evaluate(si,t,pressure,a,c,m)
        rates=mass_and_axis_rates(state,a,c,t,si)
        return [rates['equatorial_rate_m_s']/a,rates['polar_rate_m_s']/c,rates['mass_rate_kg_s']/m0]
    solution=solve_ivp(rhs,(0.,600.),[0.,0.,0.],rtol=rtol,atol=rtol*1e-2,t_eval=np.arange(0.,601.,60.),max_step=5.)
    rows=[]
    for seconds,y in zip(solution.t,solution.y.T):
        a=a0*math.exp(y[0]);c=c0*math.exp(y[1]);mass=ICE_DENSITY*4*math.pi/3*a*a*c
        state=library.evaluate(si,t,pressure,a,c,m)
        rates=mass_and_axis_rates(state,a,c,t,si)
        rows.append(dict(elapsed_s=float(seconds),mass_kg=mass,equatorial_dimension_m=2*a,polar_dimension_m=2*c,
            diameter_m=2*max(a,c),aspect_ratio=c/a,deposited_mass_kg=float(y[2])*m0,
            relative_mass_volume_residual=(mass-m0-float(y[2])*m0)/mass,
            attachment={k:state[k] for k in ('alpha_a','alpha_c','alpha_t')},rates=rates))
    complete=bool(solution.success and rows and rows[-1]['elapsed_s']==600. and all(math.isfinite(r['mass_kg']) for r in rows))
    return dict(datum_id=datum['id'],nominal_temperature_C=datum['nominal_temperature_C'],
        vapor_density_excess_g_m3=datum['vapor_density_excess_g_m3'],pressure_Pa=pressure,
        imposed_vapor_pressure_Pa=e,native_ice_supersaturation=si,initial_full_equatorial_m=2*a0,
        initial_thickness_over_span=phi,initial_mass_kg=m0,deposition_density_kg_m3=ICE_DENSITY,
        mechanism_m=m,relative_solver_tolerance=rtol,complete=complete,solver_message=solution.message,
        function_evaluations=solution.nfev,trajectory=rows,observed_diameter_mm=datum['observed_diameter_mm'],
        endpoint_diameter_mm=rows[-1]['diameter_m']*1000 if complete else None,
        endpoint_mass_kg=rows[-1]['mass_kg'] if complete else None,
        endpoint_diameter_ratio=rows[-1]['diameter_m']*1000/datum['observed_diameter_mm'] if complete else None,
        maximum_relative_mass_volume_residual=max(abs(r['relative_mass_volume_residual']) for r in rows) if rows else None)


def calculate(build_dir=BUILD):
    library=attachment.AttachmentLibrary(build_dir);identity=library.identity
    data=humidity.observations()
    # The 0.05g/m3 selection is the dense end of the cited deposition-density law;
    # this control never extrapolates it to the large-branched regime.
    chosen=[row for row in data['cases'] if row['vapor_density_excess_g_m3']<=.05]
    cases=[];refinements=[]
    for row in chosen:
        for m in (1.,10.):
            for phi in (.1,.3,1.):
                coarse=integrate(library,row,phi,m,1e-6);fine=integrate(library,row,phi,m,1e-8)
                cases.append(fine)
                refinements.append(dict(datum_id=row['id'],mechanism_m=m,initial_thickness_over_span=phi,
                    complete_pair=coarse['complete'] and fine['complete'],
                    diameter_relative_change=abs(coarse['endpoint_diameter_mm']/fine['endpoint_diameter_mm']-1) if coarse['complete'] and fine['complete'] else None,
                    mass_relative_change=abs(coarse['endpoint_mass_kg']/fine['endpoint_mass_kg']-1) if coarse['complete'] and fine['complete'] else None))
    if identity!=attachment.verify_build(build_dir):raise ValueError('Attachment build changed during integration')
    paths=[Path(__file__),Path(attachment.__file__),attachment.DRIVER,Path(flakes.__file__),Path(humidity.__file__),humidity.DATA,HERE/'facet_attachment_sources.json']
    return dict(schema='terluna.atmosphere.compact-facet-growth/1',
        producer={str(p.relative_to(ROOT)):sha(p) for p in paths},constants=constants_used(paths),
        attachment_build_record_sha256=sha(Path(build_dir)/'build.json'),attachment_source_sha256=attachment.SOURCE_SHA256,
        integration=dict(solver='scipy.solve_ivp RK45',scipy_version=scipy.__version__,numpy_version=np.__version__,
            maximum_step_s=5.,relative_tolerances=[1e-6,1e-8],absolute_tolerance_factor=.01),
        evidence='Stagnant dense-ice development control with stable spheroid capacitance; additional Earth comparison with explicit forcing/seed assumptions, not a complete growth repair.',
        equations={'source':'2021 facet coefficient routine; ZH2014/ICTG2022 axis-weighted transport and2021 dimensionally consistent mass law.',
            'mass':'dm/dt=4pi C rho_eq s_i/[1/Dv_prime +rho_eq L/(K_prime T)(L/(Rv T)-1)]',
            'shape':'dc/da=alpha_c/alpha_a; da/dt=(dm/dt)/[rho_dep (4pi/3)(2ac+a^2 alpha_c/alpha_a)]',
            'density':'rho_dep=919 kg/m3 throughout; no branching/hollowing.',
            'flow':'vapor and thermal ventilation both1; nativealpha helper has no ventilation input.',
            'printed_equation_corrections':['The facet/log-volume algebra requires Gamma/phi+2 in the denominator, not the minus2 printed in2021Eq3.',
                'The mass law requires equilibrium vapor density rho_eq (2021Eq2); ICTG2022Eq2 printed image omits it.']},
        cases=cases,refinements=refinements,
        limits=['Selected vapor-density excess <=0.05g/m3 only; no dense-ice claim for branched high-supersaturation crystals.',
            'Observed crystals were freely suspended; no-ventilation is a declared approximation, not a reconstruction of their airflow.',
            'Unknown seed shapes and mechanisms remain scenarios; no coefficients or initial shapes fitted to measured sizes.',
            'Mass and thickness are model predictions but unmeasured in this2018 comparison; size agreement cannot validate them.',
            'Only source air properties and fixed Earth pressure; no lunar atmosphere or aggregate population forecast.',
            'The source alpha approximation uses bulk thermal resistance; small-particle thermal-jump consistency is approximate, not exact implicit DiSKICE parity.',
            'Single-grain compact kinetics alone cannot close the branching/hollowing, finite-reservoir or fracture problems.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--build',type=Path,default=BUILD);parser.add_argument('--output',type=Path,default=RESULT)
    args=parser.parse_args();product=calculate(args.build);args.output.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(len(product['cases']),'refined cases;',sum(c['complete'] for c in product['cases']),'complete')
if __name__=='__main__':main()
