"""Lunar atmosphere feasibility calculations, version 1.0.

Prescribed-temperature hydrostatic/Jeans sensitivity model, not a predictive
photochemical, radiative-convective, hydrodynamic or GCM simulation.
Run: python model.py --out results
Requires numpy, scipy, matplotlib. All dimensional quantities are SI.
"""
from pathlib import Path
import argparse, json, csv, math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 1_737_400.0
GM = 4.902800118e12
G = 6.67430e-11
KB = 1.380649e-23
NA = 6.02214076e23
RG = KB * NA
YR = 365.25 * 86400
P0 = 1.2 * 101325
TS = 288.0
XO = 0.175
MN = 0.0280134
MO = 0.031998
MU = (1-XO)*MN + XO*MO
CP = 1005.0
AREA = 4*np.pi*R**2
GRAV = GM/R**2
SIGMA = 5.5e-19  # illustrative neutral collision cross section; Tucker CO proxy
AU = 149597870700.0
GMS = 1.32712440018e20
GME = 3.986004354e14
EM = 384400000.0

def jeans(r, t, n, molar_mass):
    m = molar_mass/NA
    lam = GM*m/(KB*t*r)
    flux = n*np.sqrt(KB*t/(2*np.pi*m))*(1+lam)*np.exp(-lam)
    return float(4*np.pi*r*r*m*flux), float(lam)

def atmosphere(tc=180., te=250., ph=0.1, sigma=SIGMA, points=30001,
               isothermal=False):
    # x=ln(p_surface/p), integrate d(1/r)/dx=-R_specific*T/GM.
    x = np.linspace(0, 36, points)
    p = P0*np.exp(-x)
    if isothermal:
        t = np.full_like(x,te)
    else:
        t = np.maximum(tc, TS*np.exp(-(RG/MU)/CP*x))
        xh = np.log(P0/ph)
        hot = x > xh
        t[hot] = tc+(te-tc)*(1-np.exp(-(x[hot]-xh)/3.0))
    invr = 1/R-(RG/MU)/GM*cumulative_trapezoid(t,x,initial=0)
    valid = invr > 1/(30*R)
    x,p,t,invr = [a[valid] for a in [x,p,t,invr]]
    r=1/invr
    kn = (MU/NA)*GM/(sigma*p*r*r)
    crossings = np.where(kn>=1)[0]
    if len(crossings)==0:
        return {'status':'no_exobase_below_30R','tc':tc,'te':te,'ph':ph},None
    i=int(crossings[0])
    z=np.log(kn[i-1:i+1]); w=(0-z[0])/(z[1]-z[0])
    rx=float(r[i-1]+w*(r[i]-r[i-1])); tx=float(t[i-1]+w*(t[i]-t[i-1]))
    px=float(np.exp(np.log(p[i-1])+w*(np.log(p[i])-np.log(p[i-1]))))
    x,p,t,r,kn=[a[:i+1] for a in [x,p,t,r,kn]]
    r[-1]=rx;p[-1]=px;t[-1]=tx;kn[-1]=1
    mass=float(4*np.pi/GM*np.trapezoid(r**4*p,x))
    column=float(np.trapezoid(r**2*p,x)/GM)
    n=px/(KB*tx)
    qn,ln=jeans(rx,tx,n*(1-XO),MN)
    qo,lo=jeans(rx,tx,n*XO,MO)
    q=qn+qo
    # Post-processing only: fixed total exobase n, no back-reaction on profile.
    atom_q={}
    for f in [1e-6,1e-4,.01,.1]:
        qa,_=jeans(rx,tx,n*f,0.0140067)
        atom_q[str(f)]=qa
    # Overlying molecular column; truncated at nominal exobase, adequate
    # only where optical-depth crossing is well below that boundary.
    nd=p/(KB*t)
    integ=cumulative_trapezoid(nd,r,initial=0)
    over=integ[-1]-integ
    absorption={}
    for sig in [1e-22,1e-21,1e-20]:
        # Reverse arrays to make optical depth increasing for interpolation.
        absorption[str(sig)]=float(np.interp(1,over[::-1]*sig,r[::-1])/R)
    summary=dict(status='hydrostatic_proxy',tc=tc,te=te,ph=ph,sigma=sigma,
                 mass_kg=mass,column_kg_m2=column,exo_R=rx/R,exo_T=tx,
                 exo_p_Pa=px,lambda_N2=ln,loss_N2_kg_s=qn,loss_O2_kg_s=qo,
                 loss_kg_s=q,lifetime_yr=mass/q/YR,
                 atomic_N_loss_fixed_profile=atom_q,absorption_R=absorption)
    return summary,dict(x=x,p=p,t=t,r=r,kn=kn)

def hohmann(a_outer):
    a1=AU; a2=a_outer*AU; a=(a1+a2)/2
    vo=np.sqrt(GMS/a2); vi=np.sqrt(GMS/a1)
    vto=np.sqrt(GMS*(2/a2-1/a)); vti=np.sqrt(GMS*(2/a1-1/a))
    return dict(source_AU=a_outer,departure_dv_km_s=float((vo-vto)/1000),
        arrival_vinf_km_s=float((vti-vi)/1000),
        total_endpoint_dv_km_s=float((vo-vto+vti-vi)/1000),
        flight_yr=float(np.pi*np.sqrt(a**3/GMS)/YR))

def diffusive_sensitivity(fatom=0.,tc=180.,te=250.,ph=.1):
    """Separate-species diffusive equilibrium above a prescribed lower boundary.

    Atomic fraction is imposed at ph, not predicted chemically. Temperature
    is interpolated from the reference prescribed profile; hence this is a
    stress test, not a self-consistent photochemical solution.
    """
    ref,pf=atmosphere(tc,te,ph)
    if pf is None:return dict(status='reference_profile_failed')
    rh=float(np.interp(np.log(ph),np.log(pf['p'][::-1]),pf['r'][::-1]))
    rs=np.geomspace(rh,30*R,60001)
    ts=np.interp(rs,pf['r'],pf['t'],right=te)
    integ=cumulative_trapezoid(GM/(KB*ts*rs**2),rs,initial=0)
    mols=np.array([MN,MO,MN/2,MO/2]);ms=mols/NA
    fractions=np.array([(1-fatom)*(1-XO),(1-fatom)*XO,fatom*(1-XO),fatom*XO])
    ps=ph*fractions[:,None]*np.exp(-ms[:,None]*integ)
    pt=ps.sum(axis=0);mbar=(ps*ms[:,None]).sum(axis=0)/pt
    kn=mbar*GM/(SIGMA*pt*rs**2)
    hit=np.where(kn>=1)[0]
    if len(hit)==0:
        return dict(fatom=fatom,tc=tc,te=te,ph=ph,status='no_exobase_below_30R')
    i=hit[0];r=rs[i];t=ts[i];p=ps[:,i];frac=p/p.sum()
    rates=[];lams=[]
    for mass,pp in zip(mols,p):
        q,lam=jeans(r,t,pp/(KB*t),mass);rates.append(q);lams.append(lam)
    return dict(fatom=fatom,tc=tc,te=te,ph=ph,status='diffusive_proxy',
      exo_R=float(r/R),exo_T=float(t),atomic_number_fraction_exo=float(frac[2:].sum()),
      lambda_N=float(lams[2]),loss_N2_kg_s=rates[0],loss_O2_kg_s=rates[1],
      loss_N_kg_s=rates[2],loss_O_kg_s=rates[3],total_kg_s=sum(rates),
      inventory_timescale_yr=ref['mass_kg']/sum(rates)/YR)

def growth(g,m0=1e9,start=150,end=500,cap=1e17):
    # m0 and cap are net atmospheric delivery capacity in kg/year.
    duration=end-start
    if g==0:return m0*duration,m0
    hit=np.log(cap/m0)/g
    if hit>=duration:return m0*np.expm1(g*duration)/g,m0*np.exp(g*duration)
    return (cap-m0)/g+cap*(duration-hit),cap

def savecsv(path,rows):
    keys=[k for k in rows[0] if not isinstance(rows[0][k],dict)]
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main(out):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':160,'savefig.bbox':'tight'})
    baseline,profile=atmosphere()
    hill=EM*(GM/GME/3)**(1/3)
    xi=hill/(baseline['exo_R']*R)
    tidal_k=1-1.5/xi+.5/xi**3
    lam=baseline['lambda_N2']
    tidal=dict(Hill_radius_m=hill,barrier_K=tidal_k,
      fixed_profile_N2_Jeans_multiplier=(1+tidal_k*lam)/(1+lam)*np.exp(lam*(1-tidal_k)),
      limitation='Directional restricted-three-body sensitivity, not global loss correction')
    cases=[]
    for tc in [150,180,210]:
      for te in [230,250,260,300,350,400,600]:
       for ph in [.001,.1,10]:
        s,_=atmosphere(tc,te,ph);cases.append(s)
    savecsv(out/'atmosphere_sensitivity.csv',cases)
    diffusive=[diffusive_sensitivity(f,180,t,.1)
               for t in [250,260,300,350] for f in [0,1e-8,1e-6,1e-4,.01]]
    savecsv(out/'diffusive_sensitivity.csv',diffusive)
    mass=baseline['mass_kg'];wn=(1-XO)*MN/MU;wo=XO*MO/MU
    inventory=dict(thin_mass_kg=AREA*P0/GRAV,mass_kg=mass,N2_kg=mass*wn,O2_kg=mass*wo,
       density_kg_m3=P0*MU/(RG*TS),scale_height_m=RG*TS/(MU*GRAV),
       surface_g=GRAV,mean_molar_mass=MU,column_kg_m2=baseline['column_kg_m2'],
       airship_helium_net_kg_m3=P0*(MU-.0040026)/(RG*TS),
       flight_speed_ratio=np.sqrt((GRAV/9.80665)/(P0*MU/(101325*.028965))),
       water_10m_kg=AREA*1e4,water_100m_kg=AREA*1e5)
    logistics=[]
    for years in [200,300,500]:
      qn=mass*wn/(years*YR);qo=mass*wo/(years*YR)
      for dv,ve in [(10000,10000),(10000,30000),(10000,50000),(20000,30000),(20000,50000)]:
        prop=np.expm1(dv/ve)
        e=.5*ve**2*prop/.7
        logistics.append(dict(years=years,total_rate_kg_s=qn+qo,N2_rate_kg_s=qn,
            O2_rate_kg_s=qo,dv_km_s=dv/1000,ve_km_s=ve/1000,
            propellant_per_delivered_kg=prop,electric_energy_MJ_kg=e/1e6,
            transport_power_PW=qn*e/1e15,
            oxygen_power_24kWh_PW=qo*24.3*3.6e6/1e15))
    savecsv(out/'logistics.csv',logistics)
    routes=[hohmann(a) for a in [2.77,5.2,9.58,30.1,39.5]]
    savecsv(out/'transfer_routes.csv',routes)
    lifetimes=[]
    for life in [1e5,1e6,1e7,1e8,1e9,1e10]:
        lifetimes.append(dict(depletion_timescale_yr=life,loss_kg_s=mass/life/YR,
            build_throughput_fraction=500/life,inventory_replacements_per_Gyr=1e9/life))
    savecsv(out/'maintenance.csv',lifetimes)
    # Independent irradiation-to-escape energy ceiling; no inference of T_exo.
    energy=[]
    for rr in [1.5,2,3,5]:
      for eta in [.01,.1,.3]:
        q=eta*np.pi*(rr*R)**3*4.64e-3/GM
        energy.append(dict(absorption_radius_R=rr,efficiency=eta,unfiltered_kg_s=q,
          transmission_for_10kg_s=min(1,10/q),transmission_for_100kg_s=min(1,100/q)))
    savecsv(out/'energy_limited_sensitivity.csv',energy)
    # Optical screen geometry and thrust envelope, not a solved trajectory.
    shields=[]
    for rprot in [2,3,5]:
      for d in [60e6,1.5e9]:
       rad=R*rprot+d*.00465;A=np.pi*rad**2
       for areal in [.001,.01,.1]:
        m=A*areal
        shields.append(dict(protected_radius_R=rprot,distance_m=d,areal_kg_m2=areal,
          screen_diameter_km=2*rad/1000,screen_area_m2=A,mass_kg=m,
          replacement_100yr_kg_s=m/(100*YR),
          full_absorber_solar_acceleration=1361/(299792458*areal),
          thrust_1e4ms2_N=m*1e-4,propellant_ve50kms_kg_s=m*1e-4/50000))
    savecsv(out/'shields.csv',shields)
    # Industrial growth: solve cumulative output, with 100 quadrillion kg/yr cap.
    grows=[]
    for q0 in [1e6,1e9,1e12]:
      for start in [100,150,200]:
       g=brentq(lambda rate:growth(rate,q0,start)[0]-mass,.00001,.2)
       total,peak=growth(g,q0,start)
       grows.append(dict(initial_capacity_kg_yr=q0,start_year=start,
         continuous_growth_rate=g,equivalent_annual_growth=np.expm1(g),
         doubling_time_yr=np.log(2)/g,peak_rate_kg_yr=peak,total_mass_kg=total))
    savecsv(out/'growth.csv',grows)
    # Independent basic thermal response proxy, linear energy-balance equation.
    thermal=[]
    for colfrac in [.1,1]:
      for water in [0,10,50]:
        c=colfrac*CP*baseline['column_kg_m2']+water*1000*4180
        for b in [1,2,4]:
          amp=200/np.sqrt(b*b+(2*np.pi/(29.53059*86400)*c)**2)
          thermal.append(dict(coupled_atmosphere_fraction=colfrac,ocean_m=water,
             OLR_slope_W_m2_K=b,forcing_amplitude_W_m2=200,
             response_amplitude_K=amp,radiative_timescale_days=c/b/86400))
    savecsv(out/'thermal_response.csv',thermal)
    # Verification: units, independent exact cases, convergence and mass integral.
    b2,p2=atmosphere(points=60001)
    relmass=abs(baseline['mass_kg']/b2['mass_kg']-1)
    relloss=abs(baseline['loss_kg_s']/b2['loss_kg_s']-1)
    direct=4*np.pi*np.trapezoid(profile['r']**2*profile['p']*MU/(RG*profile['t']),profile['r'])
    intdiff=abs(direct/mass-1)
    iso,ip=atmosphere(te=200,isothermal=True)
    exact=1/(1/R-(RG/MU)*200/GM*ip['x'])
    isodiff=float(np.max(abs(exact/ip['r']-1)))
    # Endpoint interpolated using ln(Kn), excluded from analytic grid comparison.
    isodiff=float(np.max(abs(exact[:-1]/ip['r'][:-1]-1)))
    assert relmass<1e-5 and relloss<1e-3 and intdiff<1e-5 and isodiff<1e-8
    assert abs(hohmann(1)['total_endpoint_dv_km_s'])<1e-10
    assert mass>inventory['thin_mass_kg']
    verification=dict(mass_grid_relative_difference=relmass,loss_grid_relative_difference=relloss,
      independent_mass_integral_relative_difference=intdiff,isothermal_analytic_relative_error=isodiff,
      zero_transfer_dv_check=True,tests_passed=True)
    result=dict(baseline=baseline,tidal=tidal,diffusive=diffusive,inventory=inventory,verification=verification,
      routes=routes,growth=grows,logistics=logistics,energy=energy,
      lifetime=lifetimes,thermal=thermal,
      constants=dict(R_m=R,GM_m3_s2=GM,pressure_Pa=P0,surface_T_K=TS,
       nitrogen_mole_fraction=1-XO,oxygen_mole_fraction=XO,
       collision_cross_section_m2=SIGMA,cp_J_kg_K=CP))
    (out/'results.json').write_text(json.dumps(result,indent=2))
    # Figures all label the approximation.
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for tc in [150,180,210]:
      s,pf=atmosphere(tc,250,.1)
      ax[0].plot(pf['t'],pf['r']/R,label=f'Middle atmosphere {tc} K')
      ax[1].semilogx(pf['p'],pf['r']/R)
    ax[0].set(xlabel='Prescribed temperature (K)',ylabel='Radius / lunar radius',ylim=(1,7))
    ax[1].set(xlabel='Pressure (Pa)',ylabel='Radius / lunar radius',ylim=(1,7));ax[1].invert_xaxis()
    ax[0].legend(fontsize=8);fig.suptitle('Hydrostatic sensitivity profiles; exobase target 250 K')
    fig.tight_layout();fig.savefig(out/'atmosphere_profiles.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.3))
    for tc in [150,180,210]:
      ts=np.arange(210,401,5);vals=[]
      for te in ts:
        s,_=atmosphere(tc,te,.1,points=12001)
        vals.append(s.get('loss_kg_s',np.nan))
      ax.semilogy(ts,vals,label=f'Middle atmosphere {tc} K')
    ax.axhline(mass/(1e9*YR),color='black',ls='--',lw=1,label='One inventory / Gyr')
    ax.axhline(.1*mass/(1e9*YR),color='gray',ls=':',lw=1,label='10% inventory / Gyr')
    ax.set(xlabel='Asymptotic thermosphere temperature (K)',ylabel='Molecular Jeans escape (kg/s)',
      title='Prescribed-profile sensitivity, not a predicted atmospheric lifetime')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'escape_sensitivity.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.2))
    years=np.linspace(150,500,500)
    for g in [.02,.04,.05,.06]:
      cumulative=[growth(g,1e9,150,y)[0] for y in years]
      ax.semilogy(years,np.maximum(cumulative,1),label=f'{100*np.expm1(g):.1f}% annual growth')
    ax.axhline(mass,color='black',ls='--',label='Atmosphere inventory')
    ax.set(xlabel='Years after 2026',ylabel='Cumulative delivered gas (kg)',ylim=(1e9,1e20),
       title='Conditional capacity growth from 10^9 kg/year at project year 150')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'industrial_growth.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.2))
    tau=np.logspace(4,10,100);ax.loglog(tau,500/tau,color='#176a7a')
    ax.set(xlabel='Protected depletion timescale (years)',ylabel='Maintenance / construction mass throughput',
      title='Replacement capacity can be small while cumulative resource demand is large')
    ax.grid(which='major',alpha=.2);fig.tight_layout();fig.savefig(out/'maintenance_ratio.png');plt.close(fig)
    print(json.dumps(dict(baseline=baseline,inventory=inventory,verification=verification,
                         growth=grows[4]),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results')
    main(Path(p.parse_args().out))
