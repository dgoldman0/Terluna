"""Lunar protection design calculations. SI unless a field says otherwise.

These calculations close selected component budgets, not a coupled atmospheric,
plasma or flexible-structure model. See report.md for validity boundaries.
"""
from pathlib import Path
import csv, json, math, hashlib, re
import numpy as np
from scipy.optimize import differential_evolution

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
MU0=4*np.pi*1e-7; C=299792458.; SIG=5.670374419e-8
GME=3.986004418e14; GMM=4.902800118e12; GMS=1.32712440018e20
AU=149597870700.; RMOON=1737400.; AMOON=384400000.
NE=np.sqrt(GMS/AU**3); YEAR=365.25*86400.; SOLAR=1361.
ALPHA=695700000./AU; RP=3*RMOON
ATM=3.14134948e18; NITROGEN=2.52867008e18

def save_csv(name, rows):
    rows=list(rows)
    with (OUT/(name+'.csv')).open('w',newline='') as f:
        fields=list(dict.fromkeys(k for row in rows for k in row))
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def titanium_data():
    lines=(ROOT/'sources/TiO2_Siefke.yml').read_text().splitlines()
    rows=[]
    for line in lines:
        z=line.split()
        if len(z)==3:
            try:rows.append([float(x) for x in z])
            except ValueError:pass
    return np.array(rows)
TI=titanium_data()

def ti_nk(w_um):
    w=np.asarray(w_um)
    if w.min()<TI[0,0] or w.max()>TI[-1,0]:raise ValueError('TiO2 data range')
    return np.interp(w,TI[:,0],TI[:,1])+1j*np.interp(w,TI[:,0],TI[:,2])

def si_n(w_um):
    """Malitson dispersion; used within the 0.21--3.71 um measured interval.
    Silica absorption is neglected here; UV loss is dominated by TiO2.
    """
    w=np.asarray(w_um);v=w*w
    return np.sqrt(1+0.6961663*v/(v-0.0684043**2)+0.4079426*v/(v-0.1162414**2)+0.8974794*v/(v-9.896161**2))

def mix_n(n1,n2,f):
    """Lorentz-Lorenz effective-medium approximation, not measured mixed film."""
    l1=(n1*n1-1)/(n1*n1+2);l2=(n2*n2-1)/(n2*n2+2)
    l=(1-f)*l1+f*l2
    return np.sqrt((1+2*l)/(1-l))

def stack_rt(w,indices,thickness_um):
    """Stable normal-incidence Fresnel recursion, complex n uses +ik."""
    ns=[np.ones_like(w,dtype=complex)]+indices+[np.ones_like(w,dtype=complex)]
    r=(ns[-2]-ns[-1])/(ns[-2]+ns[-1]);t=2*ns[-2]/(ns[-2]+ns[-1])
    for j in range(len(indices)-1,-1,-1):
        phase=np.exp(1j*2*np.pi*ns[j+1]*thickness_um[j]/w)
        rj=(ns[j]-ns[j+1])/(ns[j]+ns[j+1]);tj=2*ns[j]/(ns[j]+ns[j+1])
        den=1+rj*r*phase**2
        t=tj*t*phase/den;r=(rj+r*phase**2)/den
    return abs(r)**2,abs(t)**2

def planck_weight(w):
    wm=w*1e-6
    return wm**-5/np.expm1(0.01438776877/(5772*wm))

def optical_stack(w, ar, silica_um=10., titanium_um=1.):
    s=si_n(w).astype(complex); ti=ti_nk(w)
    p=mix_n(np.ones_like(s),s,.5)
    m=mix_n(s,ti,.5)
    return stack_rt(w,[p,s,m,ti,m,s,p],[ar[0],ar[1],ar[2],titanium_um,ar[3],silica_um,ar[4]])

def xray_mu(energy, formula, molar_mass, density):
    """Independent atom estimate from CXRO f2. Validated study range >=50 eV."""
    f2=np.zeros_like(np.asarray(energy,dtype=float))
    for element,count in formula.items():
        a=np.loadtxt(ROOT/'sources'/f'{element}.nff',skiprows=1)
        f2 += count*np.exp(np.interp(np.log(energy),np.log(a[:,0]),np.log(a[:,2])))
    n=density/molar_mass*6.02214076e23
    wavelength=1.239841984e-6/np.asarray(energy)
    return 2*2.8179403262e-15*wavelength*n*f2

def spectra():
    w=np.linspace(.4,2.5,501);wt=planck_weight(w)
    def objective(ar):
        r,t=optical_stack(w,ar)
        return -float(np.trapezoid(t*wt,w)/np.trapezoid(wt,w))
    fit=differential_evolution(objective,[(.01,.6)]*5,seed=704,maxiter=220,popsize=10,tol=1e-7,polish=True)
    ar=fit.x
    wf=np.linspace(.4,2.5,7001);weight=planck_weight(wf)
    rr,tt=optical_stack(wf,ar)
    r0,t0=stack_rt(wf,[ti_nk(wf),si_n(wf).astype(complex)],[1.,10.])
    def avg(a,lo,hi):
        sel=(wf>=lo)&(wf<=hi)
        return float(np.trapezoid(a[sel]*weight[sel],wf[sel])/np.trapezoid(weight[sel],wf[sel]))
    save_csv('visible_stack',({'wavelength_um':float(w),'R':float(r),'T':float(t),'absorption':float(1-r-t)} for w,r,t in zip(wf,rr,tt)))
    uv=[]
    for nm in [121.6,150,180,200,250,300,350,400,500,1000]:
        nk=ti_nk(np.array([nm/1000]))[0];od=4*np.pi*nk.imag/(nm/1000)
        uv.append({'wavelength_nm':nm,'TiO2_n':float(nk.real),'TiO2_k':float(nk.imag),'TiO2_1um_absorption_only_T':float(np.exp(-od)),'thickness_um_for_absorption_T_1e_5':float(np.log(1e5)/(4*np.pi*nk.imag/(nm/1000))) if nk.imag>0 else None})
    save_csv('uv_absorption',uv)
    xr=[]
    for kev in [.05,.1,.2,.5,1,2,5,10,20,30]:
        e=np.array([1000*kev]); ms=xray_mu(e,{'si':1,'o':2},.0600843,2200)[0];mt=xray_mu(e,{'ti':1,'o':2},.079866,3900)[0]
        xr.append({'energy_keV':kev,'wavelength_nm':1.239841984/kev,'silica_attenuation_length_um':float(1e6/ms),'T_10um_silica_1um_titania':float(np.exp(-ms*10e-6-mt*1e-6)),'T_100um_silica_1um_titania':float(np.exp(-ms*100e-6-mt*1e-6))})
    save_csv('xray_absorption',xr)
    checks={'energy_balance_min_A':float(np.min(1-rr-tt)),'vacuum_layer_transmission_error':float(np.max(abs(stack_rt(wf,[np.ones_like(wf,dtype=complex)],[10])[1]-1)))}
    hsrs=np.loadtxt(ROOT/'sources/TSIS1_HSRS_stride100.csv',delimiter=',',skiprows=1)
    sw=hsrs[:,0]/1000;flux=hsrs[:,1]
    sel=(sw>=.4)&(sw<=2.5);sr,st=optical_stack(sw[sel],ar)
    solar_t=float(np.trapezoid(st*flux[sel],sw[sel])/np.trapezoid(flux[sel],sw[sel]))
    vis=(sw>=.4)&(sw<=.7);vr,vt=optical_stack(sw[vis],ar)
    visible_t=float(np.trapezoid(vt*flux[vis],sw[vis])/np.trapezoid(flux[vis],sw[vis]))
    reliable=(sw>=.21)&(sw<=2.73);_,rt=optical_stack(sw[reliable],ar)
    admitted=float(np.trapezoid(rt*flux[reliable],sw[reliable]*1000))
    # Full-spectrum bounds: outside measured/evaluated 210--2730 nm is allowed
    # anywhere from zero to unity; a broad bound, not a measured solar T.
    incident=float(np.trapezoid(flux[reliable],sw[reliable]*1000))
    solar={'TSIS_T_400_2500nm':solar_t,'TSIS_T_400_700nm':visible_t,'incident_210_2730nm_W_m2':incident,'transmitted_210_2730nm_W_m2':admitted,'full_1361W_m2_T_lower_bound':admitted/SOLAR,'full_1361W_m2_T_upper_bound':(admitted+SOLAR-incident)/SOLAR,'quadrature_note':'HSRS sampled every 0.1 nm, treated as a broadband approximation; optical fit used Planck weighting.'}
    save_csv('solar_weighted_summary',[solar])
    return {'AR_layers_um':ar.tolist(),'AR_index_models':['50% silica / 50% void','silica','50% silica / 50% titania','50% silica / 50% titania, buried','50% silica / 50% void, rear'], 'T_weighted_400_2500nm':avg(tt,.4,2.5),'T_weighted_400_700nm':avg(tt,.4,.7),'R_weighted_400_2500nm':avg(rr,.4,2.5),'uncoated_T_weighted_400_2500nm':avg(t0,.4,2.5),'checks':checks,'solar':solar,'uv':uv,'xray':xr}

def solar_tide(r):
    sun=np.array([AU,0.,0.]);delta=sun-r
    return GMS*(delta/np.linalg.norm(delta,axis=-1)[...,None]**3-sun/AU**3)

def holding_acceleration(d, a=AMOON, count=2048):
    """Instantaneous point-mass gravity along a circular lunar phase sweep.
    Screen r = Moon r + d*Sun unit vector. Lunar acceleration is computed from
    gravity, avoiding an inconsistent assumed circular acceleration in solar tide.
    Sun direction's annual centripetal acceleration is included. No SRP here.
    """
    ph=np.arange(count)*2*np.pi/count
    rm=np.stack([a*np.cos(ph),a*np.sin(ph),np.zeros(count)],axis=1)
    rs=rm+np.array([d,0.,0.])
    u= -GME*rm/a**3+GME*rs/np.linalg.norm(rs,axis=1)[:,None]**3
    u += np.array([GMM/d**2-d*NE**2,0.,0.])
    u += solar_tide(rm)-solar_tide(rs)
    mag=np.linalg.norm(u,axis=1)
    return ph,u,mag

def mass_closure(base, mean_a, peak_a, ve=30000., kappa=300., eta=.7, buffer_days=7, cant_deg=45., storage_hours=0., storage_J_kg=1e6):
    # kappa includes electrical supply, processing, thrusters and their thermal
    # equipment; it is a system requirement, not current proven capability.
    axial=np.cos(np.deg2rad(cant_deg))
    p_fraction=peak_a*ve/(2*eta*kappa*axial)
    b_fraction=mean_a*buffer_days*86400/(ve*axial)
    storage_fraction=peak_a*ve/(2*eta*axial)*storage_hours*3600/storage_J_kg
    den=1-p_fraction-b_fraction-storage_fraction
    if den<=0:return {'closed':False,'power_mass_fraction':p_fraction,'fuel_buffer_fraction':b_fraction}
    mass=base/den;mdot=mass*mean_a/(ve*axial);p=mdot*ve**2/(2*eta);pmax=mass*peak_a*ve/(2*eta*axial)
    return {'closed':True,'cant_deg':cant_deg,'total_mass_kg':mass,'power_mass_fraction':p_fraction,'fuel_buffer_fraction':b_fraction,'energy_storage_mass_kg':mass*storage_fraction,'propellant_kg_s':mdot,'mean_power_W':p,'peak_power_W':pmax,'power_system_mass_kg':pmax/kappa,'fuel_buffer_mass_kg':mdot*buffer_days*86400,'Gyr_propellant_kg':mdot*YEAR*1e9}

def orbital_designs():
    rows=[]
    for d in np.arange(20000,180001,2000)*1000.:
        _,_,aa=holding_acceleration(d)
        area=np.pi*(RP+ALPHA*d)**2
        # Conservative SRP force allowance, 0.15 radiation-pressure coefficient
        # acting on 50 g/m2 would give 1.36e-5 m/s2 before added hardware.
        mean_a=float(aa.mean())+2e-5; peak_a=(float(aa.max())+2e-5)*1.25
        z=mass_closure(area*.05,mean_a,peak_a)
        rows.append({'distance_km':d/1000,'area_m2':area,'mean_a_m_s2':float(aa.mean()),'max_a_m_s2':float(aa.max()),'mean_design_a_m_s2':mean_a,'peak_design_a_m_s2':peak_a,**z})
    save_csv('near_lunar_distance_sweep',rows)
    best=min(rows,key=lambda x:x.get('propellant_kg_s',float('inf')))
    d=best['distance_km']*1000;ph,u,aa=holding_acceleration(d)
    save_csv('near_lunar_phase',({'phase_deg':float(p*180/np.pi),'ax_m_s2':float(v[0]),'ay_m_s2':float(v[1]),'norm_a_m_s2':float(a)} for p,v,a in zip(ph,u,aa)))
    ecc=[]
    for am in [363300000.,AMOON,405500000.]:
        _,_,a=holding_acceleration(d,am)
        ecc.append({'lunar_distance_km':am/1000,'mean_a':float(a.mean()),'max_a':float(a.max())})
    save_csv('eccentric_distance_sensitivity',ecc)
    cases=[]
    for sigma in [.025,.05,.1,.25]:
        for ve in [10000.,30000.,50000.,100000.]:
            for kappa in [50.,100.,300.,1000.]:
                z=mass_closure(best['area_m2']*sigma,best['mean_design_a_m_s2'],best['peak_design_a_m_s2'],ve,kappa)
                cases.append({'areal_mass_kg_m2':sigma,'exhaust_m_s':ve,'specific_power_W_kg':kappa,**z})
    # Fill failed cases to a uniform schema.
    keys=list(dict.fromkeys(k for r in cases for k in r))
    save_csv('propulsion_sensitivity',({k:r.get(k) for k in keys} for r in cases))
    # A fixed Sun-Earth L1 hub plus a Moon-following optical screen at 1.5M km from Moon.
    _,_,far_a=holding_acceleration(1.5e9)
    far_area=np.pi*(RP+ALPHA*1.5e9)**2
    far={'distance_km':1.5e6,'area_m2':far_area,'mean_a':float(far_a.mean()),'max_a':float(far_a.max()),**mass_closure(far_area*.05,float(far_a.mean()),float(far_a.max())*1.25)}
    # Broad fixed-aperture comparison only: not a solved orbit/structure.
    d_earth=1.5e9;limit=405500000.;margin=RP+ALPHA*(d_earth+limit)
    radius=limit+margin
    band_width=2*(limit+margin);band_height=2*(limit*np.sin(np.deg2rad(5.15))+margin)
    broad={'disk_radius_km':radius/1000,'disk_area_m2':np.pi*radius**2,'band_width_km':band_width/1000,'band_height_km':band_height/1000,'band_area_m2':band_width*band_height,'at_50g_disk_mass':.05*np.pi*radius**2,'at_50g_band_mass':.05*band_width*band_height}
    annual_A=best['area_m2']/150;gaps=[]
    for clear in [1e-2,1e-3,1e-4,1e-5]:
        gaps.append({'uncovered_fraction':clear,'system_T_if_material_T_1e_6':clear+(1-clear)*1e-6,'area_uncovered_m2':clear*best['area_m2']})
    save_csv('leakage_budget',gaps)
    sizes=[]
    for mult in [2.,3.,4.,5.]:
        area=np.pi*(mult*RMOON+ALPHA*d)**2
        sizes.append({'protected_radius_RMoon':mult,'screen_diameter_km':2*np.sqrt(area/np.pi)/1000,'area_m2':area,**mass_closure(area*.05,best['mean_design_a_m_s2'],best['peak_design_a_m_s2'])})
    save_csv('protected_radius_sensitivity',sizes)
    backup=[]
    for hours in [0.,1.,2.,4.]:
        backup.append({'storage_hours_at_peak_power':hours,'storage_J_kg':1e6,**mass_closure(best['area_m2']*.05,best['mean_design_a_m_s2'],best['peak_design_a_m_s2'],storage_hours=hours)})
    save_csv('energy_storage_sensitivity',backup)
    return {'reference':best,'far_tracking':far,'broad_aperture':broad,'eccentricity_sensitivity':ecc,'protected_radius_sensitivity':sizes,'energy_storage_sensitivity':backup,'annual_filter_area_150yr_m2':annual_A,'minimum_specific_power_for_reference_ve30km_s_no_fuel_W_kg':best['peak_design_a_m_s2']*30000/(2*.7*np.cos(np.deg2rad(45)))}

def loop_model(radius, ampturn, bundle_radius=1000., Je=1e8, specific_strength=1e6):
    log=np.log(8*radius/bundle_radius)
    inductance=MU0*radius*(log-2)
    u=.5*inductance*ampturn**2
    tension=MU0*ampturn**2/(4*np.pi)*(log-.75)
    return {'radius_km':radius/1000,'ampere_turns':ampturn,'moment_A_m2':ampturn*np.pi*radius**2,'center_B_T':MU0*ampturn/(2*radius),'equivalent_inductance_H':inductance,'stored_energy_J':u,'hoop_tension_N':tension,'ideal_support_mass_kg':2*np.pi*radius*tension/specific_strength,'conductor_mass_kg':2*np.pi*radius*ampturn/Je*6000,'terminal_current_A':50000.,'turn_count_at_50kA':ampturn/50000,'conductor_length_at_50kA_m':2*np.pi*radius*ampturn/50000,'rough_bundle_self_field_T':MU0*ampturn/(2*np.pi*bundle_radius)}

def magnetic_designs():
    pressure=[]
    moment=1.5e21
    for p in [2e-9,20e-9,100e-9]:
        b=np.sqrt(2*MU0*p);standoff=(1e-7*moment/b)**(1/3)
        pressure.append({'solar_wind_pressure_nPa':p*1e9,'balance_field_nT':b*1e9,'dipole_equatorial_standoff_km':standoff/1000,'uniform_sphere_field_energy_at_RP_J':p*4*np.pi*RP**3/3})
    save_csv('wind_balance',pressure)
    anchors=[]
    for a in [50000.,100000.,200000.,500000.]:
        turns=moment/(4*np.pi*a*a)
        z=loop_model(a,turns,bundle_radius=min(1000.,a/50))
        anchors.append({'station_count':4,**z,'total_four_station_conductor_support_kg':4*(z['conductor_mass_kg']+z['ideal_support_mass_kg'])})
    save_csv('anchor_magnets',anchors)
    reference=anchors[1]
    rigid=[]
    for energy in [.001,.01,.1,1.,10.]:
        rig=np.sqrt(energy*(energy+2*.9382720813))
        rigid.append({'proton_kinetic_GeV':energy,'rigidity_GV':float(rig),'BL_for_one_radian_T_m':float(rig/.299792458),'radius_at_71nT_km':float(rig/.299792458/(71e-9)/1000),'radius_at_50uT_km':float(rig/.299792458/(50e-6)/1000)})
    save_csv('rigidity',rigid)
    volume=[]
    for radius in [RP,2e6,4e8]:
        for b in [71e-9,5e-6,50e-6]:
            vol=4*np.pi*radius**3/3;u=b*b/(2*MU0)*vol
            volume.append({'radius_km':radius/1000,'B_T':b,'uniform_field_energy_J':u,'dipole_equatorial_vertical_cutoff_GV_diagnostic':.299792458*b*radius/4,'energy_over_1MJ_kg_support_scale_kg':u/1e6,'rebuild_power_if_energy_lost_daily_W':u/86400,'rebuild_power_if_energy_lost_yearly_W':u/YEAR})
    save_csv('field_volume',volume)
    optimized=[]
    for radius in [2e6,4e8]:
        for cutoff in [1.,5.,20.]:
            b=4*cutoff/(.299792458*radius);u=b*b/(2*MU0)*4*np.pi*radius**3/3
            optimized.append({'radius_km':radius/1000,'equatorial_vertical_cutoff_GV':cutoff,'required_equatorial_B_T':b,'uniform_field_energy_scale_J':u,'energy_over_1MJ_kg_support_scale_kg':u/1e6})
    save_csv('fixed_cutoff_scaling',optimized)
    corridor_loop=loop_model(2e6,5e-5*2*2e6/MU0)
    global_dipole_moment=4*20/(.299792458*4e8)*(4e8)**3/1e-7
    for n in [1e3,1e4,1e5]:
        # Illustrative moment sum, ignores mutual inductance and geometric cancellation.
        z=loop_model(1e6,global_dipole_moment/(n*np.pi*1e12))
        z['node_count']=n;z['summed_conductor_support_kg']=n*(z['conductor_mass_kg']+z['ideal_support_mass_kg']);z['summed_self_energy_J']=n*z['stored_energy_J'];optimized.append(z)
    return {'pressure':pressure,'anchor_reference':reference,'rigidity':rigid,'corridor_loop_50uT_isolated_center':corridor_loop,'corridor_200_station_mass_before_margin_kg':200*(corridor_loop['conductor_mass_kg']+corridor_loop['ideal_support_mass_kg']),'global_moment_for_equatorial_20GV_at400000km':global_dipole_moment,'distributed_global_moment_examples':optimized[-3:],'field_at_surface_equator_T':1e-7*moment/RMOON**3,'field_at_3R_equator_T':1e-7*moment/RP**3,'distributed_surface_geometry':anchor_field_geometry(moment)}

def anchor_field_geometry(moment):
    """Four 100 km loops approximated as dipoles, two at each lunar pole.
    Centers 200 km from spin axis. All moments parallel +z; actual terrain and
    coil supports are separate work. R_loop / nearest evaluation distance <.04.
    """
    z=np.sqrt(RMOON**2-(2e5)**2)
    centers=np.array([[2e5,0,z],[-2e5,0,z],[2e5,0,-z],[-2e5,0,-z]])
    def field(p):
        p=np.asarray(p);v=p[:,None,:]-centers[None,:,:];r=np.linalg.norm(v,axis=-1)
        m=np.array([0.,0.,moment/4])
        return np.sum(1e-7*(3*v*np.sum(v*m,axis=-1)[...,None]/r[...,None]**5-m/r[...,None]**3),axis=1)
    th=np.linspace(0,np.pi,361);az=np.arange(72)*2*np.pi/72
    dirs=np.array([[np.sin(t)*np.cos(a),np.sin(t)*np.sin(a),np.cos(t)] for t in th for a in az])
    b=np.linalg.norm(field(dirs*RP),axis=1)
    rows=[]
    for p in [2e-9,20e-9,100e-9]:
        needed=np.sqrt(2*MU0*p)
        lo=np.ones(len(dirs))*RP;hi=np.ones(len(dirs))*3e7
        for _ in range(40):
            mid=(lo+hi)/2;bm=np.linalg.norm(field(dirs*mid[:,None]),axis=1)
            lo=np.where(bm>needed,mid,lo);hi=np.where(bm>needed,hi,mid)
        radii=(lo+hi)/2
        rows.append({'pressure_nPa':p*1e9,'minimum_scalar_balance_radius_km':float(radii.min()/1000),'maximum_scalar_balance_radius_km':float(radii.max()/1000)})
    save_csv('distributed_anchor_field',rows)
    cold_area=4*(2*np.pi*1e5)*(2*np.pi*1e3)
    cryo=[]
    for heat in [.01,.05,.1]:
        for cop in [.002,.01]:
            cryo.append({'cold_surface_area_m2':cold_area,'heat_leak_W_m2':heat,'COP':cop,'cold_load_W':cold_area*heat,'electrical_W':cold_area*heat/cop})
    save_csv('magnet_cryogenics',cryo)
    return {'centers_m':centers.tolist(),'moment_total_A_m2':moment,'minimum_vacuum_B_at3R_T':float(b.min()),'maximum_vacuum_B_at3R_T':float(b.max()),'scalar_pressure_balance_radii':rows,'cryogenic_sensitivity':cryo}

def atmospheric_interfaces(orb):
    """Diagnostic energy ceilings, NOT an escape prescription or climate model."""
    ref=orb['reference'];binding=GMM/RP
    diffuse_flux=1000*1e10/4*(6.62607015e-34*C/121.6e-9)
    diffuse_power=4*np.pi*RP**2*diffuse_flux
    plume_power=.5*ref['propellant_kg_s']*30000.**2
    rows=[]
    for loss in [1.,10.,100.,1000.]:
        allowed=loss*binding/.1
        rows.append({'allowance_kg_s':loss,'efficiency_assumption':.1,'permitted_incident_power_W':allowed,'permitted_fraction_of_total_thruster_kinetic_power':allowed/plume_power,'Gyr_mass_kg':loss*YEAR*1e9})
    save_csv('plume_energy_allowances',rows)
    thermal=[]
    for absorbed in [10.,50.,100.]:
        for emissivity in [.2,.8]:
            temp=(absorbed/(2*emissivity*SIG))**.25
            earthward=.5*absorbed*ref['area_m2']/(np.pi*(ref['distance_km']*1000)**2)
            thermal.append({'absorbed_solar_W_m2':absorbed,'two_face_emissivity':emissivity,'equilibrium_K':temp,'screen_IR_irradiance_at_Moon_center_plane_W_m2':earthward})
    save_csv('screen_thermal_sensitivity',thermal)
    return {'binding_J_kg_at3R':binding,'isotropic_1000R_Lyman_alpha_flux_per_surface_W_m2':diffuse_flux,'isotropic_1000R_incident_power_at3R_W':diffuse_power,'energy_only_loss_ceiling_at_eta_01_kg_s':.1*diffuse_power/binding,'total_thruster_exhaust_kinetic_power_W':plume_power,'Moon_protected_angular_radius_from_screen_deg':float(np.rad2deg(np.arcsin(RP/(ref['distance_km']*1000)))),'plume_allowances':rows}

def maintenance(orb,mag):
    reference=orb['reference'];base=reference['area_m2']*.05
    rows=[]
    for name,mass in [('near_filter',base),('all_near_hardware_excluding_fuel',reference['total_mass_kg']-reference['fuel_buffer_mass_kg']),('lunar_magnets_10x_floor',mag['anchor_reference']['total_four_station_conductor_support_kg']*10),('broad_L1_band',orb['broad_aperture']['at_50g_band_mass']),('broad_L1_disk',orb['broad_aperture']['at_50g_disk_mass'])]:
        for life in [20.,100.,1000.]:
            gross=mass/(life*YEAR)
            rows.append({'component':name,'mass_kg':mass,'life_years':life,'gross_replacement_kg_s':gross,'fresh_at_99percent_recovery_kg_s':gross*.01,'fresh_at_99_9percent_recovery_kg_s':gross*.001,'Gyr_fresh_at_99_9percent_recovery_kg':gross*.001*YEAR*1e9,'rebuild_power_at_100MJ_kg_W':gross*1e8,'fraction_of_500yr_atmosphere_construction_rate':gross/(ATM/(500*YEAR))})
    save_csv('maintenance',rows)
    delivery=reference['propellant_kg_s']
    return {'reference_propellant_atmos_construction_rate_fraction':delivery/(ATM/(500*YEAR)),'reference_propellant_Gyr_Moon_masses':delivery*YEAR*1e9/7.342e22,'propellant_production_transport_power_at_100MJ_kg_W':delivery*1e8,'reference_filter_build_power_150yr_at_100MJ_kg_W':base*1e8/(150*YEAR),'solar_generation_area_at_300W_m2_m2':reference['peak_power_W']/300,'thruster_waste_radiator_area_600K_eps09_m2':(.3*reference['peak_power_W'])/(.9*SIG*600**4),'filter_replacement_100yr_kg_s':base/(100*YEAR)}

def main():
    op=spectra();orb=orbital_designs();mag=magnetic_designs();maint=maintenance(orb,mag)
    result={'constants':{'atmosphere_kg':ATM,'N2_kg':NITROGEN,'construction_average_500yr_kg_s':ATM/(500*YEAR),'nitrogen_300yr_kg_s':NITROGEN/(300*YEAR),'initial_transport_power_reference_W':[3.3e16,2.35e17]},'optics':op,'orbits':orb,'magnetics':mag,'maintenance':maint,'atmospheric_interfaces':atmospheric_interfaces(orb)}
    result['input_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'sources').glob('*')) if p.suffix in ['.nff','.yml','.csv']}
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
