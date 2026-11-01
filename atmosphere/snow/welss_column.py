"""Declared-source Welss ob_hp0 deposition/settling comparison.

Independent continuous-mass characteristics, without collisions or a vapour cap.
McSnow source constants are retained as reference inputs, not new Terluna defaults.
The 2024 public revision is not proved to have produced the 2023 author archive.

Geometry, atmosphere and deposition equations derived from McSnow:
Copyright (c) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss.
Authors: Slavko Brdar, Axel Seifert, Christoph Siewert, Markus Karrer,
Jan-Niklas Welss. SPDX-License-Identifier: BSD-3-Clause.
See references/mcsnow_bsd_3_clause.txt and welss_column_sources.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.io import netcdf_file
from scipy.special import roots_genlaguerre

from atmosphere.snow import boehm

HERE = Path(__file__).resolve().parent
REVISION = boehm.SOURCE_REVISION
SCHEMA = 'terluna.atmosphere.welss-column/1'
ARCHIVE_SHA256 = '119871a2d269a90321a7fba21658b6f217aa2f8641dbda37f6cd9ae6ec41a304'
SOURCE_CONSTANTS = dict(ps=101325., grav=9.81, tlr=.0062, Rd=287.102,
                        Rv=461.401, c_pa=1005., rhoi=919., T_3=273.15,
                        e_3=610.78, A_i=21.8745584, A_w=17.2693882,
                        B_i=7.66, B_w=35.86)
PI = math.pi
RHOI = SOURCE_CONSTANTS['rhoi']
ALPHA = .0028 * 10.**(2.*2.1 - 3.)
BETA = 2.1
AREA_ALPHA = .2285 * 10.**(2.*1.88 - 4.)
AREA_BETA = 1.88
D_THRESHOLD = (RHOI * PI/6 / ALPHA)**(1/(BETA-3.))
M_THRESHOLD = PI/6 * RHOI * D_THRESHOLD**3
# McSnow declares these two constants with unkinded, default-REAL literals.
SUTHERLAND_A = float(np.float32(1.458e-6))
SUTHERLAND_T = float(np.float32(110.4))
SOURCE_RATE = 16.5  # physical particles m^-3 s^-1
SOURCE_BOTTOM, SOURCE_TOP = 2700., 3000.
SEED_MEAN = 4.8e-13
CELL_WIDTH = 20.
DOMAIN_TOP = 5000.
END_TIME = 36000.
MOMENTS = ('Fn', 'Fm', 'Nd', 'Md', 'deposition')


def saturation(T):
    T = np.asarray(T, dtype=float)
    c = SOURCE_CONSTANTS
    return (c['e_3']*np.exp(c['A_i']*(T-c['T_3'])/(T-c['B_i'])),
            c['e_3']*np.exp(c['A_w']*(T-c['T_3'])/(T-c['B_w'])))


def viscosity(T):
    return SUTHERLAND_A*np.asarray(T)**1.5/(np.asarray(T)+SUTHERLAND_T)


def ideal_atmosphere(z):
    """Pinned init_atmo_ideal for h1=0,h2=600,ssat=.05,lwc0=.0002."""
    z = np.asarray(z, dtype=float)
    c = SOURCE_CONSTANTS
    T = c['T_3']-c['tlr']*z
    ice, water = saturation(T)
    pv = np.where(z > 600., 1.05*ice, water)
    # Preserve source pv_srf, which uses local psatw despite its variable name.
    pd = (c['ps']-water)*np.exp(c['grav']/(c['Rd']*c['tlr'])*np.log(T/c['T_3']))
    p = pd+pv
    rho = pd/(c['Rd']*T)+pv/(c['Rv']*T)
    qc = np.where((z >= 0.) & (z <= 600.), .0002/rho, 0.)
    qv = 1./(1.+(p/pv-1.)*c['Rv']/c['Rd'])
    return dict(T=T,p=p,pv=pv,rho=rho,eta=viscosity(T),qc=qc,qv=qv,
                psati=ice,psatw=water)


def atmosphere(z):
    """Pinned current linear interpolation and water-saturation adjustment."""
    z = np.asarray(z, dtype=float)
    if not np.isfinite(z).all() or (z < 0).any() or (z > DOMAIN_TOP).any():
        raise ValueError('height must be within the declared 0–5000 m domain')
    index = np.clip(np.floor((DOMAIN_TOP-z)/CELL_WIDTH).astype(int),0,249)
    center = DOMAIN_TOP-(index+.5)*CELL_WIDTH
    other = index+np.where(z > center,-1,1)
    other = np.where(index == 0,1,np.where(index == 249,248,other))
    center2 = DOMAIN_TOP-(other+.5)*CELL_WIDTH
    a,b = ideal_atmosphere(center),ideal_atmosphere(center2)
    f = (z-center)/(center2-center)
    T = a['T']+(b['T']-a['T'])*f
    p = a['p']+(b['p']-a['p'])*f
    total = a['qv']+a['qc']+(b['qv']+b['qc']-a['qv']-a['qc'])*f
    c = SOURCE_CONSTANTS
    ice,water = saturation(T)
    qsat = c['Rd']/c['Rv']*water/(p-(1-c['Rd']/c['Rv'])*water)
    qv = np.minimum(total,qsat)
    pv = c['Rv']/c['Rd']*qv*p/(1+(c['Rv']/c['Rd']-1)*qv)
    rho = (p-pv)/(c['Rd']*T)+pv/(c['Rv']*T)
    return dict(T=T,p=p,pv=pv,rho=rho,eta=viscosity(T),qv=qv,
                qc=np.maximum(total-qsat,0.),psati=ice,psatw=water)


def geometry(mass):
    """Pure-ice, habit0, agggeo1 geometry and the literal source C branch."""
    m = np.asarray(mass,dtype=float)
    if not np.isfinite(m).all() or (m <= 0).any():
        raise ValueError('finite positive particle masses required')
    solid = (6/PI*m*(1/RHOI))**(1/3)
    solid_area = PI/4*solid**2
    D = np.where(m < M_THRESHOLD, solid, (m*(1/ALPHA))**(1/BETA))
    A = np.where(m < M_THRESHOLD, solid_area, np.maximum(AREA_ALPHA*D**AREA_BETA,solid_area))
    void = PI/6*D**3-m*(1/RHOI)
    # The source uses the sign of this residual even for compact spheres.
    # Its floating-point behaviour is checked independently by Fortran fixtures.
    C = D*np.where(void > 0.,1/PI,.5)
    return dict(D=D,A=A,q=A/(PI/4*D**2),C=C,void=void)


def microphysics(mass,state):
    """Return source geometry, terminal speed and unbounded-vapour growth."""
    m = np.asarray(mass,dtype=float)
    g = geometry(m)
    props = boehm.particle_properties(g['D'],m,SOURCE_CONSTANTS['grav'],state,g['q'])
    T,p,eta,rho = [np.asarray(state[k]) for k in ('T','p','eta','rho')]
    c = SOURCE_CONSTANTS
    Dv = .211e-4*(T/c['T_3'])**1.94*c['ps']/p
    Ka = 418.4e-5*(5.69+.017*(T-c['T_3']))
    L = 2834100.-(T-c['T_3'])*(290.+(T-c['T_3'])*4.)
    def ventilation(X):
        return np.where(X < 1.,1.+.14*X**2,.86+.28*X)
    fv = ventilation(np.sqrt(props['Re'])*(eta/(rho*Dv))**(1/3))
    fh = ventilation(np.sqrt(props['Re'])*(eta*c['c_pa']/Ka)**(1/3))
    Dv_eff,Ka_eff = Dv*fv,Ka*fh
    ice = np.asarray(state['psati'])
    Hi = (L/(c['Rv']*T)-1.)*(L*Dv_eff*ice)/(Ka_eff*c['Rv']*T**2)
    growth = 4*PI*g['C']*Dv_eff*(np.asarray(state['pv'])-ice)/(c['Rv']*T*(1+Hi))
    return dict(**g,v=props['v'],Re=props['Re'],growth=growth,
                Dv=Dv,Ka=Ka,L=L,fv=fv,fh=fh)


def seed_quadrature(mass_nodes=16,height_nodes_per_cell=2):
    """Gamma(shape2,mean4.8e-13) × uniform source, including the full mass tail."""
    if mass_nodes < 2 or height_nodes_per_cell < 1:
        raise ValueError('at least two mass nodes and one height node required')
    x,mw = roots_genlaguerre(mass_nodes,1.)
    z,zw = np.polynomial.legendre.leggauss(height_nodes_per_cell)
    edges = np.arange(SOURCE_BOTTOM,SOURCE_TOP,CELL_WIDTH)
    heights = (edges[:,None]+CELL_WIDTH*(z+1)/2).ravel()
    height_weights = np.tile(CELL_WIDTH*zw/2,edges.size)
    return x*(SEED_MEAN/2),mw,heights,height_weights


def arrival_weight(age,mode):
    """Population weight for continuous injection starting at t=0."""
    age = np.asarray(age)
    if mode == 'last5h':
        return np.clip((END_TIME-age)/(END_TIME-18000.),0.,1.)
    if mode == 'four_snapshots':
        times = END_TIME-np.arange(3,-1,-1)*600.
        return np.mean(age[...,None] <= times,axis=-1)
    if mode == 'final':
        return (age <= END_TIME).astype(float)
    raise ValueError('unknown averaging mode')


def arrival_primitive(age,mode):
    """Integral of the arrival weight from zero to a nonnegative particle age."""
    a = np.maximum(np.asarray(age,dtype=float),0.)
    if mode == 'last5h':
        b = np.clip(a-18000.,0.,18000.)
        return np.minimum(a,18000.)+b-b*b/36000.
    if mode == 'four_snapshots':
        times = END_TIME-np.arange(3,-1,-1)*600.
        return np.mean(np.minimum(a[...,None],times),axis=-1)
    if mode == 'final':
        return np.minimum(a,END_TIME)
    raise ValueError('unknown averaging mode')


def trajectory(source_height,masses,*,rtol=1e-7,max_step_m=20.,deadline=None,rates=None):
    """Continuous log-mass and age versus descending height; masses stay positive."""
    masses = np.asarray(masses,dtype=float)
    if (masses <= 0).any() or not 0 < source_height <= DOMAIN_TOP:
        raise ValueError('positive masses and a source height in the domain required')
    n = masses.size
    rates = rates or (lambda z,m: microphysics(m,atmosphere(z)))
    def rhs(z,y):
        if deadline is not None and time.monotonic() > deadline:
            raise TimeoutError('Welss characteristic exceeded wall-time limit')
        m = np.exp(y[:n]);r = rates(z,m)
        v,G = np.asarray(r['v']),np.asarray(r['growth'])
        if (v <= 0).any() or not np.isfinite(v).all() or not np.isfinite(G).all():
            raise FloatingPointError('invalid trajectory rates')
        return np.r_[-G/(m*v),-np.broadcast_to(1/v,m.shape)]
    result = solve_ivp(rhs,(source_height,0.),np.r_[np.log(masses),np.zeros(n)],
                       method='DOP853',rtol=rtol,atol=1e-9,max_step=max_step_m,dense_output=True)
    if not result.success:
        raise RuntimeError(result.message)
    return result


def integrate(*,mass_nodes=16,height_nodes_per_cell=2,bin_nodes=3,rtol=1e-7,
              max_step_m=20.,wall_limit_s=600.,progress=False):
    """Source-quadrature expectation, binned like McSnow write_hprofiles."""
    started = time.monotonic();deadline = started+wall_limit_s
    masses,mw,heights,zw = seed_quadrature(mass_nodes,height_nodes_per_cell)
    q,qw = np.polynomial.legendre.leggauss(bin_nodes)
    modes = ('last5h','four_snapshots','final')
    profiles = {mode:np.zeros((len(MOMENTS),251)) for mode in modes}
    evaluations = 0;ground_ages=[];ground_masses=[]
    for index,(z0,zweight) in enumerate(zip(heights,zw)):
        sol = trajectory(z0,masses,rtol=rtol,max_step_m=max_step_m,deadline=deadline)
        evaluations += sol.nfev
        lower = np.arange(0.,z0,CELL_WIDTH)
        upper = np.minimum(lower+CELL_WIDTH,z0)
        zz = lower[:,None]+(upper-lower)[:,None]*(q+1)/2
        spatial_weights = (upper-lower)[:,None]*qw/(2*CELL_WIDTH)
        sampled = sol.sol(zz.ravel())
        mm,age = np.exp(sampled[:mass_nodes]),sampled[mass_nodes:]
        rates = microphysics(mm,atmosphere(zz.ravel()))
        v = rates['v']
        values = np.asarray([np.ones_like(mm),mm,1/v,mm/v,rates['growth']/v])
        lower_age = sol.sol(lower)[mass_nodes:]
        upper_age = sol.sol(upper)[mass_nodes:]
        for mode in modes:
            contributions = values*arrival_weight(age,mode)[None,:,:]*mw[None,:,None]
            point_moments = contributions.sum(axis=1).reshape(len(MOMENTS),len(lower),bin_nodes)
            cell_moments = np.sum(point_moments*spatial_weights[None,:,:],axis=-1)
            # da/dz=-1/v gives this moment exactly from endpoint ages, even
            # where tiny newborn crystals cross only a small fraction of a bin.
            cell_moments[2] = mw@(arrival_primitive(lower_age,mode)-arrival_primitive(upper_age,mode))/CELL_WIDTH
            profiles[mode][:,:len(lower)] += SOURCE_RATE*zweight*cell_moments
        at_ground = sol.sol(0.)
        ground_ages.extend(at_ground[mass_nodes:].tolist())
        ground_masses.extend(np.exp(at_ground[:mass_nodes]).tolist())
        if progress:
            print(f'Welss source height {index+1}/{len(heights)}; {time.monotonic()-started:.1f} s',flush=True)
    return dict(height=np.arange(251)*CELL_WIDTH,profiles=profiles,
                settings=dict(mass_nodes=mass_nodes,height_nodes_per_cell=height_nodes_per_cell,
                              bin_nodes=bin_nodes,rtol=rtol,max_step_m=max_step_m,
                              mass_node_range_kg=[float(masses.min()),float(masses.max())]),
                diagnostics=dict(evaluations=evaluations,wall_seconds=time.monotonic()-started,
                                 minimum_ground_age_s=min(ground_ages),maximum_ground_age_s=max(ground_ages),
                                 minimum_ground_mass_kg=min(ground_masses),maximum_ground_mass_kg=max(ground_masses),
                                 source_number_flux_m2_s=SOURCE_RATE*(SOURCE_TOP-SOURCE_BOTTOM),
                                 seed_number_normalization=float(mw.sum()),
                                 seed_mean_mass_kg=float(mw@masses)))


def lanczos5(values):
    """NCL6.6.2 five-point low-pass(fc=.2,nsigma1), opt0 missing endpoints."""
    fc=float(np.float32(.2));k=np.arange(-2,3,dtype=float)
    weights=2*fc*np.sinc(2*fc*k)*np.sinc(k/2)
    weights=(weights/weights.sum()).astype(np.float32).astype(float)
    # NCL wgt_runave divides by sum only when the rounded weights sum above1.
    weights/=max(float(weights.sum()),1.)
    a=np.asarray(values,dtype=float);out=np.full(a.shape,np.nan)
    windows=np.lib.stride_tricks.sliding_window_view(a,5,axis=-1)
    out[...,2:-2]=np.sum(windows*weights,axis=-1)
    return out


def read_archive(path):
    path=Path(path)
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != ARCHIVE_SHA256:
        raise ValueError('ob_hp0.nc hash mismatch; refusing a substituted reference')
    with netcdf_file(path,'r',mmap=False) as nc:
        arrays={}
        for key in ('height','Fn','Fm','Nd','Md'):
            v=nc.variables[key];a=np.asarray(v[:],dtype=float).copy()
            fill=v._attributes.get('_FillValue',9.96921e36)
            arrays[key]=np.where(a==fill,np.nan,a)
        attrs={k:v.decode() if isinstance(v,bytes) else v for k,v in nc._attributes.items()}
    return dict(arrays=arrays,attributes=attrs,sha256=actual,bytes=path.stat().st_size,
                member='ob_hp0.nc', doi='10.5281/zenodo.10077724',
                archive_url='https://zenodo.org/api/records/10077724/files/McSnow_paper_plot_data.zip/content',
                credit='Welss, Siewert and Seifert (2023), CC-BY-4.0; external member retained outside Git',
                evidence='Archived author numerical calculations, not observations')


def comparison(result,reference):
    z=result['height'];data=reference['arrays']
    if not np.array_equal(z,data['height']):
        raise ValueError('reference height bins differ from declared source output')
    output={}
    for mode,raw in result['profiles'].items():
        filtered=lanczos5(raw)
        comparison_rows={}
        for i,key in enumerate(MOMENTS[:4]):
            obs=data[key];model=filtered[i]
            valid=np.isfinite(obs)&np.isfinite(model)&(obs>0)&(z<SOURCE_TOP)
            domains={}
            for label,mask in [('all_populated',valid),('below_nucleation',valid&(z<SOURCE_BOTTOM))]:
                d=model[mask]-obs[mask]
                domains[label]=dict(count=int(mask.sum()),relative_L1=float(np.sum(np.abs(d))/np.sum(obs[mask])),
                                   relative_bias=float(np.sum(d)/np.sum(obs[mask])),
                                   relative_RMSE=float(np.sqrt(np.mean(d*d))/np.mean(obs[mask])))
            comparison_rows[key]=domains
        output[mode]=dict(raw_binned={k:raw[i].tolist() for i,k in enumerate(MOMENTS)},
                          filtered_binned={k:[float(x) if np.isfinite(x) else None for x in filtered[i]] for i,k in enumerate(MOMENTS)},
                          comparison=comparison_rows)
    return output


def source_identity():
    names = ['welss_column.py', 'boehm.py', 'welss_column_sources.json',
             'references/welss_microphysics_reference.json', 'references/welss_reference_sources.json',
             'references/build_welss_reference.py', 'references/mcsnow_welss_constants.inc',
             'references/mcsnow_welss_functions.inc', 'references/mcsnow_types.inc',
             'references/mcsnow_functions.inc', 'references/mcsnow_bsd_3_clause.txt',
             'references/welss_filter_reference.json', 'references/build_welss_filter_reference.py']
    return {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names}


def product(result,reference,identity):
    if source_identity() != identity:
        raise RuntimeError('benchmark source or fixture changed during integration')
    return dict(schema=SCHEMA,producer=dict(model='atmosphere/snow/welss_column.py', files=identity),
                evidence='Declared public-source deterministic deposition/settling comparison with archived author calculations; not empirical validation or exact archived-run reproduction.',
                source_revision=REVISION,source_constants=SOURCE_CONSTANTS,
                configuration=dict(case='ob_hp0',habits=False,aggregation=False,riming=False,
                    nucleation_rate_m3_s=SOURCE_RATE,source_height_m=[SOURCE_BOTTOM,SOURCE_TOP],
                    seed_gamma_shape=2.,seed_mean_mass_kg=SEED_MEAN,ventilation=True,
                    atmosphere='Pinned 2024 linear interpolator; h1=0, h2=600, ssat=.05, lwc0=.0002',
                    end_time_s=END_TIME,injection='Continuous-rate expectation of the 5 s stochastic source; no finite-step McSnow replay',
                    averaging=dict(last5h=[18000.,36000.],four_snapshots=[34200.,34800.,35400.,36000.],final=[36000.]),
                    output='20 m averages over [z,z+20), matching binned Nv, Nmv, N and Nm; one stated NCL filter pass'),
                numerical_settings=result['settings'],diagnostics=result['diagnostics'],
                dependencies=dict(numpy=np.__version__,scipy=scipy.__version__),
                units=dict(height='m',Fn='m^-2 s^-1',Fm='kg m^-2 s^-1',Nd='m^-3',Md='kg m^-3',deposition='kg m^-3 s^-1'),
                height_m=result['height'].tolist(),reference={k:v for k,v in reference.items() if k!='arrays'},
                profiles=comparison(result,reference))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mass-nodes',type=int,default=16)
    parser.add_argument('--height-nodes',type=int,default=2)
    parser.add_argument('--bin-nodes',type=int,default=3)
    parser.add_argument('--rtol',type=float,default=1e-7)
    parser.add_argument('--max-step-m',type=float,default=20.)
    parser.add_argument('--wall-limit-s',type=float,default=600.)
    args=parser.parse_args()
    reference=read_archive(args.reference)
    identity=source_identity()
    result=integrate(mass_nodes=args.mass_nodes,height_nodes_per_cell=args.height_nodes,
                     bin_nodes=args.bin_nodes,rtol=args.rtol,max_step_m=args.max_step_m,
                     wall_limit_s=args.wall_limit_s,progress=True)
    if hashlib.sha256(args.reference.read_bytes()).hexdigest() != reference['sha256']:
        raise RuntimeError('archive input changed during integration')
    args.output.write_text(json.dumps(product(result,reference,identity),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
