"""Linear finite-volume energy-balance climate on an equal-area spherical grid.

Synthetic land/water maps, prescribed albedo and linear OLR; no calibrated lunar
radiative transfer, atmosphere dynamics, moisture, ice feedback or latent heat.
Exact linear propagation for piecewise-constant insolation finds the periodic
orbit directly, rather than accepting an arbitrary spun-up initial temperature.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import eigh

DAY=86400.

@dataclass(frozen=True)
class ClimateConfig:
    nlat: int = 6
    nlon: int = 24
    steps: int = 720
    period_days: float = 29.53059
    solar_w_m2: float = 1361.
    transmission: float = .90  # bolometric input; independent of EUV screening
    albedo: float = .30
    olr_at_273_w_m2: float = 200.
    olr_slope_w_m2_k: float = 2.
    heat_transport_w_m2_k: float = 1.
    atmospheric_column_kg_m2: float = 75000.
    participating_atmosphere: float = .10
    cp_air: float = 1005.
    land_capacity_j_m2_k: float = 2e6
    water_depth_m: float = 10.
    layout: str = 'concentrated'


def laplacian(nlat,nlon):
    """Conservative equal-area finite-volume angular Laplacian on a sphere."""
    if nlat<2 or nlon<4: raise ValueError('Grid too small')
    nx=nlat*nlon; dx=2/nlat; dl=2*np.pi/nlon
    x=-1+(np.arange(nlat)+.5)*dx
    L=np.zeros((nx,nx))
    def edge(i,j,w):
        L[i,i]-=w;L[j,j]-=w;L[i,j]+=w;L[j,i]+=w
    for a in range(nlat):
        for b in range(nlon):
            i=a*nlon+b
            edge(i,a*nlon+(b+1)%nlon,1/((1-x[a]**2)*dl**2))
            if a<nlat-1:
                face=-1+(a+1)*dx
                edge(i,(a+1)*nlon+b,(1-face**2)/dx**2)
    return x,L


def solve_climate(cfg=ClimateConfig()):
    if cfg.steps<16 or cfg.period_days<=0 or not 0<=cfg.albedo<1 or not 0<cfg.transmission<=1:
        raise ValueError('Invalid forcing')
    if min(cfg.olr_slope_w_m2_k,cfg.atmospheric_column_kg_m2,cfg.cp_air,cfg.land_capacity_j_m2_k)<=0 or not 0<=cfg.participating_atmosphere<=1 or min(cfg.heat_transport_w_m2_k,cfg.water_depth_m)<0:
        raise ValueError('Invalid thermal inputs')
    x,L=laplacian(cfg.nlat,cfg.nlon);lon=2*np.pi*np.arange(cfg.nlon)/cfg.nlon
    water=np.zeros((cfg.nlat,cfg.nlon),bool)
    if cfg.layout=='concentrated': water[:,:cfg.nlon//2]=True
    elif cfg.layout=='distributed': water[:,::2]=True
    elif cfg.layout=='dry': pass
    else: raise ValueError('Unknown synthetic geography')
    water=water.ravel(); n=len(water)
    capacity=np.full(n,cfg.atmospheric_column_kg_m2*cfg.cp_air*cfg.participating_atmosphere+cfg.land_capacity_j_m2_k)
    capacity+=water*cfg.water_depth_m*1000*4180
    root_c=np.sqrt(capacity)
    A=cfg.heat_transport_w_m2_k*L-cfg.olr_slope_w_m2_k*np.eye(n)
    lam,V=eigh(A/root_c[:,None]/root_c[None,:])
    if max(lam)>=0: raise RuntimeError('Expected dissipative operator')
    dt=cfg.period_days*DAY/cfg.steps
    alpha=np.exp(lam*dt); beta=np.expm1(lam*dt)/lam
    phase=2*np.pi*(np.arange(cfg.steps)+.5)/cfg.steps
    cosine=np.sqrt(1-x[:,None]**2)[None,:,:]*np.maximum(0,np.cos(lon[None,None,:]-phase[:,None,None]))
    # Analytically global mean is 1/4. Normalize finite-grid/time quadrature only.
    correction=.25/float(cosine.mean())
    incident=cfg.solar_w_m2*cfg.transmission*cosine.reshape(cfg.steps,n)*correction
    q=incident*(1-cfg.albedo)-cfg.olr_at_273_w_m2
    f=(q/root_c)@V
    z=np.zeros(n)
    for fk in f:z=alpha*z+beta*fk
    initial=z/(-np.expm1(lam*cfg.period_days*DAY))
    z=initial.copy(); temps=[]; max_residual=0.
    for fk,qk in zip(f,q):
        theta=V@z/root_c
        integ_z=beta*z+(beta-dt)/lam*fk
        integral_theta=V@integ_z/root_c
        znext=alpha*z+beta*fk
        nexttheta=V@znext/root_c
        storage=float(np.mean(capacity*(nexttheta-theta)))
        balance=float(dt*qk.mean()-cfg.olr_slope_w_m2_k*integral_theta.mean())
        max_residual=max(max_residual,abs(storage-balance)/dt)
        # Report interval means, appropriate forcing for subsequent budget screens.
        temps.append(integral_theta/dt+273.15)
        z=znext
    temperature=np.asarray(temps)
    result=dict(layout=cfg.layout,water_fraction=float(water.mean()),
        participating_atmosphere=cfg.participating_atmosphere,
        olr_at_273_w_m2=cfg.olr_at_273_w_m2,olr_slope_w_m2_k=cfg.olr_slope_w_m2_k,
        heat_transport_w_m2_k=cfg.heat_transport_w_m2_k,water_depth_m=cfg.water_depth_m,
        global_mean_K=float(temperature.mean()),minimum_K=float(temperature.min()),
        maximum_K=float(temperature.max()),
        grid_time_fraction_273_313=float(np.mean((temperature>=273.15)&(temperature<=313.15))),
        freezing_flag=bool(temperature.min()<273.15),
        radiative_extrapolation_flag=bool(temperature.min()<200 or temperature.max()>340),
        max_discrete_energy_residual_W_m2=max_residual,
        periodic_residual_K=float(np.max(abs(V@(z-initial)/root_c))),
        quadrature_normalization=correction,
        status='UNCALIBRATED_LINEAR_THERMAL_SCREEN_NOT_HABITABILITY_PROOF')
    for mask,label in ((water,'water'),(~water,'land')):
        if mask.any():
            result[label+'_min_K']=float(temperature[:,mask].min())
            result[label+'_max_K']=float(temperature[:,mask].max())
    return result,dict(temperature_K=temperature,incident_W_m2=incident,
                      sin_latitude=x,longitude_rad=lon,water_mask=water,capacity=capacity)
