import math
import numpy as np
import pytest
from atmosphere.snow.spatial_crystal_growth import MaskDomain, GrowingCrystal, Combined, grid
from atmosphere.snow.spatial_diffusion import Domain, Attachment
from atmosphere.snow.facet_growth import ICE_DENSITY

@pytest.mark.parametrize('rim',[1.,.7,.2])
def test_vectorized_stationary_geometry_matches_checked_solver(rim):
 d=Domain(2.,resolution=8,extent=20.,cavity_radius=1-rim,cavity_depth=1. if rim<1 else 0.)
 m=MaskDomain(d.rf,d.zf,d.solid)
 difference=(d.matrix-m.matrix).data
 assert not len(difference) or np.max(abs(difference))<1e-9
 laws=[Attachment('nucleation',.2,.002),Attachment('nucleation',1.,.007)]
 reference=d.solve(.05,100,laws,heat_parameter=.5,tolerance=1e-11)
 answer=m.coupled(.05,100,laws,.5,tolerance=1e-11)
 assert answer['half_flux']==pytest.approx(reference['half_flux'],rel=1e-8)
 assert answer['relative_mass_residual']<1e-8


def test_surface_diffusion_needs_an_accessible_kink():
 rf=np.arange(16.);zf=np.arange(16.);rr,zz=np.meshgrid((rf[:-1]+rf[1:])/2,(zf[:-1]+zf[1:])/2,indexing='ij')
 solid=((rr<5)&(zz<5))|((rr<2)&(zz<7))
 d=MaskDomain(rf,zf,solid);n=10.;e=d.terrace_enhancement(n)
 top=(d.ice[:,3]==1)&(d.ice[:,5]==7)
 assert np.all(e[top]==0)
 lower=(d.ice[:,3]==1)&(d.ice[:,5]==5)
 coords=d.ice[lower,4];enhance=e[lower]
 np.testing.assert_allclose(enhance,np.exp(-(coords-min(coords))/n),rtol=1e-13)


def test_curvature_flux_has_a_growth_threshold():
 law=Combined(Attachment('constant',.3),np.array([0.,.5,1.]),np.array([.01,.02,.03]))
 p,dp=law.flux(np.array([.005,.02,.1]))
 np.testing.assert_allclose(p,[0.,0.,.07]);np.testing.assert_allclose(dp,[0.,0.,1.])


@pytest.mark.parametrize('dx',[2e-6,1e-6,.5e-6])
def test_seed_volume_exactly_matches_given_mass(dx):
 a=30e-6;c=40e-6;mass=.35*2*math.pi*a*a*c*ICE_DENSITY
 g=GrowingCrystal(dx,a,c,200,500,extent=80.)
 g.measured_seed(a,c,mass,2e-6)
 assert g.mass()==pytest.approx(mass,rel=1e-12)
 assert np.all((g.fraction>=0)&(g.fraction<=1))
 # Force expansion with a partial boundary cell inside the old uniform core.
 g.fraction[g.radial_core-3,0]=.3
 before=g.mass();g.expand()
 assert g.mass()==pytest.approx(before,rel=1e-12)


@pytest.mark.parametrize('heating',[False,True])
def test_growing_cylinder_converges_to_independent_analytic_time(heating):
 # Reflecting top and bottom make a cylinder infinite in z. The radial
 # diffusion/kinetic resistance can be integrated exactly as R grows.
 air=dict(ice_saturation_pressure_Pa=100.,vapor_gas_constant_J_kg_K=460.,vapor_molecular_speed_m_s=500.,
          latent_sublimation_J_kg=2.8e6,diffusivity_m2_s=2e-5,conductivity_W_m_K=.02)
 T=240.;si=.1;alpha=.2;R0=5e-6;R1=10e-6;Rout=100e-6
 rhoeq=air['ice_saturation_pressure_Pa']/(460*T);D=air['diffusivity_m2_s'];X0=D/125
 F=lambda r:.5*r*r*math.log(Rout/r)+.25*r*r
 L=air['latent_sublimation_J_kg'];beta=(L/(460*T)-1)/T
 chi=beta*L*rhoeq*D/air['conductivity_W_m_K'] if heating else 0.
 exact=ICE_DENSITY/(rhoeq*D*si)*((1+chi)*(F(R1)-F(R0))+X0/alpha*(R1-R0))
 errors=[]
 for dx in (1.25e-6,.625e-6,.3125e-6):
  g=GrowingCrystal(dx,R0,dx,64,1,extent=20.,adaptive=False)
  g.rf=grid(64,Rout/dx);g.zf=np.array([0.,1.]);g.far=Rout/dx
  g.volume=math.pi*np.diff(g.rf**2)[:,None];g.fraction=((g.rf[:-1]+g.rf[1:])[:,None]/2*dx<R0).astype(float)
  g.solid=g.fraction>0;g.reflect_top=True;g.initial_mass=g.mass()
  while g.state()['a_m']<R1-dx/10:
   g.advance(air,T,si,[Attachment('constant',alpha)]*2,1e10,kinks=False,heating=heating,tolerance=1e-11)
  errors.append(abs(g.time/exact-1))
  assert g.state()['mass_ledger_relative_residual']<1e-10
 assert errors[2]<.03
 assert errors[2]<errors[1]<errors[0]
 print('cylinder temporal errors',errors)


def test_tilted_surface_kinetic_uptake_converges_to_cone_area():
 # Independent continuum surface: a right cone with radius=height=1 has
 # lateral area sqrt(2)*pi. Unweighted staircase faces overcount it by sqrt(2).
 errors=[]
 for n in (20,40,80):
  edges=np.arange(3*n+1)/n;centers=(edges[:-1]+edges[1:])/2
  r,z=np.meshgrid(centers,centers,indexing='ij')
  d=MaskDomain(edges,edges,r+z<1-1e-12)
  neighbors=np.bincount(d.ice[:,0].astype(int),minlength=d.cells)
  weights=1/np.sqrt(neighbors[d.ice[:,0].astype(int)])
  laws=[Combined(Attachment('constant',.3),0.,kinetic_weight=weights[d.ice[:,3]==k]) for k in (0,1)]
  answer=d.coupled(.1,1e-6,laws,tolerance=1e-11)
  inferred_area=answer['half_flux']/(.3*.1*1e-6)
  errors.append(abs(inferred_area/(math.sqrt(2)*math.pi)-1))
 assert errors[-1]<.04
 assert errors[2]<errors[1]<errors[0]


def test_sealed_air_is_excluded_from_ambient_supply_without_becoming_ice():
 edges=np.arange(11.)
 solid=np.zeros((10,10),bool);solid[:5,:5]=True;solid[2,2]=False
 closed=MaskDomain(edges,edges,solid)
 filled=solid.copy();filled[2,2]=True
 reference=MaskDomain(edges,edges,filled)
 assert closed.closed_air.sum()==1
 assert not closed.solid[2,2]
 assert closed.index[2,2]==-1
 assert (closed.matrix-reference.matrix).count_nonzero()==0
 laws=[Attachment('nucleation',1.,.006),Attachment('nucleation',1.,.002)]
 got=closed.coupled(.05,100.,laws,.5)
 expected=reference.coupled(.05,100.,laws,.5)
 assert got['half_flux']==pytest.approx(expected['half_flux'],rel=1e-13)
 assert got['relative_mass_residual']<1e-8


def test_one_cell_open_channel_retains_its_gas_supply():
 edges=np.arange(11.)
 solid=np.zeros((10,10),bool);solid[:5,:5]=True;solid[2,2:6]=False
 d=MaskDomain(edges,edges,solid)
 assert not np.any(d.closed_air)
 assert d.exterior_air[2,2]
 assert d.index[2,2]>=0
 result=d.coupled(.05,100.,[Attachment('constant',.2)]*2)
 assert result['half_flux']>0
 assert result['relative_mass_residual']<1e-8


def test_reflecting_top_is_not_an_external_humidity_reservoir():
 from atmosphere.snow.spatial_crystal_growth import exterior_air_masks
 solid=np.zeros((8,8),bool);solid[2:6,2:8]=True;solid[3,3:8]=False
 open_air,open_closed=exterior_air_masks(solid,reflect_top=False)
 reflected_air,reflected_closed=exterior_air_masks(solid,reflect_top=True)
 assert open_air[3,3] and not open_closed[3,3]
 assert reflected_closed[3,3] and not reflected_air[3,3]


def test_growth_freezes_sealed_partial_ice_and_bounds_its_finite_vapour():
 air=dict(ice_saturation_pressure_Pa=400.,vapor_gas_constant_J_kg_K=460.,
          vapor_molecular_speed_m_s=560.,latent_sublimation_J_kg=2.8e6,
          diffusivity_m2_s=2e-5,conductivity_W_m_K=.024)
 g=GrowingCrystal(1e-6,5e-6,5e-6,10,10,extent=4.,adaptive=False)
 g.solid[2,2]=False;g.fraction[2,2]=.3;g.initial_mass=g.mass()
 initial=g.mass()
 result=g.advance(air,268.,.05,[Attachment('constant',.2)]*2,1.,kinks=False)
 assert g.fraction[2,2]==.3
 assert not g.solid[2,2]
 assert g.mass()>initial
 assert g.state()['mass_ledger_relative_residual']<1e-10
 volume=2*math.pi*(3**2-2**2)*g.dx**3*.7
 assert result['closed_gas_volume_m3']==pytest.approx(volume)
 assert result['trapped_vapor_upper_bound_kg']==pytest.approx(400/(460*268)*1.05*volume)
 assert result['closed_air_cells']==1
