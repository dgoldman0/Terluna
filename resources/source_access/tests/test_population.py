"""Disjoint population, conserved materials and conditional characterization."""
from copy import deepcopy
import json,math
import numpy as np
import pytest
from resources.source_access import population as p


def inputs():return json.loads((p.HERE/'population_inputs.json').read_text())

def header():
 return {'reference_diameter_km':1.,'reference_number':813.,'log10_diameter_knots_km':[-3.,-.58,-.08,.116,1.48],'slopes':[-2.5,-1.5,-1.7,-2.6]}

def row(number=999999,H=18.,diameter=.4,a=1.,e=.1,i=3.):
 return {'spkid':str(20000000+number),'number':number,'full_name':str(number)+' (2099 AA)','designation':'2099 AA','H':str(H) if H is not None else None,'diameter':str(diameter) if diameter is not None else None,'diameter_sigma':'0','albedo':'.2','spec_B':'S','spec_T':None,'a':str(a),'e':str(e),'i':str(i)}


def test_catalogue_diameter_precedence_and_unmeasured_objects_retained():
 cfg=inputs();cfg['draws']=16;cfg['density_bounds_kg_m3']=[[1000,1000]]*3
 c=p.Characterization(cfg,header());records=[row(H=30),row(999998,H=18,diameter=None),row(999997,H=None,diameter=None)]
 known,audit=p.catalogued(records,{},cfg,c,np.random.default_rng(5))
 assert audit['known_measured_diameter_count']==1 and audit['known_inferred_diameter_count']==1
 assert audit['unquantified_no_H_or_diameter_ids']==[records[2]['spkid']]
 expected=math.pi/6*400**3*1000
 assert np.all(known['narrow'][:,0]>expected)
 assert np.array_equal(known['narrow'],known['all_nea'])


def test_science_exclusion_is_not_reintroduced_as_missing_population():
 cfg=inputs();cfg['draws']=16;protected=row(433,H=18,diameter=2)
 other=row(999998,H=18,diameter=.1,a=2,e=.5)
 c=p.Characterization(cfg,header());known,audit=p.catalogued([protected,other],{},cfg,c,np.random.default_rng(5))
 assert protected['spkid'] in audit['protected_ids_excluded']
 bins=[{'H_low':17.75,'H_high':18.25,'count_low':2.,'count_central':2.,'count_high':2.}]
 orbits=[{'H_low':17.75,'H_high':18.25,'category_probabilities':[.5,0,.5]}]
 result=p.residual_counts([protected,other],bins,orbits,cfg['geometry'])[0]
 assert result['known_count_by_category']==[1,0,1]
 assert result['missing_mean_by_branch']==[[0.,0.,0.]]*3
 assert np.all(known['narrow']==0)


def test_metals_and_water_are_fractional_and_jointly_bounded():
 cfg=inputs();c=p.Characterization(cfg,header());n=10000
 v=c.materials(np.ones(n),np.full(n,.05),np.ones(n,dtype=int),np.random.default_rng(9))
 assert np.all(v>=0) and np.all(v[:,1]+v[:,2]<=v[:,0]*(1+1e-14))
 assert np.quantile(v[:,1]/v[:,0],.5)<.10
 assert np.quantile(v[:,2]/v[:,0],.5)<.05
 assert .577e12*math.pi/6<=v[:,0].min() and v[:,0].max()<=3.595e12*math.pi/6


def test_missing_population_branches_and_pools_are_coupled_and_nested():
 cfg=inputs();cfg.update(draws=24,moment_samples=512,explicit_count_threshold=128)
 c=p.Characterization(cfg,header())
 bins=[{'H_low':18.,'H_high':18.5,'missing_mean_by_branch':[[1.,2.,3.],[2.,3.,4.],[3.,4.,5.]]}]
 v,n,_=p.undiscovered(bins,cfg,c,np.random.default_rng(3))
 for j in range(3):
  assert np.all(v[j]['narrow']<=v[j]['wider'])
  assert np.all(v[j]['wider']<=v[j]['all_nea'])
  assert np.array_equal(v[j]['all_nea'],v[j]['all_no_don_quixote'])
  if j:assert np.all(v[j-1]['all_nea']<=v[j]['all_nea'])
 again=p.undiscovered(bins,cfg,c,np.random.default_rng(3))[0]
 assert np.array_equal(v[-1]['all_nea'],again[-1]['all_nea'])


def test_conditional_albedo_prior_is_finite_and_higher_albedo_means_smaller_body():
 cfg=inputs();c=p.Characterization(cfg,header());H=np.full(10000,18.)
 dark=c.albedo(H,np.ones(10000,dtype=int),np.random.default_rng(1))
 bright=c.albedo(H,np.zeros(10000,dtype=int),np.random.default_rng(1))
 assert np.all(np.isfinite(dark)) and np.all(dark>0)
 assert np.median(dark)<np.median(bright)
 assert np.median(1329*10**(-H/5)/np.sqrt(dark))>np.median(1329*10**(-H/5)/np.sqrt(bright))


def test_dense_compound_sum_preserves_material_constraints():
 cfg=inputs();cfg['draws']=128;c=p.Characterization(cfg,header())
 sample=c.synthetic(np.full(10000,23.),np.random.default_rng(7));moments=(sample.mean(0),np.cov(sample,rowvar=False))
 v,n=p.aggregate_population(1000,moments,c,22.75,23.25,cfg,np.random.default_rng(8))
 assert np.all(v>=0) and np.all(v[:,1:].sum(1)<=v[:,0]*(1+1e-14))
 assert abs(n.mean()-1000)<20


def test_missing_count_scenarios_have_no_invented_gaussian_probability():
 cfg=inputs();assert 'not draws' in cfg['population_count_rule']
 assert len(cfg['population_count_bins'])==25
 assert cfg['population_count_bins'][-1]['published_cumulative_central']==9120000


def test_changed_external_bytes_are_rejected(tmp_path):
 f=tmp_path/'data';f.write_bytes(b'new')
 with pytest.raises(ValueError,match='identity mismatch'):p.checked(f,'0'*64)
