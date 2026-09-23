"""Audit A3 implementation, comparisons and still-open accuracy gates.

Failed inherited-A2 agreement checks remain failures. They are not relabeled as
passing implementation tests. This report distinguishes a usable benchmark
from qualification of the saved A2 field or the production renderer.
"""
from __future__ import annotations
import hashlib,json,math,re,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import coupled_transport as t

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 benchpath=t.ROOT/'data/a3/generated/benchmark.json';B=read(benchpath);tests=read(t.ROOT/'validation/a3/transport-tests.json');extra=read(t.ROOT/'validation/a3/extension-tests.json');diagnosis=read(t.ROOT/'validation/a3/clear-diagnosis.json');medium=read(t.ROOT/'validation/a3/clear-refine-medium.json')
 densepath=t.ROOT/'validation/a3/clear-refine-dense.json';dense=read(densepath) if densepath.exists() else None
 impl=[r for r in tests['checks'] if not r['name'].startswith('zero-cloud full MC reproduces A2')]+extra['checks'];inherited=[r for r in tests['checks'] if r['name'].startswith('zero-cloud full MC reproduces A2')]
 tap=(t.ROOT/'validation/a3/unit-tests.tap').read_text();p=re.search(r'^# pass (\d+)$',tap,re.M);f=re.search(r'^# fail (\d+)$',tap,re.M);node=dict(pass_count=int(p[1]) if p else 0,fail_count=int(f[1]) if f else None)
 gates=[]
 def add(name,ok,**data):gates.append(dict(name=name,state='PASS' if ok else 'FAIL',passed=bool(ok),**data))
 add('Analytic, phase, geometry, limiting-case and mixed-transmission checks',all(c['passed'] for c in impl),passed_checks=sum(c['passed'] for c in impl),total_checks=len(impl))
 add('Inherited JavaScript regression suite',node['pass_count']==159 and node['fail_count']==0,**node)
 add('All 48 wavelengths from 360 to 830 nm',np.array_equal(B['wavelength_nm'],np.arange(360,831,10)))
 add('Full benchmark has no unresolved safety-cap or majorant failures',B['unresolved_paths']==0,unresolved=B['unresolved_paths'])
 current=t.load_inputs();matched=all(B['identity'][key]==current[field] for key,field in [('profile','profile_sha256'),('field','field_sha256'),('molecular','molecular_field_sha256')]);add('Profile/cloud/A2 identities remain matched',matched)
 source_match=all(sha(t.ROOT/f)==s for f,s in B['identity']['sources'].items());add('Executed transport kernel and spectral source hashes match delivered files',source_match)
 zmax=0.;closure=0.;seed_tests=0;statistics=[]
 for case in B['cases']:
  for name,s in case['treatments'].items():
   for b in s['batches_by_wavelength']:
    v0,v1=b;den=math.hypot(v0['standard_error'][0],v1['standard_error'][0]);z=abs(v0['mean'][0]-v1['mean'][0])/den if den>0 else 0.;zmax=max(zmax,z);seed_tests+=1
   if name in ['coupled_mc','clear_mc','clear_source_mc']:
    m=np.array(s['mean']);closure=max(closure,float(np.max(abs(m[:,0]-m[:,1:].sum(axis=1))/np.maximum(1e-15,m[:,0]))))
  full=case['treatments']['coupled_mc']['color'];fixed=case['treatments']['clear_source_mc']['color'];a2fixed=case['treatments']['a2_frozen_source']['color'];clear=case['treatments']['clear_mc']['color'];rep=case['deterministic_color']['matched_production_expression'];oldclear=case['deterministic_color']['a2_clear'];se=math.hypot(full['XYZ_standard_error'][1],fixed['XYZ_standard_error'][1]);delta=full['XYZ'][1]-fixed['XYZ'][1]
  statistics.append(dict(case=case['id'],coupled_Y=full['XYZ'][1],coupled_Y_standard_error=full['XYZ_standard_error'][1],fixed_clear_source_Y=fixed['XYZ'][1],feedback_delta_Y=delta,feedback_delta_Y_standard_error=se,feedback_fraction_of_fixed=delta/fixed['XYZ'][1],feedback_z=delta/se,mixed_scattering_fraction=full['components_Y'][2]/full['XYZ'][1],matched_production_expression_fractional_difference=rep['XYZ'][1]/full['XYZ'][1]-1,matched_expression_vs_same_a2_source_fractional_difference=rep['XYZ'][1]/a2fixed['XYZ'][1]-1,a2_fixed_source_fractional_difference=a2fixed['XYZ'][1]/fixed['XYZ'][1]-1,a2_clear_fractional_difference=oldclear['XYZ'][1]/clear['XYZ'][1]-1,cloud_view_optical_depth=case['cloud_optical_depth']))
 add('Independent seed batches are statistically consistent at the declared 6-SE check',zmax<=6.,maximum_batch_z=zmax,comparisons=seed_tests)
 add('Gas-only, cloud-only and mixed path contributions close to total radiance',closure<1e-11,largest_relative_residual=closure)
 sq=max(r['single_quadrature_spectral_L1'] for r in B['cases']);rq=max(r['replay_quadrature_spectral_L1'] for r in B['cases'])
 add('Independent coupled first-order integration converges under Gauss refinement',sq<.005,maximum_spectral_L1=sq)
 add('Matched production-expression replay converges under Gauss refinement',rq<.002,maximum_spectral_L1=rq)
 implementation=all(g['passed'] for g in gates)
 gates.append(dict(name='Saved A2 multiple-scattering field agrees with independent clear-gas MC',state='PASS' if all(r['passed'] for r in inherited) else 'FAIL — selected blue wavelengths',passed=all(r['passed'] for r in inherited),kind='inherited_reference_qualification',checks=inherited))
 audit=[]
 for wave in [400,460,580]:
  rr=[r for r in diagnosis['results'] if r['wavelength_nm']==wave and r['mode']==1]
  # The two estimates use independent seeds and two roulette survival settings.
  # Equal-N pooled average; SE of their independent average.
  mean=sum(r['mean'] for r in rr)/len(rr);se=math.sqrt(sum(r['SE']**2 for r in rr))/len(rr);refined=None
  reference=dense or medium
  if wave in reference['lam']:refined=reference['zenith'][reference['lam'].index(wave)]
  audit.append(dict(wavelength_nm=wave,saved_a2=rr[0]['A2'],mc_mean=mean,mc_se=se,mc_samples=sum(r['samples'] for r in rr),refined=refined,relative_mc_minus_saved=(mean-rr[0]['A2'])/rr[0]['A2'],refined_fractional_difference=None if refined is None else (refined-mean)/mean))
 gates.append(dict(name='Independent dense molecular refinement',state='COMPLETED — two diagnostic wavelengths' if dense else 'PARTIAL — medium refinement complete; dense run pending',kind='reference_diagnostic',data=audit))
 for text in ['Coupled transport convergence beyond the frozen grid and crop','Production GPU closure, cache, reflections and complete image comparison','Reflecting-ground and finite-Sun boundary experiments']:
  gates.append(dict(name=text,state='OPEN',kind='future_accuracy_gate'))
 source_files=['tools/a3/'+n for n in ['coupled_transport.py','deterministic.py','run_benchmark.py','test_transport.py','test_extensions.py','diagnose_clear.py','refine_molecular.py','validate_a3.py']]
 out=dict(schema='open-moon-a3-validation/1',generated_utc=datetime.now(timezone.utc).isoformat(),status='coupled-benchmark-implemented; inherited-A2-blue-agreement-failed',implementation_passed=implementation,all_accuracy_gates_passed=False,production_replacement_qualified=False,benchmark_sha256=sha(benchpath),source_sha256={f:sha(t.ROOT/f) for f in source_files},gates=gates,implementation_checks=impl,case_statistics=statistics,clear_audit=audit,medium_refinement=medium,dense_refinement=dense,benchmark_samples=B['samples_per_wavelength_per_treatment'],main_wavelength_histories=len(B['cases'])*4*48*B['samples_per_wavelength_per_treatment'],browser_validation='Recorded separately in validation/a3/lab-browser.json after exact-build rendering.',limitations=['Monte Carlo standard errors quantify sampling, not physical-model error.','The finite cloud grid/crop remains a specified benchmark input.','A2 blue-band all-order disagreement is retained; no gain or color correction is applied.','Matched production-expression comparison uses spectral shared-profile inputs and black ground. It is not the unchanged production renderer output.'])
 (t.ROOT/'A3_VALIDATION.json').write_text(json.dumps(out,indent=2)+'\n')
 # Small actual-field cross-sections for the diagnostic viewer.
 c=t.load_inputs();nx,ny,nz=map(int,c['dims']);ext=c['ext'].reshape(nz,ny,nx);slices={}
 for r in B['cases']:
  z=float(r['origin_m'][2]);iz=int(np.clip(round((z-c['lo'][2])/(c['hi'][2]-c['lo'][2])*(nz-1)),0,nz-1));v=ext[iz]
  slices[r['id']]=dict(width=nx,height=ny,z_index=iz,maximum=float(c['ext'].max()),values=v.ravel().tolist())
 display=dict(slices=slices,audit_rows=audit,audit_note='The independent zero-cloud check exposes a higher-order blue-band convergence gap in the saved A2 field. Full coupled MC and its fixed-clear-source MC control avoid that grid; A2 remains an explicitly measured comparator.',audit_foot='Spectral radiance in W m⁻² sr⁻¹ nm⁻¹, at the +45° Sun / zenith diagnostic. MC uses two independent 1,048,576-history runs with different roulette settings. '+('The dense ' if dense else 'The medium ')+'deterministic refinement changes numerical grids and order convergence, not the physical atmosphere.')
 (t.ROOT/'data/a3/generated/display.json').write_text(json.dumps(display,separators=(',',':'))+'\n')
 print(json.dumps(dict(implementation_passed=implementation,gates=[dict(name=g['name'],state=g['state']) for g in gates],case_statistics=statistics),indent=2))
 return implementation
if __name__=='__main__':raise SystemExit(0 if main() else 1)
