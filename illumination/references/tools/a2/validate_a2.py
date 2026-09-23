"""Assemble and validate the A2 atmospheric/cloud accuracy checkpoint."""
from __future__ import annotations
import hashlib, importlib.util, json, math, subprocess, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
GEN=ROOT/'data/a2/generated'
# The A2 lab page and its checker live with the other labs in visualization/labs.
LABS=ROOT.parents[1]/'visualization'/'labs'
LAB_FILES={'tools/a2/check_lab.py','src/a2-accuracy-lab.js','index.a2-accuracy-lab.html','build_a2_accuracy.py'}
def locate(name):return (LABS if name in LAB_FILES else ROOT)/name
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def main():
 start=time.monotonic();checks=[]
 def gate(name,passed,detail):checks.append({'name':name,'passed':bool(passed),'detail':detail});return passed
 # Existing JS regression suite, including A1 source-level tests.
 proc=subprocess.run(['node','--test',*map(str,sorted((ROOT/'tests').glob('*.test.cjs')))],cwd=ROOT,text=True,capture_output=True)
 test_lines=proc.stdout+proc.stderr
 import re
 def grab(label):
  m=re.search(r'^# '+re.escape(label)+r'\s+(\d+)\s*$',test_lines,re.M);return int(m.group(1)) if m else None
 tests,passed,failed=grab('tests'),grab('pass'),grab('fail');gate('JavaScript regression suite',proc.returncode==0 and failed==0,{'tests':tests,'pass':passed,'fail':failed})
 # Molecular full-spectrum references.
 standard=load(GEN/'moon-fair-molecular-standard-full.json');fine=load(GEN/'moon-fair-molecular-fine-full-chunked.json');atlas=load(GEN/'moon-fair.multiple-sky.json');atlas_np=np.load(GEN/'moon-fair.multiple-sky.npz');a1_np=np.load(ROOT/'data/a1/generated/moon-fair.spectra.npz')
 gate('Fine chunks cover exact 360–830 nm / 10 nm grid',fine['wavelength_nm']==list(map(float,range(360,831,10))),{'bands':len(fine['wavelength_nm'])})
 max_chunk_increment=max(c['final_increment_fraction'] for c in fine['chunks']);gate('Fine spectral chunks converge by order increment',max_chunk_increment<5e-4,{'max_final_increment_fraction':max_chunk_increment,'orders':[c['orders_computed'] for c in fine['chunks']]})
 conv=fine['convergence'];gate('Standard/fine molecular spectral convergence',conv['overall_standard_to_fine_spectral_L1']<.01 and conv['max_ray_standard_to_fine_luminance_Y_relative']<.02,conv)
 # Full-atlas first-order cross-check against A1.
 single=atlas_np['single'];ref=a1_np['radiance'];delta=np.abs(single-ref);den=np.abs(ref).sum(axis=-1);global_l1=float(delta.sum()/max(1e-30,np.abs(ref).sum()));bright=den>0.01*den.max();per=np.divide(delta.sum(axis=-1),den,out=np.zeros_like(den),where=den>0);bright_stats={'threshold_fraction_of_max_spectral_sum':.01,'pixels':int(bright.sum()),'median':float(np.median(per[bright])),'p95':float(np.quantile(per[bright],.95)),'max':float(per[bright].max())}
 gate('Full-atlas first-order reproduction',global_l1<.005 and bright_stats['p95']<.02,{'global_spectral_L1':global_l1,'bright_pixel_relative_spectral_L1':bright_stats})
 # Added scattering must remain nonnegative.
 multiple=atlas_np['multiple'];min_added=float((multiple-single).min());gate('Higher scattering orders add nonnegative radiance',min_added>=-1e-12,{'minimum_multiple_minus_single':min_added})
 # Surface downward energy diagnostic for sun well above horizon.
 sp=importlib.util.spec_from_file_location('a2mol',ROOT/'tools/a2/molecular_multiple_scattering.py');M=importlib.util.module_from_spec(sp);sp.loader.exec_module(M)
 state,rows,oz,shells,profile_sha=M.a1.unpack_profile(ROOT/'data/a1/generated/moon-fair.profile.json');spec=load(ROOT/'data/a1/spectral-inputs.json');lam=atlas_np['wavelength_nm'];beta=1.24062e-6*(lam/1000.)**-4;sigma=np.array(spec['ozone_m2'],float);solar=np.array(spec['solar_W_m2_nm'],float);nodes,w=np.polynomial.legendre.leggauss(24);mu=np.sin(np.deg2rad(atlas_np['elevations_deg']));phi=np.deg2rad(atlas_np['azimuths_deg']);energy=[]
 for si,sun in enumerate(atlas_np['suns_deg']):
  if sun<15:continue
  diffuse=np.empty(len(lam))
  for li in range(len(lam)):
   azint=np.trapezoid(multiple[si,:,:,li],phi,axis=1)*2;diffuse[li]=np.trapezoid(azint*mu,mu)
  smu=math.sin(math.radians(float(sun)));ca,co,blocked=M.a1.columns(state['planet']['radius'],state['upper']['top_m'],2.,smu,1e99,rows,oz,shells,nodes,w);direct=solar*np.exp(-beta*ca-sigma*M.DU*co)*smu;ratio=(direct+diffuse)/(solar*smu);energy.append({'sun_deg':float(sun),'min_ratio':float(ratio.min()),'median_ratio':float(np.median(ratio)),'max_ratio':float(ratio.max())})
 gate('High-Sun downward energy remains below incident horizontal spectral flux',max(e['max_ratio'] for e in energy)<=1.01,energy)
 # Cloud capture identity and graphics integrity.
 cap_ids=['grid-coarse','grid-a1','grid-fine','crop-075','crop-150','crop-200'];caps={i:load(GEN/f'capture-{i}.json') for i in cap_ids};profile_expected=next(x for x in load(ROOT/'data/a1/generated/profiles.json')['profiles'] if x['id']=='moon-fair')['profile_sha256']
 cap_ok=True
 for i,c in caps.items():
  e=c['entry'];field=GEN/e['field_file'];cap_ok &= c['profile_sha256']==profile_expected and sha(field)==e['field_sha256'] and not e['gpu']['context_lost'] and e['field_parity']['maximum_absolute_extinction_error_m1']<1e-7
 gate('All cloud capture identities and GPU/CPU sampler checks pass',cap_ok,{i:{'field_sha256':c['entry']['field_sha256'],'sampler_max_abs_m1':c['entry']['field_parity']['maximum_absolute_extinction_error_m1'],'context_lost':c['entry']['gpu']['context_lost']} for i,c in caps.items()})
 # Continuous-field discretization trend from common deterministic probes.
 vox={i:caps[i]['entry']['voxelization']['sum_absolute_difference_over_sum_continuous'] for i in ['grid-coarse','grid-a1','grid-fine']};gate('Cloud voxelization improves monotonically with resolution',vox['grid-coarse']>vox['grid-a1']>vox['grid-fine'],vox)
 cloud=load(GEN/'cloud-convergence.json');g=cloud['grid_convergence'];gate('A1 cloud grid approaches fine first-order transport',g['grid-a1']['radiance_L1_vs_grid_fine']<.10 and g['grid-a1']['transmittance_L1_vs_grid_fine']<.01,g)
 c=cloud['crop_convergence'];gate('A1 crop is converged for central first-order transport',c['grid-a1']['radiance_L1_vs_crop_200']<1e-4 and c['grid-a1']['transmittance_L1_vs_crop_200']<1e-4,c)
 gate('Cloud Monte Carlo paths resolve without safety-cap failures',cloud['mc']['unresolved_total']==0,{'unresolved_total':cloud['mc']['unresolved_total'],'samples_per_ray':cloud['mc']['samples_per_ray']})
 # MC mean differences are reported as sensitivity, with standard errors retained rather than hard-gated as a deterministic error bound.
 # Source identity and production boundary.
 source_files=['tools/a2/molecular_multiple_scattering.py','tools/a2/build_multiple_atlas.py','tools/a2/merge_molecular_chunks.py','tools/a2/capture_cloud_convergence.py','tools/a2/cloud_convergence.py','tools/a2/validate_a2.py','tools/a2/check_lab.py','tests/a2_frozen_probe.js','src/a2-accuracy-lab.js','index.a2-accuracy-lab.html','build_a2_accuracy.py','A2_METHODS.md']
 result={'schema':'open-moon-a2-validation/1','generated_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'status':'reference-checkpoint-passed' if all(x['passed'] for x in checks) else 'failed','scope':'A2 adds shared-profile molecular multiple scattering, a full standard-grid sky atlas, fine-grid spectral/ray validation, and separated cloud voxel/crop convergence. Production shoreline rendering remains unchanged.','checks':checks,'unit_tests':{'tests':tests,'pass':passed,'fail':failed},'molecular':{'profile_sha256':profile_sha,'standard_sha256':sha(GEN/'moon-fair-molecular-standard-full.json'),'standard_field_sha256':sha(GEN/'moon-fair-molecular-standard-full.npz'),'fine_chunked_sha256':sha(GEN/'moon-fair-molecular-fine-full-chunked.json'),'atlas_sha256':sha(GEN/'moon-fair.multiple-sky.json'),'atlas_spectra_sha256':sha(GEN/'moon-fair.multiple-sky.npz'),'global_first_order_L1_vs_A1':global_l1,'bright_pixel_stats':bright_stats,'fine_convergence':conv,'fine_chunk_final_increment_max':max_chunk_increment,'high_sun_energy':energy,'selected_ray_multiple_to_single_luminance':[r['multiple_to_single_luminance_Y'] for r in fine['rays']]},'cloud':{'capture_sha256':{i:sha(GEN/f'capture-{i}.json') for i in cap_ids},'field_voxelization_relative_L1':vox,'transport_convergence':cloud['grid_convergence'],'crop_convergence':cloud['crop_convergence'],'mc_grid_mean_L1_vs_fine':cloud['mc']['grid_mean_L1_vs_fine'],'mc_crop_mean_L1_vs_crop_200':cloud['mc']['crop_mean_L1_vs_crop_200'],'mc_samples_per_ray':cloud['mc']['samples_per_ray']},'source_sha256':{p:sha(locate(p)) for p in source_files},'limitations':['Molecular atlas uses the standard angular/radial grid; fine-grid convergence is sampled on seven rays across the full spectrum.','The atmospheric profile, ozone and upper extension remain specified inputs rather than climate-equilibrium predictions.','Cloud-grid comparison remains a frozen fair-weather field; continuous-field error is sampled at deterministic probes rather than bounded everywhere.','Cloud multiple-scattering crop/grid sensitivity is Monte Carlo and retains sampling uncertainty.','Production atmospheric/cloud lighting, cache interpolation and display calibration remain unchanged and unvalidated against A2.'],'elapsed_s':time.monotonic()-start}
 (ROOT/'A2_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'checks':[(x['name'],x['passed']) for x in checks],'elapsed_s':result['elapsed_s']},indent=2));return 0 if result['status'].endswith('passed') else 1
if __name__=='__main__':raise SystemExit(main())
