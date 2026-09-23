"""Check identities and numeric evidence, then write the compact A1 validation record.

This command validates the reference checkpoint. It explicitly leaves production
integration, field discretization convergence and global physical accuracy open.
"""
from __future__ import annotations
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,re,subprocess,importlib.metadata
ROOT=Path(__file__).resolve().parents[2]
GEN=ROOT/'data/a1/generated'
REPO=ROOT.parents[1]
LABS=REPO/'visualization'/'labs'
COLUMN=REPO/'atmosphere'/'column'
VIEWER=REPO/'immersion'/'light-cycle'/'shoreline'
BASELINE_HTML='c7ead1fd83e0d55518d2951e08cd46d31f323280885125382a8013b5d74636d9'
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p:Path):return json.loads(p.read_text())
def main():
    manifest=read(GEN/'profiles.json');sky_checks=[]
    for path,digest in manifest['sources'].items():assert sha(REPO/path)==digest, 'Stale profile source: '+path
    for p in manifest['profiles']:
        actual=sha(GEN/p['file']);assert actual==p['profile_sha256'];s=read(GEN/(p['id']+'.sky.json'));assert s['profile_sha256']==actual
        assert s['solver_sha256']==sha(ROOT/'tools/a1/reference_optics.py');assert s['spectral_sha256']==sha(ROOT/'data/a1/spectral-inputs.json')
        assert all(x['converged'] for x in p['independent_js_paths'])
        col=max(max(x['relative_to_scipy']) for x in s['column_validation'] if x['relative_to_scipy'])
        js=max(abs(v-x['scipy'][i])/max(1,x['scipy'][i]) for x in s['column_validation'] if x['scipy'] for i,v in enumerate(x['js']))
        ray=max(x['relativeSpectralL1'] for x in s['ray_convergence'])
        assert col<1e-6 and js<1e-6 and ray<.01,(p['id'],col,js,ray)
        sky_checks.append({'id':p['id'],'profile_sha256':actual,'gauss_vs_scipy_max_scaled_relative_column_error':col,'js_vs_scipy_max_scaled_relative_column_error':js,'six_ray_max_relative_spectral_L1':ray,'numerical_gate_passed':True})
    gpu=read(GEN/'moon-fair.gpu.json');mc=read(GEN/'moon-fair.mc.json');field=read(GEN/'moon-fair.cloud.json');field_sha=sha(GEN/'moon-fair.cloud.json')
    assert gpu['passed'] and mc['passed'] and not gpu['errors'] and not gpu['external_requests'];assert gpu['field_sha256']==mc['field_sha256']==field_sha
    assert mc['solver_sha256']==sha(ROOT/'tools/a1/cloud_reference.py');assert field['provenance']['shared_atmospheric_profile_sha256']==manifest['profiles'][0]['profile_sha256']
    assert gpu['html_sha256']==BASELINE_HTML and not gpu['gpu']['context_lost'];assert mc['unresolved_paths']==0
    browser=read(LABS/'validation/a1/lab-browser.json');assert browser['passed'] and not browser['errors'] and not browser['external_requests'];assert browser['html_sha256']==sha(LABS/'Open_Moon_A1_Accuracy_Lab.html')
    assert sha(VIEWER/'Open_Moon_Shoreline.html')==BASELINE_HTML
    tap=(ROOT/'validation/a1/unit-tests.tap').read_text();counts={k:int(re.search(r'^# '+k+r' (\d+)$',tap,re.M).group(1)) for k in ('tests','pass','fail','skipped')};assert counts['tests']==counts['pass'] and counts['fail']==counts['skipped']==0
    paths=[LABS/'build_accuracy.py',LABS/'index.accuracy-lab.html',LABS/'src/accuracy-lab.js',COLUMN/'atmospheric-profile.js',ROOT/'src/profile-optics.js',ROOT/'src/frozen-cloud-field.js',COLUMN/'weather-column.js',VIEWER/'src/cloud-renderer.js',ROOT/'tests/a1_frozen_probe.js',ROOT/'tests/atmospheric-profile.test.cjs',ROOT/'tests/frozen-cloud-field.test.cjs',ROOT/'data/a1/spectral-inputs.json',ROOT/'data/a1/BRUNETON_LICENSE.txt',ROOT/'A1_METHODS.md',ROOT/'A1_THIRD_PARTY_NOTICES.md',ROOT/'requirements.txt',*sorted((ROOT/'tools/a1').glob('*.py')),*sorted((ROOT/'tools/a1').glob('*.cjs'))]
    evidence=ROOT/'validation/a1'
    for filename,obj in [('frozen-field-gpu.json',gpu),('cloud-reference.json',mc),('profile-manifest.json',manifest),('clear-sky-numerics.json',[{'id':p['id'],**{k:read(GEN/(p['id']+'.sky.json'))[k] for k in ('profile_sha256','solver_sha256','spectral_sha256','model','budgets','column_validation','ray_convergence')}} for p in manifest['profiles']])]:
        (evidence/filename).write_text(json.dumps(obj,indent=2)+'\n')
    output={'schema':'open-moon-a1-validation/1','generated_utc':datetime.now(timezone.utc).isoformat(),'status':'reference-checkpoint-passed','scope':'Additive A1 reference implementation, offline inspector and cloud-only transport benchmark. Original shoreline rendering remains unchanged.',
    'publication':{'remote_write_performed':False,'baseline':'Saved revision-04 source built on 032765c51b1aeb140c732e8f484e47efda15cf4d','patch_application':'Apply only after saved revision 04, after git apply --check; do not replace newer files.'},
    'baseline_html_sha256':BASELINE_HTML,'source_sha256':{str(p.relative_to(REPO)):sha(p) for p in paths},'runtime':{'node':subprocess.check_output(['node','--version'],text=True).strip(),'python_packages':{x:importlib.metadata.version(x) for x in ['numpy','scipy','numba','playwright']}},
    'unit_tests':counts,'clear_sky':{'states':9,'tested_column_geometries_per_state':4,'sampled_radiance_convergence_rays_per_state':6,'checks':sky_checks,'scattering_orders':1,'maximum_tested_relative_spectral_L1':max(x['six_ray_max_relative_spectral_L1'] for x in sky_checks),'whole_atlas_error_bound':None,'physical_error_bound':None},
    'fair_moon_column':{'vertical_reference_density_metres':manifest['profiles'][0]['vertical']['air_m'],'legacy_exponential_metres':1.2*48400,'relative_change':manifest['profiles'][0]['vertical']['air_m']/(1.2*48400)-1,'upper_sensitivity':manifest['profiles'][0]['upperSensitivity']},
    'frozen_cloud':{'field_sha256':field_sha,'dimensions':field['dimensions'],'gpu_backend':gpu['gpu']['renderer'],'gpu_cpu_max_absolute_extinction_error_m1':gpu['field_parity']['maximum_absolute_extinction_error_m1'],'voxelization_24_probe_relative_L1':gpu['voxelization']['sum_absolute_difference_over_sum_continuous'],'voxelization_gate_closed':False,'samples_per_ray_per_order':mc['samples_per_ray_per_order'],'ray_count':len(mc['rays']),'analytic_checks':mc['analytic_tests'],'unresolved_paths':mc['unresolved_paths'],'first_order_gpu_mc_agreement':all(x['single_order_agrees_with_GPU_within_6SE_plus_1e_minus5'] for x in mc['rays']),'physical_comparison_scope':mc['model'],'production_radiance_error_bound':None},
    'lab_browser':browser,'open_gates':['Higher-order molecular scattering and independent full-sky calibration','Species-resolved optical inputs and spectral-bin convergence','Cloud voxel and crop convergence','Full production cloud closure versus coupled gas/cloud reference','Interactive angular-cache, time and position error','Production integration of shared atmospheric state','Entraining, mixed-phase and conserved-condensate cloud model','Sustained traversal and consumer-hardware measurements'],
    'initial_setup_failure':'Before building the required baseline HTML, the initial 137-test run failed only the missing-build fixture; 136 passed. The final 159-test suite passes after the prescribed build.',
    'generated_artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in [LABS/'Open_Moon_A1_Accuracy_Lab.html',GEN/'profiles.json',GEN/'moon-fair.gpu.json',GEN/'moon-fair.mc.json',GEN/'moon-fair.cloud.json']}}
    (ROOT/'A1_VALIDATION.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'status':output['status'],'tests':counts,'sky_max_tested_relative_L1':output['clear_sky']['maximum_tested_relative_spectral_L1'],'cloud_grid_gate_closed':False},indent=2))
if __name__=='__main__':main()
