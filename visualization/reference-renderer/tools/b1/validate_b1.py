"""Check executed B1 records without turning execution checks into error bounds."""
from __future__ import annotations
import hashlib,json,re,subprocess
from pathlib import Path
import numpy as np
import render_scene as r
ROOT=r.ROOT
EXPECTED={(c,float(s)) for c in ['clear','fair'] for s in [60,30,12,6,0,-2,-6]}|{('high',6.),('high',-2.),('thin',6.),('thick',6.)}
BASELINE_HTML='c7ead1fd83e0d55518d2951e08cd46d31f323280885125382a8013b5d74636d9'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 v=ROOT/'validation/b1';v.mkdir(parents=True,exist_ok=True);gen=ROOT/'data/b1/generated'
 report=dict(schema='open-moon-b1-validation/1',publication_parent='9ad75cfee4fee5497cf6cce08947dca8cc3ef35f',scope='Additive offline scene painter; production shoreline unchanged.',gates=[],frames=[],hemispheres=[],maps=[],controls=[],open_gates=[])
 def gate(name,passed,**details):report['gates'].append(dict(name=name,passed=bool(passed),**details))
 inputs=read(v/'input-integrity.json');gate('Inherited generated input integrity',inputs['passed'] and all(sha(r.b.locate(x['path']))==x['sha256'] for x in inputs['files']),files=len(inputs['files']));report['input_dependencies']=inputs['files']
 tests=read(v/'boundary-tests.json');gate('Boundary/transport implementation checks',tests['failed']==0 and tests['passed']==23,passed_tests=tests['passed'],failed_tests=tests['failed'])
 text=(v/'inherited-tests.tap').read_text();passn=int(re.findall(r'^# pass (\d+)$',text,re.M)[-1]);failn=int(re.findall(r'^# fail (\d+)$',text,re.M)[-1]);gate('Inherited JavaScript regression',passn==159 and failn==0,passed_tests=passn,failed_tests=failn)
 seen={'frames':set(),'hemispheres':set()};total=0;all_hash=True;identity=True;finite=True;closure=0.;unresolved=0;relative_errors=[];dimensions=True;input_identity=True;prepared={}
 for folder in ['final','controls']:
  for p in sorted((gen/folder).glob('*.json')):
   m=read(p);d=np.load(p.with_suffix('.npz'));projection=m['camera']['projection'];kind='controls' if folder=='controls' else ('frames' if projection=='perspective' else 'hemispheres')
   if kind in seen:seen[kind].add((m['scenario'],float(m['sun_elevation_deg'])))
   total+=m['total_paths'];unresolved+=m['unresolved_paths'];identity&=(m['sources']==r.raw_identity())
   immutable={k:val for k,val in m.items() if k not in ['identity_sha256','unresolved_paths','elapsed_s','total_paths','radiance_units','display','sampling','outputs']}
   identity&=hashlib.sha256(json.dumps(immutable,sort_keys=True).encode()).hexdigest()==m['identity_sha256']
   dims=(m['height'],m['width'],3)
   dimensions&=d['XYZ'].shape==dims and d['unit_transfer'].shape==(m['height'],m['width'],m['wavelength_groups']) and m['total_paths']==m['height']*m['width']*m['samples_per_pixel_per_band']*m['wavelength_groups']
   if m['scenario'] not in prepared:
    c,*_=r.prepare(m['scenario']);prepared[m['scenario']]=(c['profile_sha256'],c['field_sha256'],hashlib.sha256(c['ext'].tobytes()).hexdigest(),c['lo'].tolist(),c['hi'].tolist())
   input_identity&=prepared[m['scenario']]==(m['profile_sha256'],m['base_cloud_sha256'],m['cloud_extinction_sha256'],m['cloud_bounds_min_m'],m['cloud_bounds_max_m'])
   all_hash&=all((p.parent/name).is_file() and sha(p.parent/name)==value for name,value in m['outputs'].items())
   finite&=all(np.isfinite(d[k]).all() for k in d.files)
   y=d['XYZ'][:,:,1];surf=d['surface_path_XYZ'][:,:,1];se=d['XYZ_standard_error'][:,:,1]
   closure=max(closure,float(np.max(np.maximum(0,surf-y))),float(np.max(np.maximum(0,-surf))))
   pix=r.display(d['filtered_linear_srgb']);rawpix=r.display(d['linear_srgb'])
   sample=dict(name=p.stem,scenario=m['scenario'],sun_elevation_deg=m['sun_elevation_deg'],width=m['width'],height=m['height'],samples_per_group=m['samples_per_pixel_per_band'],groups=m['wavelength_groups'],spectral_histories=m['total_paths'],unresolved=m['unresolved_paths'],metadata_sha256=sha(p),radiance_archive_sha256=sha(p.with_suffix('.npz')),median_Y=float(np.median(y)),median_Y_sampling_SE=float(np.median(se)),median_per_pixel_relative_Y_SE=float(np.median(se/np.maximum(y,1e-10))),surface_path_fraction_of_sum_Y=float(surf.sum()/max(1e-30,y.sum())),zero_EV_filtered_white_fraction=float(np.mean(np.all(pix>=254,axis=2))),zero_EV_raw_white_fraction=float(np.mean(np.all(rawpix>=254,axis=2))),elapsed_s=m['elapsed_s'])
   relative_errors.append(sample['median_per_pixel_relative_Y_SE']);report[kind].append(sample)
 gate('Complete perspective and hemisphere state sets',all(x==EXPECTED for x in seen.values()),perspective_count=len(report['frames']),hemisphere_count=len(report['hemispheres']),expected_each=len(EXPECTED))
 gate('Black-boundary control present',len(report['controls'])==1 and report['controls'][0]['scenario']=='fair' and report['controls'][0]['sun_elevation_deg']==12)
 gate('Frame source and state identities',identity and input_identity);gate('Raw array dimensions and path counts',dimensions);gate('Frame output hashes',all_hash);gate('Finite raw numerical arrays',finite);gate('Surface-path contribution is a nonnegative subset',closure<1e-7,maximum_violation_cd_m2=closure)
 mapseen=set();mapfinite=True;maphash=True;mapidentity=True;mapclosure=0.;maprange=True
 for p in sorted((gen/'maps').glob('*.json')):
  m=read(p);d=np.load(p.with_suffix('.npz'));mapseen.add((m['scenario'],float(m['sun_elevation_deg'])));n=m['size']**2*m['samples_per_band']*m['spectral_groups'];total+=n;unresolved+=m['unresolved_paths']
  maphash&=sha(p.with_suffix('.npz'))==m['output_sha256'];mapidentity&=all(sha(r.b.locate(name))==value for name,value in m['source_hashes'].items())
  mapfinite&=all(np.isfinite(d[k]).all() for k in d.files)
  mapclosure=max(mapclosure,float(np.max(np.abs(d['total_XYZ']-d['direct_XYZ']-d['diffuse_XYZ']))))
  vis=d['direct_beam_visibility'];maprange&=bool(vis.min()>=-1e-12 and vis.max()<=1+1e-12)
  report['maps'].append(dict(name=p.stem,scenario=m['scenario'],sun_elevation_deg=m['sun_elevation_deg'],size=m['size'],spectral_histories=n,metadata_sha256=sha(p),irradiance_archive_sha256=sha(p.with_suffix('.npz')),median_direct_lux=float(np.median(d['direct_XYZ'][:,:,1])),median_diffuse_lux=float(np.median(d['diffuse_XYZ'][:,:,1])),median_diffuse_sampling_SE_lux=float(np.median(d['diffuse_XYZ_standard_error'][:,:,1])),visibility_min=float(vis.min()),visibility_max=float(vis.max()),unresolved=m['unresolved_paths'],elapsed_s=m['elapsed_s']))
 gate('Complete irradiance-map state set',mapseen==EXPECTED,count=len(mapseen),expected=len(EXPECTED));gate('Map source/output identity and finite values',maphash and mapidentity and mapfinite);gate('Direct plus diffuse irradiance closes',mapclosure<1e-6,maximum_residual=mapclosure);gate('Direct visibility within [0,1]',maprange);gate('All accepted image/map histories complete',unresolved==0,spectral_histories=total,unresolved=unresolved)
 spectra=read(v/'scene-spectral-check.json');gate('Selected 48-band versus 12-group diagnostic',all(x['provisional_gate'] for x in spectra['cases']),cases=len(spectra['cases']),status='reference-limited by Monte Carlo uncertainty; no whole-image spectral qualification')
 report['spectral_comparisons']=[{k:c[k] for k in ['scenario','sun_deg','pixel_fraction','Y_relative_difference','Y_relative_two_SE','provisional_gate','gate_definition']} for c in spectra['cases']]
 report['spectral_check_histories']=spectra['total_paths'];report['accepted_image_map_histories']=total
 browser=read(v/'browser.json');gate('Exact-build browser/pixel checks',browser['passed'] and browser['html_sha256']==sha(ROOT/'Open_Moon_B1_Scene_Painter.html'),state_checks=browser['state_checks'],layouts=len(browser['layouts']),errors=browser['errors'],console_errors=browser['console_errors'],external_requests=browser['external_requests'])
 report['browser_navigation']=browser['navigation_attempts'];report['browser_scope']=browser['load_method'];report['html_sha256']=sha(ROOT/'Open_Moon_B1_Scene_Painter.html')
 viewer=ROOT.parents[1]/'immersion'/'dist'/'experiences'/'shoreline'/'index.html'  # the experience B1 was built beside (pinned build; see README)
 gate('Production shoreline build preserved',viewer.is_file() and sha(viewer)==BASELINE_HTML,sha256=sha(viewer) if viewer.is_file() else None)
 baseline=read(ROOT/'B1_BASELINE_MANIFEST.json');changed=[]
 for name,info in baseline['files'].items():
  p=ROOT/name
  if not p.is_file() or sha(p)!=info['sha256']:changed.append(name)
 gate('Inherited source and evidence snapshot preserved',not changed,files_checked=len(baseline['files']),changed=changed)
 pilots=[]
 for p in sorted((gen/'refinement-pilots').glob('*.npz')):
  new=gen/'final'/p.name
  if not new.exists():continue
  lo=np.load(p);hi=np.load(new);m=read(new.with_suffix('.json'));pm=read(p.with_suffix('.json'))
  pilots.append(dict(name=p.stem,pilot_spp=pm['samples_per_pixel_per_band'],delivered_spp=m['samples_per_pixel_per_band'],pilot_median_Y_SE=float(np.median(lo['XYZ_standard_error'][:,:,1])),delivered_median_Y_SE=float(np.median(hi['XYZ_standard_error'][:,:,1])),raw_XYZ_L1_difference_divided_by_delivered_sum=float(np.abs(lo['XYZ']-hi['XYZ']).sum()/np.maximum(hi['XYZ'].sum(),1e-30)),interpretation='Nested random histories; image-change diagnostic, not an independent error bound.'))
 report['sampling_refinement']=pilots;report['image_pixel_relative_SE_median_range']=[min(relative_errors),max(relative_errors)]
 report['open_gates']=['Pixel/path-sample convergence and rare specular-event variance.','Complete-image 12-group versus 48-bin spectral convergence, including low Sun.','8-by-8 primary-Sun angular quadrature and eight-point direct-map disk quadrature refinement.','Coupled cloud voxel/crop and terrain-resolution/domain sensitivity.','Global energy closure and physical cloud, ozone, upper atmosphere and surface-spectrum assumptions.','Rough/finite-depth water, genuine fog/layered/convection fields, evolving weather.','Production GPU/cache integration and sustained real-time performance.']
 report['implementation_and_delivery_passed']=all(g['passed'] for g in report['gates']);report['whole_image_accuracy_qualified']=False
 (ROOT/'B1_VALIDATION.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(passed=report['implementation_and_delivery_passed'],gates=len(report['gates']),failed=[g['name'] for g in report['gates'] if not g['passed']],frames=len(report['frames']),hemispheres=len(report['hemispheres']),maps=len(report['maps']),controls=len(report['controls']),histories=total),indent=2))
 return 0 if report['implementation_and_delivery_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
