"""Build a self-contained, offline viewer from calculated skies and local geometry."""
import base64,gzip,json,hashlib,html
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def build():
 worlds={};angles=None;packing=[]
 for key in ['earth','moon','moon_no_ozone']:
  p=ROOT/'data'/f'{key}_atlas.npz';d=np.load(p);rgb=d['rgb'].astype(float);n,h,w,_=rgb.shape
  if angles is None:angles=d['suns'].tolist()
  else:assert angles==d['suns'].tolist()
  assert angles[0]==-90 and angles[-1]==90
  scale=np.maximum(np.max(np.abs(rgb),axis=(1,2,3))/1000,1e-20)
  packed=np.ones((n,h,w,4),np.float16);packed[...,:3]=(rgb/scale[:,None,None,None]).astype(np.float16)
  rec=packed[...,:3].astype(float)*scale[:,None,None,None]
  mask=np.abs(rgb)>np.max(np.abs(rgb),axis=(1,2,3))[:,None,None,None]*1e-6
  err=float(np.max((np.abs(rec-rgb)/np.maximum(np.abs(rgb),1e-8))[mask]));assert err<.001
  meta=json.loads(str(d['metadata']));atm=meta['atmosphere']
  worlds[key]=dict(width=w,height=h,layers=n,scales=scale.tolist(),payload=base64.b64encode(gzip.compress(packed.tobytes(),compresslevel=9,mtime=0)).decode(),direct=d['direct'].tolist(),direct_horizontal=d['direct_horizontal'].tolist(),diffuse=d['diffuse'].tolist(),sh=d['sh'].tolist(),cloud_direct=d['cloud_direct'].tolist(),cloud_diffuse=d['cloud_diffuse'].tolist(),radius=atm['radius_m'],ray_beta_rgb=(1.24062e-6*(np.array([680.,550.,440.])/1000.)**-4*atm['density_scale']).tolist(),atmosphere=atm,calculation=meta)
  packing.append(dict(world=key,maximum_significant_relative_quantization=err,atlas_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 landscape={name:base64.b64encode((ROOT/'data'/'landscape'/f'{name}.png').read_bytes()).decode() for name in ['albedo','normal','geometry','shadow']}
 landscape['scene']=json.loads((ROOT/'data'/'landscape'/'scene.json').read_text())
 checks_path=ROOT/'validation'/'full_cycle_checks.json';checks=json.loads(checks_path.read_text()) if checks_path.exists() else {}
 summary=('The full-cycle optical fields were recomputed at standard source-grid resolution, with the scattering-increment check extended through +90°. The clear-sky output has 174 Sun angles and 81 × 65 directional samples per angle. The earlier numerical study reported a lunar global energy-closure discrepancy up to 4.6%; this issue remains unresolved. Earlier Monte Carlo spot comparisons concern the molecular solver, not the new weather display. New full-cycle and browser checks are recorded in the source bundle. The deeper-night source grid remains exploratory. No empirical weather or future-lunar-atmosphere validation is claimed.')
 obj=dict(format_version=2,sun_angles=angles,worlds=worlds,landscape=landscape,check_summary=summary,validation=checks)
 text=(ROOT/'cycle.template.html').read_text();text=text.replace('__VERTEX__',(ROOT/'scene.vert').read_text()).replace('__FRAGMENT__',(ROOT/'scene.frag').read_text()).replace('__APP__',(ROOT/'cycle.js').read_text()).replace('__LICENSES__',html.escape((ROOT/'LICENSES.txt').read_text())).replace('__DATASET__',json.dumps(obj,separators=(',',':')))
 out=ROOT/'Open_Moon_Full_Cycle.html';out.write_text(text);(ROOT/'validation'/'packing_checks.json').write_text(json.dumps(packing,indent=2))
 rows=['phase_from_noon,lunar_elapsed_hours,earth_elapsed_hours,sun_elevation_deg,earth_clear_lux,moon_300DU_clear_lux,moon_0DU_clear_lux']
 arrays={k:(np.array(v['direct_horizontal'])+np.array(v['diffuse']))@np.array([.2126,.7152,.0722]) for k,v in worlds.items()}
 for phase in np.linspace(0,1,721):
  a=float(np.rad2deg(np.arcsin(np.cos(phase*2*np.pi))));lux=[float(np.interp(a,angles,arrays[k])) for k in worlds];rows.append(','.join(map(str,[phase,phase*29.53059*24,phase*24,a,*lux])))
 (ROOT/'data'/'full_cycle_brightness.csv').write_text('\n'.join(rows)+'\n');print(out,out.stat().st_size,'bytes',flush=True)
 return out
if __name__=='__main__':build()
