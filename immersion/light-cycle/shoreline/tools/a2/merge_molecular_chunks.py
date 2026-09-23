"""Merge independent A2 wavelength chunks and quantify standard/fine convergence."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
CHUNKS=['360-410','420-470','480-530','540-590','600-650','660-710','720-770','780-830']
def l1(a,b):return float(np.sum(np.abs(a-b))/max(1e-30,np.sum(np.abs(b))))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=Path,default=ROOT/'data/a2/generated');ap.add_argument('--standard',type=Path,default=ROOT/'data/a2/generated/moon-fair-molecular-standard-full.json');ap.add_argument('--output',type=Path,default=ROOT/'data/a2/generated/moon-fair-molecular-fine-full-chunked.json');a=ap.parse_args()
 docs=[]
 for c in CHUNKS:
  p=a.directory/f'mol-fine-{c}.json';d=json.loads(p.read_text());docs.append((p,d,hashlib.sha256(p.read_bytes()).hexdigest()))
 profile={d['profile_sha256'] for _,d,_ in docs};specsha={d['spectral_inputs_sha256'] for _,d,_ in docs}
 if len(profile)!=1 or len(specsha)!=1:raise ValueError('Chunk identity mismatch')
 nr=len(docs[0][1]['rays'])
 if any(len(d['rays'])!=nr for _,d,_ in docs):raise ValueError('Ray count mismatch')
 lam=[]
 for _,d,_ in docs:lam.extend(d['wavelength_nm'])
 if lam!=sorted(lam) or len(set(lam))!=len(lam):raise ValueError('Wavelength chunks overlap or unordered')
 spec=json.loads((ROOT/'data/a1/spectral-inputs.json').read_text());allw=spec['wavelength_nm']
 if lam!=allw:raise ValueError('Chunk merge does not cover the exact A1 spectral grid')
 cie=np.array(spec['cie1931_xyz'],float);M=np.array(spec['xyz_to_linear_srgb'],float);weights=np.full(len(lam),10.);weights[[0,-1]]=5.;xyz_weights=cie*weights[:,None]*683.
 standard=json.loads(a.standard.read_text());si={w:i for i,w in enumerate(standard['wavelength_nm'])}
 rays=[];std_num=std_den=0.
 for ri in range(nr):
  singles=[];multis=[];stdmulti=[];validation=[];meta=None
  for p,d,sha in docs:
   r=d['rays'][ri];cur=(r['sun_deg'],r['view_elevation_deg'],r['view_azimuth_deg'])
   if meta is None:meta=cur
   elif meta!=cur:raise ValueError('Ray geometry mismatch')
   singles.extend(r['single']);multis.extend(r['multiple_through_final_order']);stdmulti.extend(standard['rays'][ri]['multiple_through_final_order'][si[w]] for w in d['wavelength_nm']);validation.append({'chunk':p.name,'max_bandset_L1_vs_A1':r['single_relative_spectral_L1_vs_a1_validation_bands']})
  single=np.array(singles);multi=np.array(multis);sm=np.array(stdmulti);single_xyz=single@xyz_weights;multi_xyz=multi@xyz_weights;std_xyz=sm@xyz_weights;single_rgb=single_xyz@M.T;multi_rgb=multi_xyz@M.T;std_rgb=std_xyz@M.T
  err=l1(sm,multi);std_num+=np.sum(np.abs(sm-multi));std_den+=np.sum(np.abs(multi))
  rays.append({'sun_deg':meta[0],'view_elevation_deg':meta[1],'view_azimuth_deg':meta[2],'single':single.tolist(),'multiple':multi.tolist(),'single_XYZ':single_xyz.tolist(),'multiple_XYZ':multi_xyz.tolist(),'single_linear_srgb':single_rgb.tolist(),'multiple_linear_srgb':multi_rgb.tolist(),'multiple_to_single_luminance_Y':float(multi_xyz[1]/single_xyz[1]) if single_xyz[1]>1e-30 else None,'standard_multiple_XYZ':std_xyz.tolist(),'standard_to_fine_spectral_L1':err,'standard_to_fine_luminance_Y_relative':float(abs(std_xyz[1]-multi_xyz[1])/max(1e-30,abs(multi_xyz[1]))),'validation':validation})
 out={'schema':'open-moon-a2-molecular-multiple-chunked/1','profile_sha256':next(iter(profile)),'spectral_inputs_sha256':next(iter(specsha)),'wavelength_nm':lam,'chunks':[{'file':p.name,'sha256':sha,'orders_computed':d['orders_computed'],'final_increment_fraction':d['order_history'][-1]['max_radiance_moment_increment_fraction'],'quality':d['quality']} for p,d,sha in docs],'merge':'Wavelengths are uncoupled in scalar elastic Rayleigh transport; chunk concatenation is algebraically identical to one solve on the same numerical grids.','rays':rays,'convergence':{'standard_file':a.standard.name,'standard_sha256':hashlib.sha256(a.standard.read_bytes()).hexdigest(),'overall_standard_to_fine_spectral_L1':float(std_num/std_den),'max_ray_standard_to_fine_spectral_L1':max(r['standard_to_fine_spectral_L1'] for r in rays),'max_ray_standard_to_fine_luminance_Y_relative':max(r['standard_to_fine_luminance_Y_relative'] for r in rays)},'limitations':['Fine-grid spectral solve is assembled from independent six-wavelength chunks to avoid array-scaling stalls; elastic wavelengths do not couple.','The shortest-wavelength chunk remains the slowest to converge in scattering order.','The profile remains a specified atmospheric state rather than a climate prediction.']}
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'output':str(a.output),'convergence':out['convergence'],'luminance_ratios':[r['multiple_to_single_luminance_Y'] for r in rays],'chunk_orders':[x['orders_computed'] for x in out['chunks']]},indent=2))
if __name__=='__main__':main()
