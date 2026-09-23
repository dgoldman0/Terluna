"""Build the A2 full angular sky atlas from a saved molecular scattering field."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sp=importlib.util.spec_from_file_location('a2mol',ROOT/'tools/a2/molecular_multiple_scattering.py');MOL=importlib.util.module_from_spec(sp);sp.loader.exec_module(MOL)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--field',type=Path,default=ROOT/'data/a2/generated/moon-fair-molecular-standard-full.npz');ap.add_argument('--profile',type=Path,default=ROOT/'data/a1/generated/moon-fair.profile.json');ap.add_argument('--a1-spectra',type=Path,default=ROOT/'data/a1/generated/moon-fair.spectra.npz');ap.add_argument('--output',type=Path,default=ROOT/'data/a2/generated/moon-fair.multiple-sky.json');ap.add_argument('--spectra-output',type=Path,default=ROOT/'data/a2/generated/moon-fair.multiple-sky.npz');a=ap.parse_args();start=time.monotonic()
 state,rows,oz,_,profile_sha=MOL.a1.unpack_profile(a.profile);d=np.load(a.field);lam=d['lam'];rg=d['rg'];ag=d['ag'];srg=d['srg'];sag=d['sag'];beam=d['beam'];mom=d['source_mom'];beta=d['beta'];sigma=d['sigma'];zero=np.zeros_like(mom)
 if str(d['profile_sha256'])!=profile_sha:raise ValueError('Field/profile mismatch')
 spec=json.loads((ROOT/'data/a1/spectral-inputs.json').read_text());w=np.array(spec['wavelength_nm'],float)
 if not np.array_equal(lam,w):raise ValueError('Field does not cover exact spectral input grid')
 cie=np.array(spec['cie1931_xyz'],float);rgbmat=np.array(spec['xyz_to_linear_srgb'],float);dw=np.full(len(lam),10.);dw[[0,-1]]=5.;xyz_w=cie*dw[:,None]*683.
 suns=np.array([-6.,-2.,0.,2.,6.,15.,45.,90.]);elev=90*(np.arange(25)/24)**2;az=np.linspace(0,180,33);single=np.zeros((len(suns),len(elev),len(az),len(lam)));multiple=np.zeros_like(single);shells=MOL.shell_altitudes(state['upper']['top_m'],160,oz,state['sounding']['end_m'])
 for si,sun in enumerate(suns):
  for ei,el in enumerate(elev):
   for ai,azi in enumerate(az):
    args=(math.radians(el),math.radians(azi),math.radians(sun),state['planet']['radius'],state['upper']['top_m'],rows,oz,beta,sigma,beam,rg,ag,srg,sag)
    single[si,ei,ai]=MOL.eval_ray(*args,zero,shells);multiple[si,ei,ai]=MOL.eval_ray(*args,mom,shells)
  print(json.dumps({'sun_deg':sun,'elapsed_s':time.monotonic()-start}),flush=True)
 a1=np.load(a.a1_spectra);ref=a1['radiance']
 if ref.shape!=single.shape or not np.array_equal(a1['wavelength_nm'],lam) or str(a1['profile_sha256'])!=profile_sha:raise ValueError('A1 comparison mismatch')
 delta=np.abs(single-ref);global_l1=float(delta.sum()/max(1e-30,np.abs(ref).sum()));den=np.abs(ref).sum(axis=-1);per=np.divide(delta.sum(axis=-1),den,out=np.zeros_like(den),where=den>1e-12);valid=den>1e-8
 single_xyz=np.einsum('...l,lc->...c',single,xyz_w);multiple_xyz=np.einsum('...l,lc->...c',multiple,xyz_w);single_rgb=single_xyz@rgbmat.T;multiple_rgb=multiple_xyz@rgbmat.T
 ratio=np.divide(multiple_xyz[...,1],single_xyz[...,1],out=np.full(single_xyz.shape[:-1],np.nan),where=single_xyz[...,1]>1e-9)
 np.savez_compressed(a.spectra_output,wavelength_nm=lam,suns_deg=suns,elevations_deg=elev,azimuths_deg=az,single=single,multiple=multiple,profile_sha256=np.array(profile_sha),field_sha256=np.array(hashlib.sha256(a.field.read_bytes()).hexdigest()))
 out={'schema':'open-moon-a2-multiple-sky/1','profile_sha256':profile_sha,'field_file':a.field.name,'field_sha256':hashlib.sha256(a.field.read_bytes()).hexdigest(),'solver_sha256':hashlib.sha256((ROOT/'tools/a2/molecular_multiple_scattering.py').read_bytes()).hexdigest(),'atlas_builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'suns_deg':suns.tolist(),'elevations_deg':elev.tolist(),'azimuths_deg':az.tolist(),'single_linear_srgb':single_rgb.tolist(),'multiple_linear_srgb':multiple_rgb.tolist(),'single_XYZ':single_xyz.tolist(),'multiple_XYZ':multiple_xyz.tolist(),'comparison_to_a1_first_order':{'global_spectral_L1':global_l1,'median_valid_pixel_spectral_L1':float(np.median(per[valid])),'p95_valid_pixel_spectral_L1':float(np.quantile(per[valid],.95)),'max_valid_pixel_spectral_L1':float(per[valid].max()),'valid_pixel_threshold_sum_radiance':1e-8},'multiple_to_single_luminance_Y':{'median_for_single_Y_gt_1e_minus_9':float(np.nanmedian(ratio)),'p05':float(np.nanquantile(ratio,.05)),'p95':float(np.nanquantile(ratio,.95)),'max':float(np.nanmax(ratio))},'model':'Standard-grid spherical scalar Rayleigh successive-orders field, 48 wavelengths, black surface, A1 shared profile and ozone, no aerosol/refraction/cloud.','limitations':['Standard-grid atlas; fine-grid validation is supplied on seven selected rays in the chunked fine reference.','Displayed RGB requires a separate exposure/tone curve; saved values remain linear.','Specified atmospheric state remains conditional.'],'elapsed_s':time.monotonic()-start}
 a.output.write_text(json.dumps(out,separators=(',',':'))+'\n');print(json.dumps({'output':str(a.output),'spectra':str(a.spectra_output),'comparison':out['comparison_to_a1_first_order'],'luminance_ratio':out['multiple_to_single_luminance_Y'],'elapsed_s':out['elapsed_s']},indent=2))
if __name__=='__main__':main()
