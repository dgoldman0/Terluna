"""Quantify A2 cloud voxel and finite-crop transport convergence.

Input fields are exact GPU readbacks of the revision-04 production density at
controlled resolutions/crops. Deterministic first-order transport is paired with
Monte Carlo multiple-scattering checks on four common central rays.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, time
from pathlib import Path
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[2]
sp=importlib.util.spec_from_file_location('a1_cloud_reference',ROOT/'tools/a1/cloud_reference.py');A1=importlib.util.module_from_spec(sp);sp.loader.exec_module(A1)
@njit(cache=True)
def deterministic_single(origin,d,sun,lo,hi,dims,ext,ice,primary,secondary):
    a,b=A1.box(origin,d,lo,hi)
    if b<=a:return 0.,1.
    ds=(b-a)/primary;tr=1.;L=0.
    for i in range(primary):
        p=origin+d*(a+(i+.5)*ds);sigma,icef=A1.voxel(p,lo,hi,dims,ext,ice);opacity=1.-math.exp(-sigma*ds)
        if opacity>1e-12:
            sa,sb=A1.box(p,sun,lo,hi);tau=0.
            if sb>sa:
                sds=(sb-sa)/secondary
                for j in range(secondary):
                    ss,_=A1.voxel(p+sun*(sa+(j+.5)*sds),lo,hi,dims,ext,ice);tau+=ss*sds
            L+=tr*opacity*A1.phase(np.dot(d,sun),icef)*math.exp(-tau)
        tr*=1.-opacity
    return L,tr
def field(path):
    s,lo,hi,dims,ext,ice,majorant,albedo,sha=A1.load_field(path);return {'path':path,'spec':s,'lo':lo,'hi':hi,'dims':dims,'ext':ext,'ice':ice,'majorant':majorant,'albedo':albedo,'sha':sha}
def l1(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);return float(np.sum(np.abs(a-b))/max(1e-30,np.sum(np.abs(b))))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=Path,default=ROOT/'data/a2/generated');ap.add_argument('--mc-samples',type=int,default=8192);ap.add_argument('--output',type=Path,default=ROOT/'data/a2/generated/cloud-convergence.json');a=ap.parse_args()
    names=['grid-coarse','grid-a1','grid-fine','crop-075','crop-150','crop-200'];F={n:field(a.directory/f'moon-fair-{n}.cloud.json') for n in names};base=F['grid-a1'];E=base['spec']['provenance']['a2']['baseExtent_m'];base_y=base['lo'][1];sun=np.array([math.sqrt(.5),math.sqrt(.5),0.]);rays=[]
    for iz in range(3):
      for ix in range(4):rays.append((np.array([-E+(ix+.5)/4*(2*E),base_y-100,-E+(iz+.5)/3*(2*E)]),np.array([0.,1.,0.])))
    start=time.monotonic();out={'schema':'open-moon-a2-cloud-convergence/1','solver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'field_sha256':{n:F[n]['sha'] for n in names},'base_extent_m':E,'sun_direction':sun.tolist(),'deterministic_budget':{'primary':2048,'secondary':512},'fields':{},'grid_convergence':{},'crop_convergence':{},'mc':{'samples_per_ray':a.mc_samples,'rays':[]},'limitations':['Frozen-grid comparisons inherit the finite sample/capture assumptions recorded in each field.','First-order convergence uses deterministic quadrature over the frozen fields, not the continuous production density.','Multiple-scattering convergence is Monte Carlo with reported standard errors and only four central rays.']}
    for n,f in F.items():
        vals=[]
        for o,d in rays:
            L,tr=deterministic_single(o,d,sun,f['lo'],f['hi'],f['dims'],f['ext'],f['ice'],2048,512);vals.append({'radiance':L,'transmittance':tr})
        out['fields'][n]={'dimensions':f['dims'].tolist(),'bounds_min_m':f['lo'].tolist(),'bounds_max_m':f['hi'].tolist(),'horizontal_spacing_m':float((f['hi'][0]-f['lo'][0])/(f['dims'][0]-1)),'vertical_spacing_m':float((f['hi'][1]-f['lo'][1])/(f['dims'][1]-1)),'first_order':vals}
    ref=out['fields']['grid-fine']['first_order'];rr=[x['radiance'] for x in ref];rt=[x['transmittance'] for x in ref]
    for n in ['grid-coarse','grid-a1','grid-fine']:
        v=out['fields'][n]['first_order'];out['grid_convergence'][n]={'radiance_L1_vs_grid_fine':l1([x['radiance'] for x in v],rr),'transmittance_L1_vs_grid_fine':l1([x['transmittance'] for x in v],rt)}
    ref=out['fields']['crop-200']['first_order'];rr=[x['radiance'] for x in ref];rt=[x['transmittance'] for x in ref]
    for n in ['crop-075','grid-a1','crop-150','crop-200']:
        v=out['fields'][n]['first_order'];out['crop_convergence'][n]={'radiance_L1_vs_crop_200':l1([x['radiance'] for x in v],rr),'transmittance_L1_vs_crop_200':l1([x['transmittance'] for x in v],rt)}
    selected=[rays[i] for i in [1,4,6,9]]
    for ri,(o,d) in enumerate(selected):
        row={'ray_index':[1,4,6,9][ri],'origin':o.tolist(),'direction':d.tolist(),'grid':{},'crop':{}}
        for group,gnames in [('grid',['grid-coarse','grid-a1','grid-fine']),('crop',['crop-075','grid-a1','crop-150','crop-200'])]:
            for j,n in enumerate(gnames):
                f=F[n];mean,se,fail=A1.estimate(o,d,sun,1.,0.,f['albedo'],f['lo'],f['hi'],f['dims'],f['ext'],f['ice'],f['majorant'],True,a.mc_samples,8675309+ri*100+j)
                row[group][n]={'multiple_mean':mean,'standard_error':se,'unresolved':int(fail)}
        out['mc']['rays'].append(row);print(json.dumps({'ray':ri,'elapsed_s':time.monotonic()-start}),flush=True)
    # Summaries compare means; uncertainty remains attached to every point.
    def mc_l1(group,n,refn):return l1([r[group][n]['multiple_mean'] for r in out['mc']['rays']],[r[group][refn]['multiple_mean'] for r in out['mc']['rays']])
    out['mc']['grid_mean_L1_vs_fine']={n:mc_l1('grid',n,'grid-fine') for n in ['grid-coarse','grid-a1','grid-fine']};out['mc']['crop_mean_L1_vs_crop_200']={n:mc_l1('crop',n,'crop-200') for n in ['crop-075','grid-a1','crop-150','crop-200']};out['mc']['unresolved_total']=sum(r[g][n]['unresolved'] for r in out['mc']['rays'] for g in ['grid','crop'] for n in r[g]);out['elapsed_s']=time.monotonic()-start;a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'output':str(a.output),'grid':out['grid_convergence'],'crop':out['crop_convergence'],'mc_grid':out['mc']['grid_mean_L1_vs_fine'],'mc_crop':out['mc']['crop_mean_L1_vs_crop_200'],'elapsed_s':out['elapsed_s']},indent=2))
if __name__=='__main__':main()
