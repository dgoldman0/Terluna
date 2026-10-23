"""Independent analytic and source-geometry controls for spatial diffusion."""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
import math
import numpy as np
from atmosphere.snow.spatial_diffusion import Domain, Attachment
from atmosphere.snow.hollow_shape import ShapeLibrary

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULT=HERE/'results/spatial_diffusion_assessment.json'


def file_id(path):
    data=path.read_bytes()
    return dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))


def corner(domain,values,face):
    ids=np.where(domain.ice[:,3]==face)[0]
    k=ids[np.argmin((domain.ice[ids,4]-1)**2+(domain.ice[ids,5]-domain.aspect)**2)]
    return float(values[k])


def run():
    sources=ShapeLibrary(); spheres=[]; cylinders=[]; cavities=[]; rings=[]
    for extent in (20.,80.):
        # Exact spherical Robin capacities at the inscribed/circumscribed
        # outer spheres bracket the cylindrical-domain result in the continuum.
        for n in (16,32,64):
            d=Domain(resolution=n,extent=extent,shape='sphere')
            cap=d.dirichlet()
            for lam in (.1,10.,1000.):
                answer=d.solve(.1,lam,[Attachment('constant')]*2,tolerance=1e-12)
                normalized=answer['half_flux']/(2*math.pi*.1)
                bounds=[1/(1+1/lam-1/(math.sqrt(2)*extent)),1/(1+1/lam-1/extent)]
                spheres.append(dict(resolution=n,extent=extent,lambda_=lam,
                    normalized_flux=normalized,analytic_outer_sphere_bounds=bounds,
                    lies_within_continuum_bounds=bounds[0]<=normalized<=bounds[1],
                    capacity=cap,relative_mass_residual=answer['relative_mass_residual'],
                    residual=answer['nonlinear_relative_residual']))
    for phi in (.25,1.,2.5,10.):
        h=sources.evaluate(phi,'cylinder'); levels=[]
        for n in (32,64,128):
            d=Domain(aspect=phi,resolution=n,extent=80.)
            rows=[]
            for face in (0,1):
                solved=d.imposed_flux(prism=1-face,basal=face)
                v=solved['surface']; mids={}
                for name,r,z in [('prism_middle',1.,0.),('basal_middle',0.,phi)]:
                    k=np.argmin((d.ice[:,4]-r)**2+(d.ice[:,5]-z)**2)
                    mids[name]=float(v[k])
                rows.append(dict(face=face,corner_sample=corner(d,v,face),middle_samples=mids,
                    relative_mass_residual=solved['relative_mass_residual']))
            levels.append(dict(resolution=n,cells=d.cells,faces=rows))
        # Reentrant exterior angle 3*pi/2 gives leading corner dependence h^(2/3).
        # This is an own asymptotic extrapolation, retained separately from raw
        # samples and checked by the preceding refinement pair.
        extrapolations=[]
        for low,high in zip(levels[:-1],levels[1:]):
            extrapolations.append([(2**(2/3)*high['faces'][f]['corner_sample']-low['faces'][f]['corner_sample'])/(2**(2/3)-1) for f in (0,1)])
        reference=[h['prism_edge'],phi*h['basal_edge']]
        cylinders.append(dict(aspect=phi,outer_extent=80.,levels=levels,corner_extrapolations=extrapolations,
            published_fitted_edge_lengths=reference,
            relative_extrapolated_difference=[extrapolations[-1][f]/reference[f]-1 for f in (0,1)],
            published_middle_lengths=[[h['prism_middle'],h['prism_at_basal_middle']],
                [phi*h['basal_at_prism_middle'],phi*h['basal_middle']]]))
    for rim in (.0625,.125,.25,.5):
        levels=[]
        for n in (32,64,128):
            d=Domain(resolution=n,extent=80.)
            support=(d.ice[:,3]==1)&(d.ice[:,4]>1-rim)
            answer=d.imposed_flux(face_flux=support.astype(float))
            levels.append(dict(resolution=n,corner_sample=corner(d,answer['surface'],1),
                full_area_over_radius_squared=2*float(np.sum(d.ice[support,1])),
                relative_mass_residual=answer['relative_mass_residual']))
        extrap=(2**(2/3)*levels[-1]['corner_sample']-levels[-2]['corner_sample'])/(2**(2/3)-1)
        asym=sources.evaluate(1.,'cylinder')['basal_edge']-sources.evaluate(1/(1-rim),'cylinder')['basal_edge']
        rings.append(dict(rim_fraction=rim,levels=levels,extrapolated_edge_length=extrap,
            original_thin_ring_length=asym,thin_ring_to_spatial=asym/extrap,
            within_previous_declared_thin_ring_regime=rim<=.1,
            geometry='flat masked end face; cavity depth is zero, not a growing hollow profile'))
    laws=[Attachment('nucleation',.2,.002),Attachment('nucleation',1.,.007)]
    for rim in (.1,.3,1.):
        levels=[]
        for n in (12,24,48):
            d=Domain(aspect=2.,resolution=n,extent=40.,cavity_radius=1-rim,cavity_depth=1. if rim<1 else 0.)
            answer=d.solve(.05,1000.,laws,heat_parameter=.5,tolerance=1e-12)
            levels.append({k:v for k,v in answer.items() if not isinstance(v,np.ndarray)})
        cavities.append(dict(rim_fraction=rim,aspect=2.,cavity_depth_over_radius=1. if rim<1 else 0.,
            lambda_=1000.,supersaturation=.05,heat_parameter=.5,levels=levels,
            fine_to_medium_flux_change=levels[-1]['half_flux']/levels[-2]['half_flux']-1))
    return dict(schema='terluna.atmosphere.spatial-diffusion-assessment/1',
        evidence='numerical analytic and source-geometry controls; no empirical growth validation',
        inputs=[file_id(HERE/'spatial_diffusion.py'),file_id(Path(__file__).resolve()),file_id(HERE/'spatial_growth_sources.json')],
        sphere_controls=spheres,cylinder_controls=cylinders,ring_controls=rings,cavity_controls=cavities)


if __name__=='__main__':
    result=run();RESULT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(sphere_cases=len(result['sphere_controls']),
        refined_spheres_in_bounds=all(x['lies_within_continuum_bounds'] for x in result['sphere_controls'] if x['resolution']==64),
        cylinder_max_fit_difference=max(abs(y) for x in result['cylinder_controls'] for y in x['relative_extrapolated_difference']),
        cavity_max_flux_refinement=max(abs(x['fine_to_medium_flux_change']) for x in result['cavity_controls']),
        ring_ratios=[x['thin_ring_to_spatial'] for x in result['ring_controls']]),indent=2))
