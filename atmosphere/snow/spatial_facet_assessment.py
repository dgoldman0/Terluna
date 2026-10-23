"""Earth comparisons for spatial facet transport and published 2023 kinetics.

Observed shapes are supplied. Conditional rate integrals are not autonomous
trajectories, and a valid algebraic source balance does not guarantee feasible
uptake over the rest of a facet. No parameters are fitted to these observations.
"""
from __future__ import annotations
import itertools
import json
import math
from pathlib import Path
import numpy as np
from atmosphere.snow import spatial_facet_growth as model, spatial_diffusion as spatial
from atmosphere.snow import facet_attachment, facet_growth, hollow_growth_assessment as hollow, hollow_shape
from atmosphere.snow import flakes, takahashi_growth as lab, crystal_transport_replay as replay

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RESULT=HERE/'results/spatial_facet_assessment.json'
SOURCES=HERE/'spatial_growth_sources.json'


def laws_for(state,tc):
    return dict(inherited=[spatial.Attachment('harrington',critical=state['characteristic_supersaturation_a']),
                           spatial.Attachment('harrington',critical=state['characteristic_supersaturation_c'])],
        dislocation_prism=[spatial.Attachment('harrington',critical=state['characteristic_supersaturation_a'],mechanism=1),
                           spatial.Attachment('harrington',critical=state['characteristic_supersaturation_c'])],
        M1=model.libbrecht_laws(tc),broad=model.libbrecht_laws(tc,narrow=False),
        M1_natural_log=model.libbrecht_laws(tc,logarithm='ln'))


def geometry(a,c,rim,depth,n=32,extent=80.):
    return model.FacetGeometry(spatial.Domain(c/a,resolution=n,extent=extent,mesh='edge',
        cavity_radius=1-rim,cavity_depth=depth*c/a if rim<1 else 0.))


def calculated(g,state,a,t,si,laws,**kwargs):
    try:return dict(complete=True,calculation=g.solve(state,a,t,si,laws,**kwargs))
    except (ValueError,ArithmeticError) as exc:return dict(complete=False,error=str(exc))


def rates(answer):
    return np.array([answer['prism_rate_m_s'],answer['basal_rate_m_s'],answer['mass_rate_kg_s']])


def warm_comparison(library):
    rows=[];inputs=[]
    tc=-5.3;t=flakes.T0+tc
    for age in (4.5,8.,27.):
        obs=lab.observed(tc,age);a=obs['equatorial_dimension_m']/2;c=obs['polar_dimension_m']/2
        state=library.evaluate(.1,t,101000.,a,c,10.)
        si=flakes.liquid_saturation_pressure(t)/state['ice_saturation_pressure_Pa']-1
        nm,na,nc=replay.fitted_exponents(tc,age);reference=np.array([na*a,nc*c,nm*obs['mass_kg']])/(60*age)
        air=replay.air_state(tc);nu=air['viscosity_Pa_s']/air['density_kg_m3']
        fv=replay.jw_ventilation(lab.fall_speed(tc,age)*2*a/nu,.63,replay.branch_habit(tc,age))
        inputs.append(dict(age_min=age,temperature_C=tc,pressure_Pa=101000.,observed_fit=obs,
            fitted_rates=reference.tolist(),state=state,ice_supersaturation=si,JW_factor=fv))
        for rim in (1.,.3,.1,.03,.01):
            g=geometry(a,c,rim,.5)
            for name,law in laws_for(state,tc).items():
                for sources in itertools.product(('edge','middle'),repeat=2):
                    answer=calculated(g,state,a,t,si,law,sources=sources,vapor_ventilation=fv,heat_ventilation=fv)
                    row=dict(age_min=age,rim_fraction=rim,cavity_depth_fraction=.5 if rim<1 else 0.,
                        law=name,sources=list(sources),resolution=32,extent=80.,**answer)
                    if answer['complete']:row['rates_to_fitted_derivatives']=(rates(answer['calculation'])/reference).tolist()
                    rows.append(row)
    return dict(inputs=inputs,cases=rows,
        limits=['Rims and cavity depth are unmeasured sensitivity inputs; no choice is a reconstruction.',
            'Independent mass and axis fits have no supplied derivative errors.',
            'JW is a common conductance multiplier; the flow around the actual hollow crystal is not resolved.',
            'Cylinder normal rates are compared with fitted half-maximum spans; hexagonal coordinate and profile effects remain.',
            'Changing Harrington mechanism while retaining its critical supersaturation is an explicit sensitivity.'])


def chamber_series(library):
    members=hollow_shape.archive_members();records=[];integrals=[]
    for case in hollow.chamber_controls(members)['cases']:
        _,data=hollow.parse_chamber(members[case['member']]);tc=case['temperature_K']-flakes.T0
        t=case['temperature_K'];switch=case['switch_elapsed_s']
        for n in (16,32,64):
            for row in data:
                a,c,rim=row[1]*1e-6,row[2]*1e-6,row[7]/row[1]
                g=geometry(a,c,rim,.5,n=n)
                si=.2 if switch is not None and row[0]>=switch else .48
                state=library.evaluate(si,t,case['pressure_Pa'],a,c,10.)
                for name in ('inherited','M1','broad'):
                    law=laws_for(state,tc)[name]
                    for sources in (('edge','edge'),('middle','edge'),('middle','middle')):
                        answer=calculated(g,state,a,t,si,law,sources=sources)
                        records.append(dict(member=case['member'],elapsed_s=float(row[0]),resolution=n,
                            a_m=a,c_m=c,rim_fraction=rim,ice_supersaturation=si,law=name,sources=list(sources),**answer))
            for win in case['windows']:
                start,end=win['requested_interval_min']
                selected=data[(data[:,0]>=60*start)&(data[:,0]<60*end)]
                # These source windows lie wholly on one side of the switch.
                # Integrate to their first/last actual observation, not nominal endpoints.
                for name,sources in itertools.product(('inherited','M1','broad'),(('edge','edge'),('middle','edge'),('middle','middle'))):
                    rr=[x for x in records if x['member']==case['member'] and x['resolution']==n and
                        x['law']==name and x['sources']==list(sources) and selected[0,0]<=x['elapsed_s']<=selected[-1,0]]
                    item=dict(member=case['member'],law=name,sources=list(sources),resolution=n,
                        requested_interval_min=[start,end],actual_interval_s=[float(selected[0,0]),float(selected[-1,0])],
                        observed_endpoint_increment_um=(selected[-1,1:3]-selected[0,1:3]).tolist(),
                        endpoint_error_bound_um=(selected[-1,5:7]+selected[0,5:7]).tolist(),
                        complete=all(x['complete'] for x in rr),states=len(rr))
                    if item['complete']:
                        values=np.array([rates(x['calculation'])[:2] for x in rr]);times=np.array([x['elapsed_s'] for x in rr])
                        integral=np.trapezoid(values,times,axis=0)*1e6
                        item.update(conditional_increment_um=integral.tolist(),
                            increment_to_observed=(integral/np.array(item['observed_endpoint_increment_um'])).tolist(),
                            feasible_states=sum(x['calculation']['passes_uniform_local_uptake_conditions'] for x in rr),
                            minimum_surface_supersaturation=min(x['calculation']['minimum_surface_supersaturation'] for x in rr))
                    integrals.append(item)
    return dict(cases=records,conditional_integrals=integrals,
        source_members=[dict(member=x['member'],sha256=x['member_sha256'],temperature_K=x['temperature_K'],
            pressure_Pa=x['pressure_Pa'],switch_elapsed_s=x['switch_elapsed_s']) for x in hollow.chamber_controls(members)['cases']],
        limits=['Observed dimensions and rim are supplied at each time; the calculation predicts neither geometry nor rim evolution.',
            'The fixed cavity depth is half the observed semilength; cavity interior is inactive.',
            'Trapezoidal integration uses observation times; mesh refinement is separate from unknown temporal interpolation error.',
            'Supplied endpoint-error sums bound measurement perturbations, not forcing, projection or substrate uncertainty.',
            'Isothermal free-crystal heat loss excludes substrate conduction; no air circulation is modeled.'])


def refinements(library):
    rows=[]
    # Follow the previously nearly feasible M1 rim; refine its corners explicitly.
    tc=-5.3;age=8.;t=flakes.T0+tc;obs=lab.observed(tc,age)
    a=obs['equatorial_dimension_m']/2;c=obs['polar_dimension_m']/2
    state=library.evaluate(.1,t,101000.,a,c,10.);si=flakes.liquid_saturation_pressure(t)/state['ice_saturation_pressure_Pa']-1
    nm,na,nc=replay.fitted_exponents(tc,age);ref=np.array([na*a,nc*c,nm*obs['mass_kg']])/(60*age)
    air=replay.air_state(tc);fv=replay.jw_ventilation(lab.fall_speed(tc,age)*2*a*air['density_kg_m3']/air['viscosity_Pa_s'],.63,replay.branch_habit(tc,age))
    for rim in (.1,.03,.01):
        for n in (16,32,64):
            g=geometry(a,c,rim,.5,n=n)
            for name,sources in [('M1',('edge','edge')),('inherited',('middle','middle')),('broad',('middle','edge'))]:
                answer=calculated(g,state,a,t,si,laws_for(state,tc)[name],sources=sources,vapor_ventilation=fv,heat_ventilation=fv)
                if answer['complete']:answer['rates_to_fit']=(rates(answer['calculation'])/ref).tolist()
                rows.append(dict(rim_fraction=rim,law=name,sources=list(sources),resolution=n,cells=g.domain.cells,**answer))
    return rows


def run():
    sources=json.loads(SOURCES.read_text())
    for x in sources['retrieved']:
        path=ROOT/x['cache_path']
        if hollow_shape.sha(path)!=x['sha256'] or path.stat().st_size!=x['bytes']:raise ValueError('Source identity mismatch: '+x['name'])
    replay.verify_sources();hollow.read_sources()
    paths=[Path(__file__),Path(model.__file__),Path(spatial.__file__),Path(facet_growth.__file__),Path(lab.__file__),Path(replay.__file__),Path(hollow.__file__),SOURCES]
    producer={str(p.relative_to(ROOT)):hollow_shape.sha(p) for p in paths}
    library=facet_attachment.AttachmentLibrary(facet_growth.BUILD)
    warm=warm_comparison(library);print('Warm snapshots complete',flush=True)
    cold=chamber_series(library);print('Chamber prescribed-shape series complete',flush=True)
    refinement=refinements(library)
    complete=[x for x in warm['cases']+cold['cases']+refinement if x['complete']]
    summary=dict(warm_requests=len(warm['cases']),cold_requests=len(cold['cases']),refinement_requests=len(refinement),
        all_completed=len(complete),failed_requests=sum(not x['complete'] for x in warm['cases']+cold['cases']+refinement),
        feasible_warm=sum(x['complete'] and x['calculation']['passes_uniform_local_uptake_conditions'] for x in warm['cases']),
        feasible_cold=sum(x['complete'] and x['calculation']['passes_uniform_local_uptake_conditions'] for x in cold['cases']),
        maximum_source_balance_residual=max(x['calculation']['source_balance_relative_residual'] for x in complete),
        maximum_transport_mass_residual=max(x['calculation']['transport_relative_mass_residual'] for x in complete))
    result=dict(schema='terluna.atmosphere.spatial-facet-assessment/1',producer=producer,sources=sources,
        attachment_build=facet_attachment.verify_build(facet_growth.BUILD),
        warm=warm,cold=cold,refinements=refinement,summary=summary,
        evidence='Prescribed-shape Earth tests and transport feasibility; no fitted parameters, autonomous habit validation or lunar size forecast.',
        limitations=['Uniform local uptake is a hypothesis. A failed check does not exclude uniform facet advance sustained by nonlocal surface transport.',
            'M1 is a proposed narrow-facet law, not a complete width-dependent model. The broad limit and printed-log ambiguity are kept separate.',
            'The 2023 photos are selected examples, not an unbiased population distribution; no pixel data are inferred here.'])
    for name,digest in producer.items():
        if hollow_shape.sha(ROOT/name)!=digest:raise RuntimeError('Producer changed during run: '+name)
    RESULT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))
    return result


if __name__=='__main__':run()
