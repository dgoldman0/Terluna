"""Balance checks and archive-version audit for the implicit facet component."""
from __future__ import annotations
import argparse
import hashlib
import io
import itertools
import json
from pathlib import Path
import zipfile
import numpy as np
import scipy
from atmosphere.snow import facet_attachment as attachment, facet_growth as baseline, facet_coupling as coupled

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CACHE=ROOT/'research/runs/snow_scale/facet_coupling'
SOURCES=HERE/'facet_coupling_sources.json'
RESULT=HERE/'results/facet_coupling_audit.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_sources():
    record=json.loads(SOURCES.read_text())
    for item in record['retrieved']:
        path=ROOT/item['cache_path']
        if path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:
            raise ValueError('Coupling source identity differs: '+str(path))
    return record


def archive_array(archive,suffix,**kwargs):
    names=[name for name in archive.namelist() if name.endswith('/'+suffix)]
    if len(names)!=1:
        raise ValueError('Archive member is absent or ambiguous: '+suffix)
    data=archive.read(names[0])
    return np.loadtxt(io.BytesIO(data),**kwargs),dict(member=names[0],bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def archive_audit():
    """Diagnose the archived curve without silently changing its stated mechanism."""
    record=verify_sources()
    item=next(x for x in record['retrieved'] if x['cache_path'].endswith('.zip'))
    rows=[]
    with zipfile.ZipFile(ROOT/item['cache_path']) as archive:
        for axis,critical_percent in [('a',3.),('c',6.)]:
            data,identity=archive_array(archive,'Figure1-'+axis+'-axis.output',delimiter=',')
            if data.shape!=(1000,5):
                raise ValueError('Unexpected archive curve dimensions')
            # The archive calls this Figure1; the later manuscript's Figure2
            # caption supplies the matching critical supersaturations 3%/6%.
            selected=data[500:950]
            selected=selected[(selected[:,1] < .995*selected[:,4]) & (selected[:,1]>0)]
            mechanisms=[]
            for mechanism in (10.,15.):
                inferred=np.array([(ambient/surface-1)/coupled.efficiency(surface,critical_percent,mechanism)
                                   for _,surface,_,_,ambient in selected])
                mechanisms.append(dict(mechanism_m=mechanism,count=len(inferred),
                    implied_resistance_min=float(inferred.min()),implied_resistance_max=float(inferred.max()),
                    implied_resistance_median=float(np.median(inferred)),
                    implied_resistance_relative_span=float((inferred.max()-inferred.min())/np.median(inferred))))
            rows.append(dict(axis=axis,identity=identity,characteristic_supersaturation_percent=critical_percent,
                critical_value_max_relative_check=float(np.max(abs(data[500:,4]/data[500:,0]/critical_percent-1))),
                diagnostics=mechanisms))
    return dict(rows=rows,
        equation='If s_surface[1+B alpha(s_surface)]=s_ambient, B=(s_ambient/s_surface-1)/alpha.',
        evidence='Archive consistency diagnostic, not an empirical measurement or an adopted replacement parameter.',
        interpretation='M=15 yields nearly constant inferred resistance while M=10 does not. The later manuscript labels M=10. Archive versus manuscript version/parameter mismatch is unresolved; these curves cannot certify exact numerical parity.',
        limits=['The mathematical diagnostic does not establish the unarchived generating executable, exact heat treatment or full forcing conventions.',
                'Original archive retained unchanged; no curve relabeled or forced into agreement.'])


def calculate():
    sources=verify_sources()
    library=attachment.AttachmentLibrary(baseline.BUILD)
    original=library.identity
    cases=[]
    for tc,pressure,radius,phi,mechanism,si in itertools.product(
            (-5.3,-14.4,-20.1,-35.),(101000.,55000.),(5e-6,50e-6,500e-6),
            (.1,1.,10.),(1.,10.),(.001,.01,.1)):
        t=273.15+tc;c=radius*phi
        state=library.evaluate(si,t,pressure,radius,c,mechanism)
        native=baseline.mass_and_axis_rates(state,radius,c,t,si)
        implicit=coupled.solve(state,radius,c,t,si,mechanism)
        ventilated=coupled.solve(state,radius,c,t,si,mechanism,vapor_ventilation=2.,heat_ventilation=2.)
        cases.append(dict(temperature_C=tc,pressure_Pa=pressure,equatorial_radius_m=radius,
            thickness_over_span=phi,mechanism_m=mechanism,ice_supersaturation=si,
            native_alpha_a=state['alpha_a'],native_alpha_c=state['alpha_c'],
            native_alpha_floor_active=min(state['alpha_a'],state['alpha_c'])<=1.000001e-6,
            implicit=implicit,imposed_factor_two=ventilated,
            implicit_over_native_mass_rate=implicit['mass_rate_kg_s']/native['mass_rate_kg_s'],
            implicit_over_native_facet_ratio=(implicit['alpha_c']/implicit['alpha_a'])/(state['alpha_c']/state['alpha_a']),
            ventilation_mass_ratio=ventilated['mass_rate_kg_s']/implicit['mass_rate_kg_s']))
    if attachment.verify_build(baseline.BUILD)!=original:
        raise ValueError('Attachment helper changed during audit')
    summaries={}
    for label,subset in [('all',cases),('native_floor_inactive',[x for x in cases if not x['native_alpha_floor_active']])]:
        summaries[label]=dict(count=len(subset))
        for key in ('implicit_over_native_mass_rate','implicit_over_native_facet_ratio','ventilation_mass_ratio'):
            values=np.array([x[key] for x in subset])
            summaries[label][key]=dict(min=float(values.min()),median=float(np.median(values)),max=float(values.max()),
                p05=float(np.quantile(values,.05)),p95=float(np.quantile(values,.95)))
    balance=max(abs(value) for x in cases for kind in ('implicit','imposed_factor_two')
                for value in [*x[kind]['relative_surface_balance_residuals'],x[kind]['relative_thermal_balance_residual']])
    paths=[Path(__file__),Path(coupled.__file__),Path(baseline.__file__),Path(attachment.__file__),attachment.DRIVER,SOURCES]
    return dict(schema='terluna.atmosphere.facet-coupling-audit/1',
        evidence='Numerical balance verification and source-version diagnosis; no empirical validation or lunar growth forecast.',
        producer={str(path.relative_to(ROOT)):sha(path) for path in paths},
        source_record_sha256=sha(SOURCES),verified_sources=sources,
        attachment_build_record_sha256=sha(baseline.BUILD/'build.json'),
        numerics=dict(scipy_version=scipy.__version__,numpy_version=np.__version__,root_tolerance=1e-12),
        cases=cases,summary=summaries,maximum_relative_balance_residual=balance,
        archived_surface_curve_audit=archive_audit(),
        limits=['The factor-two ventilation case is an imposed conductance sensitivity, not a Reynolds/drag prediction.',
                'Characteristic supersaturations, source air properties and spheroid geometry retain their empirical and approximation limits.',
                'The implicit calculation solves the stated linearized heat and facet balances; it is not a spatially resolved crystal-flow solution.',
                'No artificial alpha floor is imposed. Extremely weak growth can differ greatly from the native floor.',
                'No branching-density law, seed selection, finite reservoir or aggregation is selected here.'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=RESULT)
    args=p.parse_args();result=calculate()
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(result['cases']),summary=result['summary'],max_residual=result['maximum_relative_balance_residual']),indent=2))

if __name__=='__main__':main()
