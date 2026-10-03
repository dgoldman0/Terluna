"""Verify complete SWASH records and compare explicit shore-profile cases."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from climate.waves.model import ROOT, sha256
from climate.waves.runup import RUNS, PROFILES
from climate.waves.swash_build import RUNS as BUILD

OUTPUT=Path(__file__).parent/'results/runup.json'


def verify_run(directory,record):
    for filename,digest in {**record['identity']['inputs'],**record['output_sha256']}.items():
        if sha256(directory/filename)!=digest:
            raise ValueError(f'Wave-resolving bytes changed: {directory/filename}')


def verify(path):
    record=json.loads(path.read_text())
    if record['schema']!='terluna.climate.runup-case/1':
        raise ValueError('Unsupported wave-resolving shore case')
    verify_run(path.parent,record['run'])
    if sha256(path.parent/'producer_source.py')!=record['producer_sha256']:
        raise ValueError('Recorded wave-resolving producer changed')
    if sha256(PROFILES)!=record['profile_product_sha256']:
        raise ValueError('The shore-profile product changed')
    if 'water_scenario' in record and sha256(ROOT/record['water_scenario']['path'])!=record['water_scenario']['sha256']:
        raise ValueError('The water-density scenario changed')
    source=ROOT/record['boundary']['source_path']
    if sha256(source)!=record['boundary']['source_sha256']:
        raise ValueError('The forcing spectrum changed')
    if ('source_record_sha256' in record['boundary'] and
            sha256(source.parent/'product.json')!=record['boundary']['source_record_sha256']):
        raise ValueError('The forcing spectrum manifest changed')
    upstream=json.loads((source.parent/'product.json').read_text())
    record['purpose']=('Numerical reference using a stationary coastal sea state'
        if upstream['schema']=='terluna.climate.shore-case/1' else
        'Selected sea state from the evolving second-cycle coast')
    r=np.loadtxt(path.parent/'runup.tbl')
    settings=record['settings']
    record['boundary_formulation']=dict(
        velocity_profile=('Discrete first-order Stokes profile for the uniform Keller-box layers'
            if settings['layers']<=4 else 'Continuum hyperbolic-cosine profile'),
        reading_rule='SWASH 12.01 changes the incident velocity reconstruction above four uniform Keller-box layers. A four-to-six-layer comparison therefore includes the boundary formulation.')
    np.testing.assert_allclose(r[[0,-1],0],[settings['spinup_s'],settings['seconds']],atol=.11,rtol=0)
    if np.any(np.diff(r[:,0])<=0) or not np.isfinite(r).all() or np.any(r[:,1]<-100):
        raise ValueError('Run-up time history is incomplete or invalid')
    np.testing.assert_allclose(record['statistics']['runup_time_p98_m'],np.percentile(r[:,1],98),atol=1e-10)
    reach=record['statistics']['breaking_x_range_m']
    if reach is not None:
        bed=np.loadtxt(path.parent/'statistics.tbl')
        shore=settings['geometry']['still_water_shore_m']
        record['derived_breaking']=dict(
            offshore_edge_depth_m=float(np.interp(reach[0],bed[:,0],bed[:,1])),
            offshore_edge_hs_m=float(np.interp(reach[0],bed[:,0],bed[:,2])),
            offshore_edge_distance_from_shore_m=reach[0]-shore,
            onshore_edge_distance_from_shore_m=reach[1]-shore,
            span_m=reach[1]-reach[0],
            definition='Seaward and landward edges of the native nodes flagged in at least 1% of the one-second samples; Hs is the full reported-window surface standard-deviation diagnostic.')
    else:
        record['derived_breaking']=None
    return record


def comparison(a,b,label):
    measures=('runup_time_p98_m','runup_peak_p98_m','runup_max_m','runup_mean_m','maximum_profile_hs_m')
    errors={key:abs(a['statistics'][key]-b['statistics'][key])/max(abs(b['statistics'][key]),.01) for key in measures}
    return dict(reference=a['name'],test=b['name'],label=label,relative_differences=errors,
                numerical_bounds=dict(runup_time_p98=.05,runup_peak_p98=.05,runup_max=.1),
                within_numerical_bounds=bool(errors['runup_time_p98_m']<.05 and errors['runup_peak_p98_m']<.05 and errors['runup_max_m']<.1))


def analyse(selection=None,output=OUTPUT):
    records={}
    for path in sorted(RUNS.glob('*/product.json')):
        case=verify(path)
        case.update(path=str(path.relative_to(ROOT)),sha256=sha256(path))
        records[case['name']]=case
    if not records:
        raise ValueError('No completed shore-profile cases')
    controls={}
    for name in ('controls.json','shore_controls.json'):
        path=RUNS/name
        controls[name]=dict(path=str(path.relative_to(ROOT)),sha256=sha256(path),record=json.loads(path.read_text()))
    for name,case in controls['controls.json']['record']['cases'].items():
        verify_run(RUNS/f'control_{name}_ilu',case['run'])
    shore=controls['shore_controls.json']['record']
    for spacing,case in shore['still_water'].items():
        verify_run(RUNS/f'rest_fourier_dx{spacing.replace(".","p")}',case['run'])
    verify_run(RUNS/'nonlinear_four_gravity',shore['nonlinear_gravity']['run'])
    verify_run(RUNS/'reflection_wall',shore['reflection']['run'])
    verify_run(RUNS/'irregular_input_flat',shore['irregular_input']['run'])
    for case in shore['irregular_input_refinements'].values():
        verify_run(RUNS/case['directory'],case['run'])
    pairs=[('calibration_west_rock_dx0p5_k4_ilu','calibration_west_rock_dx0p25_k4','Rock spacing 0.5 to 0.25 m'),
           ('calibration_west_rock_dx0p5_k4_ilu','calibration_west_rock_dx0p5_k6','Four to six layers, including the incident velocity formulation')]
    pairs.extend((f'calibration_west_slope{slope}_dx1_k4',f'calibration_west_slope{slope}_settled',
                  f'1:{slope} beach, 600 to 4200 s initialization with the same forcing phases') for slope in (20,50,100))
    pairs.append(('calibration_west_slope50_settled','calibration_west_slope50_dx0p5_settled',
                  'Settled 1:50 beach spacing 1 to 0.5 m'))
    pairs.append(('calibration_west_slope50_settled','calibration_west_slope50_later',
                  '1:50 beach, 4200 to 7800 s initialization with the same forcing phases'))
    pairs.append(('history_west_rock','history_west_rock_fine','Evolving-history peak, rock spacing 0.5 to 0.25 m'))
    coast_beach='history_west_slope50_dx1' if 'history_west_slope50_dx1' in records else 'history_west_slope50'
    pairs.append((coast_beach,'history_west_slope50_fine',
                  'Evolving-history peak, 1:50 beach spacing 1 to 0.5 m'))
    comparisons=[comparison(records[a],records[b],label) for a,b,label in pairs if a in records and b in records]
    sensitivities=[]
    a,b='history_east_bandlimited_k4','history_east_bandlimited_k6'
    if a in records and b in records:
        if records[a]['run']['identity']['inputs']['spectrum.dat']!=records[b]['run']['identity']['inputs']['spectrum.dat']:
            raise ValueError('The eastern layer comparison requires identical incoming spectrum bytes')
        band=comparison(records[a],records[b],
            'Eastern rock, four to six layers with the same restricted frequency band')
        band.update(incoming_spectrum_sha256=records[a]['run']['identity']['inputs']['spectrum.dat'],
            removed_variance_fraction=records[a]['boundary']['removed_variance_fraction'],
            reading_rule='An explicit common-band sensitivity includes the vertical discretization and the incident velocity formulation. The relaxed 5% variance guard is specific to these labelled cases; the main transfers retain the 1% guard.')
        comparisons.append(band)
    if 'history_east_face_rock' in records and b in records:
        baseline=records['history_east_face_rock']['statistics']
        values=records[b]['statistics']
        sensitivities.append(dict(reference='history_east_face_rock',test=b,
            parameter='Eastern six-layer case with its frequency ceiling reduced to the four-layer limit',
            changes_m={key:values[key]-baseline[key] for key in ('runup_time_p98_m','runup_max_m','runup_mean_m')},
            reading_rule='The frequency tail and finite realization change together; this complements the identical-spectrum four-to-six-layer comparison.'))
    if 'calibration_west_rock_buffer200' in records:
        baseline=records['calibration_west_rock_dx0p5_k4_ilu']['statistics']
        values=records['calibration_west_rock_buffer200']['statistics']
        sensitivities.append(dict(reference='calibration_west_rock_dx0p5_k4_ilu',
            test='calibration_west_rock_buffer200',parameter='Constant-depth wavemaker reach 100 to 200 m',
            changes_m={key:values[key]-baseline[key] for key in ('runup_time_p98_m','runup_max_m','runup_mean_m')},
            reading_rule='Changing the buffer also changes frequency-dependent travel phases within the finite reporting window. This is a boundary-placement sensitivity, with those phase effects retained.'))
    if 'history_west_rock' in records:
        baseline=records['history_west_rock']['statistics']
        for name,parameter in (('history_west_rock_seed2','Random-phase realization'),
                               ('history_west_rock_friction','Dimensionless bottom friction 0 to 0.005'),
                               ('history_west_rock_alpha0p5','Breaking onset coefficient 0.6 to 0.5')):
            if name in records:
                values=records[name]['statistics']
                sensitivities.append(dict(reference='history_west_rock',test=name,parameter=parameter,
                    changes_m={key:values[key]-baseline[key] for key in ('runup_time_p98_m','runup_max_m','runup_mean_m')},
                    reading_rule='A finite case comparison for the stated parameter; an ensemble probability distribution requires more realizations.'))
    build=json.loads((BUILD/'build.json').read_text())
    result=dict(schema='terluna.climate.shore-runup/1',
        evidence='Wave-resolving response of explicit rock and beach profiles to selected lunar coastal spectra.',
        reading_rule='Run-up statistics describe each finite imposed sea state and random-phase realization. The directional input is reduced to its normal incoming power. Rock information retains 118 m geographic spacing; beach slopes are explicit scenarios.',
        producer=dict(path='climate/waves/runup_checks.py',sha256=sha256(Path(__file__))),
        model=dict(name='SWASH',version=build['version'],source_url=build['source_url'],source_sha256=build['source_sha256'],
                   executable_sha256=build['executable_sha256'],source_modifications=build['source_modifications'],
                   build_record_path=str((BUILD/'build.json').relative_to(ROOT)),build_record_sha256=sha256(BUILD/'build.json')),
        geography=dict(path=str(PROFILES.relative_to(ROOT)),sha256=sha256(PROFILES),record=json.loads(PROFILES.read_text())),
        units=dict(runup='vertical metres above still water',height='m',period='s',power='W per metre of contour'),
        selected_cases=selection or [],cases=records,controls=controls,comparisons=comparisons,
        sensitivities=sensitivities,
        limitations=['The 1D reduction preserves the incoming direction projection at the source; alongshore currents and two-dimensional wave groups remain further calculations.',
                     'The irregular-spectrum flat-water controls retain amplitude and period errors; boundary reconstruction, propagation resolution and finite sampling accompany the shore-height estimates.',
                     'Interpolated flooded rock and smooth hypothetical beaches retain separate geographic evidence states.',
                     'Friction, breaking thresholds, phase realization and finite sampling remain sensitivities.',
                     'The spectral sea and the wave-resolving shore exchange forcing in one direction; reflected shore waves leave through the local boundary.',
                     'Overturning crests, entrained air, sediment evolution and lunar physical calibration remain open.'])
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(completed_cases=len(records),selected_cases=result['selected_cases'])),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--select',nargs='*')
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args()
    analyse(args.select,args.output)
