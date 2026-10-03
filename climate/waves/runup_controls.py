"""Check still water and nonlinear gravity similarity on the rock profile."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from climate.waves.model import ROOT, sha256
from climate.waves.pilot import array_text, execute
from climate.waves.runup import PROFILES, RUNS, profile_files, summarize
from climate.waves.swash_build import build
from geography.shore_profiles import numerical_profile
from shared.constants import MOON_SURFACE_GRAVITY


def reflection_control(executable):
    """Resolve incident/reflected linear waves in front of a closed wall."""
    x=np.arange(0.,160.01,.5)
    files,outputs=profile_files(x,np.full_like(x,5.),'',layers=3,seconds=600,
        spinup=400,regular=(.02,10.),flat=True)
    files['INPUT']=files['INPUT'].replace('BOU SIDE E CCW BTYPE RADIATION\n','').replace('SPON E 40\n','')
    files['INPUT']=files['INPUT'].replace('SMOO 120 SEC','SMOO 20 SEC')
    directory=RUNS/'reflection_wall'
    record=execute(executable,directory,files,outputs,timeout_s=300,threads=1)
    data=np.loadtxt(directory/'gauges.tbl').reshape(-1,3,3)
    selected=data[:,0,0]>=400
    t=data[selected,0,0]
    basis=np.column_stack((np.cos(2*np.pi*t/10),np.sin(2*np.pi*t/10),np.ones(len(t))))
    coefficient=np.linalg.lstsq(basis,data[selected,:,2],rcond=None)[0]
    amplitude=coefficient[0]-1j*coefficient[1]
    k=brentq(lambda k:MOON_SURFACE_GRAVITY*k*np.tanh(k*5)-(2*np.pi/10)**2,1e-8,10)
    basis=np.column_stack((np.exp(-1j*k*data[0,:,1]),np.exp(1j*k*data[0,:,1])))
    incoming,reflected=np.linalg.lstsq(basis,amplitude,rcond=None)[0]
    residual=float(np.linalg.norm(basis@[incoming,reflected]-amplitude)/np.linalg.norm(amplitude))
    ratio=float(abs(reflected/incoming))
    return dict(run=record,incident_amplitude_m=float(abs(incoming)),reflected_amplitude_m=float(abs(reflected)),
                reflection_amplitude_ratio=ratio,relative_fit_residual=residual,
                passed=bool(abs(ratio-1)<.05 and residual<.03),
                reading_rule='A fully reflective vertical wall has unit amplitude reflection. Three gauges separate the two linear travelling components using the analytic lunar dispersion relation.')


def irregular_input_control(executable,*,dx=.5,layers=4):
    """Check the actual transferred coastal spectrum in a flat-water flume."""
    original=RUNS/'calibration_west_rock_dx0p5_k4_ilu'
    source=json.loads((original/'product.json').read_text())
    spectrum_path=original/'spectrum.dat'
    if sha256(spectrum_path)!=source['run']['identity']['inputs']['spectrum.dat']:
        raise ValueError('Irregular-input reference spectrum changed')
    spectrum=np.loadtxt(spectrum_path,skiprows=1)
    df=spectrum[1,0]-spectrum[0,0]
    variance=spectrum[:,1]*df
    expected_hs=float(4*np.sqrt(variance.sum()))
    expected_period=float(variance.sum()/np.sum(spectrum[:,0]*variance))
    x=np.arange(0.,400.01,dx)
    files,outputs=profile_files(x,np.full_like(x,source['boundary']['depth_m']),
        spectrum_path.read_text(),dx=dx,step=dx/20,layers=layers,seconds=4500,spinup=900,flat=True)
    files['INPUT']=files['INPUT'].replace('SPON E 40','SPON E 160')
    name=('irregular_input_flat' if dx==.5 and layers==4 else
          f'irregular_input_flat_dx{str(dx).replace(".","p")}_k{layers}')
    directory=RUNS/name
    record=execute(executable,directory,files,outputs,timeout_s=2400,threads=1)
    data=np.loadtxt(directory/'gauges.tbl').reshape(-1,3,3)
    data=data[(data[:,0,0]>=900)&(data[:,0,0]<4500)]
    expected_time=np.arange(900,4500,.1)
    # The printed clock has single-precision rounding at these elapsed times.
    np.testing.assert_allclose(data[:,0,0],expected_time,
                               atol=float(np.spacing(np.float32(4500))),rtol=0)
    eta=data[:,:,2]-data[:,:,2].mean(axis=0)
    heights=4*np.std(eta,axis=0)
    power=abs(np.fft.rfft(eta,axis=0))**2
    power[1:-1]*=2
    frequency=np.fft.rfftfreq(len(eta),.1)
    periods=power.sum(axis=0)/(frequency[:,None]*power).sum(axis=0)
    hs_error=float(np.max(abs(heights/expected_hs-1)))
    period_error=float(np.max(abs(periods/expected_period-1)))
    return dict(run=record,directory=name,dx_m=dx,layers=layers,
        source_case_sha256=sha256(original/'product.json'),
        expected_hs_m=expected_hs,expected_mean_period_s=expected_period,
        measured_hs_m=heights.tolist(),measured_mean_period_s=periods.tolist(),
        maximum_clock_rounding_s=float(np.max(abs(data[:,0,0]-expected_time))),
        maximum_breaking_flag=float(np.loadtxt(directory/'field.tbl',usecols=[5]).max()),
        maximum_hs_relative_difference=hs_error,maximum_period_relative_difference=period_error,
        passed=bool(hs_error<.05 and period_error<.05),
        reading_rule='The transferred irregular spectrum drives constant-depth water with an absorbing downstream reach. Three gauges cover one complete 3600 s phase period after 900 s initialization. Nonlinear evolution and residual reflection accompany the input-normalization check.')


def run():
    executable=build()
    profiles=json.loads(PROFILES.read_text())
    rest={}
    for dx in (.5,.25):
        x,depth,geometry=numerical_profile(profiles,'west_face',dx=dx)
        files,outputs=profile_files(x,depth,'',dx=dx,step=dx/20,seconds=120,spinup=60,regular=(0.,12.))
        files['INPUT']=files['INPUT'].replace('CON REG 0 12','CON FOURIER 0 0 0.523598775598 0')
        directory=RUNS/f'rest_fourier_dx{str(dx).replace(".","p")}'
        record=execute(executable,directory,files,outputs,timeout_s=120,threads=1)
        field=np.loadtxt(directory/'field.tbl').reshape(-1,len(x),6)
        # Cells safely inside the wet part have an exact zero free surface.
        wet=field[0,:,2]>.5
        residual=float(np.max(abs(field[:,wet,3])))
        r=np.loadtxt(directory/'runup.tbl')[:,1]
        rest[str(dx)]=dict(dx_m=dx,maximum_wet_surface_residual_m=residual,
                          still_water_runup_diagnostic_m=float(r.mean()),
                          diagnostic_range_m=float(np.ptp(r)),passed=residual<1e-5,run=record)
    original=RUNS/'calibration_west_rock_dx0p5_k4_ilu'
    source=json.loads((original/'product.json').read_text())
    files={}
    for filename,digest in source['run']['identity']['inputs'].items():
        if sha256(original/filename)!=digest:
            raise ValueError('Nonlinear gravity control source changed')
        files[filename]=(original/filename).read_text()
    settings=files['INPUT']
    changes={f'GRAV {MOON_SURFACE_GRAVITY:.12g}':f'GRAV {4*MOON_SURFACE_GRAVITY:.12g}',
             'SMOO 120 SEC':'SMOO 60 SEC','DUR 1800 SEC':'DUR 900 SEC',
             'OUTPUT 001000.000 1 SEC':'OUTPUT 000500.000 .5 SEC',
             'OUTPUT 000000.000 .1 SEC':'OUTPUT 000000.000 .05 SEC',
             'OUTPUT 001000.000 .1 SEC':'OUTPUT 000500.000 .05 SEC',
             'COMPUTE 000000.000 0.025 SEC 004000.000':'COMPUTE 000000.000 0.0125 SEC 002000.000'}
    for old,new in changes.items():
        if settings.count(old)!=1:
            raise ValueError(f'Gravity-control anchor differs: {old}')
        settings=settings.replace(old,new)
    spectrum=np.loadtxt(original/'spectrum.dat',skiprows=1)
    spectrum[:,0]*=2
    spectrum[:,1]/=2
    files['INPUT']=settings
    files['spectrum.dat']='SPEC1D\n'+array_text(spectrum)
    directory=RUNS/'nonlinear_four_gravity'
    outputs=['statistics.tbl','field.tbl','gauges.tbl','runup.tbl']
    record=execute(executable,directory,files,outputs,timeout_s=600,threads=1)
    a=np.loadtxt(original/'runup.tbl');b=np.loadtxt(directory/'runup.tbl')
    np.testing.assert_allclose(a[:,0],2*b[:,0],rtol=0,atol=.001)
    p98=np.percentile(a[:,1],98);other=np.percentile(b[:,1],98)
    field_a=np.loadtxt(original/'statistics.tbl');field_b=np.loadtxt(directory/'statistics.tbl')
    selected=field_a[:,2]>=.1
    hs_error=float(np.max(abs(field_a[selected,2]-field_b[selected,2])/field_a[selected,2]))
    similarity=dict(source_sha256=sha256(original/'product.json'),gravity_ratio=4,time_ratio=.5,
        maximum_pointwise_runup_difference_m=float(np.max(abs(a[:,1]-b[:,1]))),
        runup_p98_relative_difference=float(abs(p98-other)/p98),
        maximum_profile_hs_relative_difference=hs_error,
        passed=bool(abs(p98-other)/p98<.02 and hs_error<.02),run=record)
    result=dict(schema='terluna.climate.runup-controls/1',producer_sha256=sha256(Path(__file__)),
        evidence='Still-water balance, linear wall reflection, and nonlinear gravity/clock similarity including breaking and shoreline wetting.',
        reading_rule='The still-water RUNUP value measures the wet-depth and grid diagnostic offset. Gravity similarity tests the configured equations and executable; lunar calibration remains a separate question.',
        still_water=rest,nonlinear_gravity=similarity,reflection=reflection_control(executable),
        irregular_input=irregular_input_control(executable),
        irregular_input_refinements={
            'half_spacing':irregular_input_control(executable,dx=.25),
            'six_layers':irregular_input_control(executable,layers=6)})
    (RUNS/'shore_controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(rest={k:{p:v[p] for p in ('maximum_wet_surface_residual_m','still_water_runup_diagnostic_m')} for k,v in rest.items()},
                         similarity={k:v for k,v in similarity.items() if k!='run'})),flush=True)
    return result


if __name__=='__main__':
    run()
