"""Compare a propagating spectral pulse with its analytic lunar group delay."""
from __future__ import annotations

from datetime import timedelta
import json
from pathlib import Path

import numpy as np

from climate.waves.model import ROOT, ORIGIN, sha256
from climate.waves.pilot import array_text, execute
from climate.waves.shore import stationary_boundary
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

RUNS=ROOT/'research/runs/waves/shore_history/timing_controls'
WATER_DENSITY=json.loads((ROOT/'shared/scenarios/waves.json').read_text())['water_density_kg_m3']


def pulse_files(step=900,cells=40):
    f=np.geomspace(.03,.3e1,49)*MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    directions=np.arange(36)*10.+5
    index=30
    width=(f[index+1]-f[index-1])/2
    source=dict(locations=np.array([[0.,0.]]),frequency=f,direction=directions,
                available=np.ones((1,1),bool),spectra=np.zeros((1,1,len(f),len(directions))))
    pieces=[]
    for hour in np.arange(0,24.01,.25):
        source['spectra'][0,0,index,0]=(.5/4)**2*np.exp(-.5*((hour-6)/1.5)**2)/width/10
        text=stationary_boundary(source,0)
        header,body=text.split('-99\n',1)
        if not pieces:
            pieces.append(header.replace('SWAN 1\n','SWAN 1\nTIME\n1\n').replace('LONLAT','LOCATIONS')+'-99\n')
        pieces.append((ORIGIN+timedelta(hours=float(hour))).strftime('%Y%m%d.%H%M%S')+'\n'+body)
    start=ORIGIN.strftime('%Y%m%d.%H%M%S');end=(ORIGIN+timedelta(hours=24)).strftime('%Y%m%d.%H%M%S')
    text=f"""PROJECT 'Lunar pulse' '001'
MODE NONSTAT ONED
SET GRAV {MOON_SURFACE_GRAVITY:.12g} RHO {WATER_DENSITY:.12g}
OUTPUT OPTIONS TABLE 16
CGRID REG 0 0 0 10000 0 {cells} 0 CIRCLE 36 {f[0]:.12g} {f[-1]:.12g} 48
INPGRID BOTTOM 0 0 0 1 0 10000 1
READINP BOTTOM 1 'depth.bot' 1 0 FREE
BOUNDSPEC SIDE W CONSTANT FILE 'pulse.spc'
GEN3 KOMEN DRAG WU
OFF QUAD
OFF WCAP
BREAKING CONSTANT 1 .73
PROP BSBT
INIT ZERO
QUANTITY TSEC REF {start}
POINTS 'gauges' 0 0 2000 0 5000 0 8000 0
TABLE 'gauges' NOHEAD 'pulse.tbl' TSEC XP HS TM01 QB OUTPUT {start} 15 MIN
COMPUTE NONSTAT {start} {step} SEC {end}
STOP
"""
    return dict(INPUT=text,**{'depth.bot':'1000 1000\n','pulse.spc':''.join(pieces)}),float(f[index])


def run(executable):
    cases={}
    for step in (900,450,150):
        files,frequency=pulse_files(step)
        directory=RUNS/f'dt{step}'
        result=execute(executable,directory,files,['pulse.tbl'],timeout_s=300,threads=1)
        data=np.loadtxt(directory/'pulse.tbl').reshape(-1,4,5)
        t=data[:,0,0]/3600
        if not np.all(data[:,:,4]==0):
            raise ValueError('The deep-water timing control must stay free of breaking')
        energy=data[:,:,2]**2
        centroid=np.trapezoid(t[:,None]*energy,t,axis=0)/np.trapezoid(energy,t,axis=0)
        lag=centroid-centroid[0]
        speed=MOON_SURFACE_GRAVITY/(4*np.pi*frequency)*np.cos(np.deg2rad(5.))
        expected=data[0,:,1]/speed/3600
        difference=lag-expected
        half_rise=[]
        for i in range(4):
            e=energy[:,i];k=np.flatnonzero(e>=.5*e.max())[0]
            half_rise.append(float(t[k-1]+(t[k]-t[k-1])*(.5*e.max()-e[k-1])/(e[k]-e[k-1])))
        front_error=np.array(half_rise)-half_rise[0]-expected
        cases[str(step)]=dict(step_s=step,frequency_hz=frequency,group_speed_eastward_m_s=speed,
            x_m=data[0,:,1].tolist(),analytic_delay_hours=expected.tolist(),measured_delay_hours=lag.tolist(),
            delay_error_hours=difference.tolist(),maximum_delay_error_minutes=float(abs(difference).max()*60),
            half_peak_rise_error_minutes=(front_error*60).tolist(),
            pulse_height_ratio=(data[:,:,2].max(axis=0)/data[:,0,2].max()).tolist(),
            passed=bool(abs(difference).max()*60<5),run=result)
    record=dict(schema='terluna.climate.coastal-timing-controls/1',cases=cases,
        producer_sha256=sha256(Path(__file__)),
        evidence='A linear narrowband pulse in a flat lunar flume, checked against analytic group velocity.',
        reading_rule='Centroids measure travel delay; peak-height reduction records numerical spreading. Real coastal event onsets also depend on local growth and changing direction.')
    (RUNS/'controls.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:{q:v[q] for q in ('maximum_delay_error_minutes','pulse_height_ratio','passed')} for k,v in cases.items()}),flush=True)
    return record


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    run(parser.parse_args().executable)
