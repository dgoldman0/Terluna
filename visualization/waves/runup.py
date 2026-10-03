"""Render actual SWASH shoreline snapshots and recorded run-up histories."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from climate.waves.model import ROOT, sha256

SOURCE=ROOT/'climate/waves/results/runup.json'
OUTPUT=Path(__file__).parent/'results'


def render(source=SOURCE,output=OUTPUT):
    record=json.loads(source.read_text())
    if record['schema']!='terluna.climate.shore-runup/1' or not record['selected_cases']:
        raise ValueError('Select completed shore cases for the figure')
    selected=[record['cases'][key] for key in record['selected_cases']]
    fig,axes=plt.subplots(2,len(selected),figsize=(4.2*len(selected),7.2),squeeze=False,layout='constrained')
    checks=[]
    for column,case in enumerate(selected):
        path=ROOT/case['path']
        if sha256(path)!=case['sha256']:
            raise ValueError('Shore-case record changed')
        for name in ('statistics.tbl','field.tbl','runup.tbl'):
            if sha256(path.parent/name)!=case['run']['output_sha256'][name]:
                raise ValueError(f'Shore display source changed: {name}')
        statistics=np.loadtxt(path.parent/'statistics.tbl')
        field=np.loadtxt(path.parent/'field.tbl').reshape(-1,len(statistics),6)
        runup=np.loadtxt(path.parent/'runup.tbl')
        geometry=case['settings']['geometry']
        x=statistics[:,0]-geometry['still_water_shore_m']
        peak_time=runup[np.argmax(runup[:,1]),0]
        k=np.argmin(abs(field[:,0,0]-peak_time))
        snapshot=field[k]
        wet=snapshot[:,4]>=case['settings']['wet_threshold_m']
        ax=axes[0,column]
        bottom,=ax.plot(x,-statistics[:,1],color='#66553a',lw=1.5)
        surface,=ax.plot(x,np.where(wet,snapshot[:,3],np.nan),color='#267da6',lw=1.3)
        ax.fill_between(x,-statistics[:,1],snapshot[:,3],where=wet,color='#8dc7df',alpha=.55)
        ax.fill_between(x,-statistics[:,1],-20,color='#c5b48c',alpha=.6)
        breaking=wet & (snapshot[:,5]==1)
        marks=ax.scatter(x[breaking],snapshot[breaking,3],s=9,color='#b95e29',zorder=4)
        ax.axhline(0,color='#808080',ls=':',lw=.8)
        ax.set_xlim(max(x[0],-200),max(60,float(x[wet].max())+10))
        ax.set_ylim(-8,max(4,float(np.ceil(max(case['statistics']['runup_max_m'],
                                            snapshot[wet,3].max())+.5))))
        ax.set_xlabel('Distance from still-water shoreline (m)')
        ax.set_ylabel('Elevation above still water (m)')
        title='Interpolated rock' if geometry['kind']=='rock' else f"Assumed 1:{geometry['slope_denominator']} beach"
        station={'west_face':'Western face','east_face':'Eastern face','southern_shore':'Southern shore'}[geometry['station']]
        ax.set_title(f'{station}\n{title}')
        boundary=case['boundary']
        ax.text(.03,.04,f"Incoming Hs {boundary['incoming_hs_m']:.2f} m; mean period {boundary['mean_period_s']:.1f} s\n"
                f"Surface at {snapshot[0,0]:.0f} s",transform=ax.transAxes,fontsize=9)
        ax=axes[1,column]
        line,=ax.plot((runup[:,0]-runup[0,0])/60,runup[:,1],color='#267da6',lw=.8)
        p98=case['statistics']['runup_time_p98_m']
        ax.axhline(p98,color='#9c4e28',ls='--',lw=1)
        ax.set_xlabel('Minutes in the reporting window');ax.set_ylabel('Vertical run-up (m)')
        ax.set_title(f"98th time percentile {p98:.2f} m · maximum {case['statistics']['runup_max_m']:.2f} m",fontsize=10)
        ax.grid(alpha=.2)
        np.testing.assert_array_equal(bottom.get_ydata(),-statistics[:,1])
        np.testing.assert_allclose(surface.get_ydata(),np.where(wet,snapshot[:,3],np.nan),rtol=0,atol=0,equal_nan=True)
        np.testing.assert_array_equal(line.get_ydata(),runup[:,1])
        np.testing.assert_array_equal(marks.get_offsets(),np.column_stack((x[breaking],snapshot[breaking,3])))
        checks.append(dict(case=case['name'],snapshot_time_s=float(snapshot[0,0]),runup_peak_time_s=float(peak_time),exact_array_checks=4))
    fig.suptitle('Lunar breaking and run-up across explicit shore profiles',fontsize=16)
    fig.supxlabel('Orange points: computed breaking · SWASH at lunar gravity · 1D incoming-power reduction\n'
                  'Vertical scale expanded · rock information at 118 m · prescribed beach slopes · finite realizations',fontsize=9)
    output.mkdir(parents=True,exist_ok=True)
    path=output/'runup.png'
    fig.savefig(path,dpi=160);plt.close(fig)
    sidecar=dict(schema='terluna.visualization.shore-runup/1',source_sha256=sha256(source),image_sha256=sha256(path),
                 renderer_sha256=sha256(Path(__file__)),cases=checks,
                 reading_rule='Each upper panel shows the nearest one-second free-surface output to the largest recorded run-up. Lower panels show the complete reported waterline record.')
    path.with_suffix('.json').write_text(json.dumps(sidecar,indent=2)+'\n')
    print(path)


if __name__=='__main__':
    render()
