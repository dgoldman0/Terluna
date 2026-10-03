"""Render the recorded second-cycle coastal time histories and exposure maps."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from climate.waves.cycle import window
from climate.waves.model import ROOT, sha256

SOURCE=ROOT/'climate/waves/results/coastal_history.json'
OUTPUT=Path(__file__).parent/'results'
COLORS={'west_face':'#205b9b','east_face':'#dc7b32','southern_shore':'#258879','offshore':'#7a7b83'}
LABELS={'west_face':'Western face','east_face':'Eastern face','southern_shore':'Southern shore','offshore':'Local offshore reference'}


def render(source=SOURCE,output=OUTPUT):
    record=json.loads(source.read_text())
    if record['schema']!='terluna.climate.coastal-history/1':
        raise ValueError('Unsupported coastal history')
    arrays=ROOT/record['arrays']['path']
    if sha256(arrays)!=record['arrays']['sha256']:
        raise ValueError('Coastal history arrays changed')
    with np.load(arrays,allow_pickle=False) as data:
        a={key:data[key].copy() for key in data.files}
    start,end=record['clock']['report_hours']
    t,h=window(a['time_hours'],a['stations'],start,end)
    pt,power=window(a['spectrum_time_hours'],a['power'],start,end)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,1,figsize=(12,7),sharex=True,layout='constrained')
    lines=[]
    for i,name in enumerate(record['stations']):
        line,=axes[0].plot((t-start)/24,h[:,i,3],color=COLORS[name],label=LABELS[name],lw=1.5)
        np.testing.assert_array_equal(line.get_ydata(),h[:,i,3]);lines.append(line)
        if name!='offshore':
            line,=axes[1].plot((pt-start)/24,power[:,i,0],color=COLORS[name],lw=1.5)
            np.testing.assert_array_equal(line.get_ydata(),power[:,i,0]);lines.append(line)
    axes[0].axhline(1,color='#999999',ls=':',lw=1)
    axes[0].set_ylabel('Significant wave height (m)')
    axes[0].legend(ncol=4,loc='upper right',frameon=False)
    axes[1].set_ylabel('Power travelling toward shore (W/m)')
    axes[1].set_xlabel('Earth days into the second lunar solar cycle')
    axes[1].set_xlim(0,(end-start)/24)
    for ax in axes:
        ax.grid(alpha=.2)
        ax.set_ylim(bottom=0)
    fig.suptitle('Eastern Smythii: waves arriving through a full lunar cycle',fontsize=16)
    fig.supxlabel('474 m wave grid on 118 m terrain · full first cycle of spin-up · one simulated reporting cycle',fontsize=9)
    output.mkdir(parents=True,exist_ok=True)
    path=output/'coastal_history.png'
    fig.savefig(path,dpi=170);plt.close(fig)
    sidecar=dict(schema='terluna.visualization.coastal-history/1',source_path=str(source.relative_to(ROOT)),
        source_sha256=sha256(source),arrays_sha256=sha256(arrays),image_sha256=sha256(path),
        renderer_sha256=sha256(Path(__file__)),exact_display_array_checks=len(lines),
        reading_rule=record['reading_rule'])
    path.with_suffix('.json').write_text(json.dumps(sidecar,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,5.5),layout='constrained')
    for ax,key,title in zip(axes,('map_mean_hs_m','map_max_hs_m'),('Mean significant height','Largest six-hourly significant height')):
        values=np.ma.array(a[key],mask=~a['wet'])
        mesh=ax.pcolormesh(a['longitude_deg'],a['latitude_deg'],values,cmap='viridis',vmin=0,shading='nearest')
        np.testing.assert_allclose(mesh.get_array().compressed(),values.compressed(),rtol=0,atol=0)
        ax.contour(a['longitude_deg'],a['latitude_deg'],a['wet'],levels=[.5],colors='#343434',linewidths=.7)
        for name,station in record['stations'].items():
            xy=station['longitude_latitude'];ax.scatter(*xy,s=24,c=COLORS[name],edgecolors='white',linewidths=.7)
            label={'west_face':'West','east_face':'East','southern_shore':'South','offshore':'Offshore'}[name]
            ax.annotate(label,xy,xytext=(-5,5) if name=='west_face' else (5,5),
                        textcoords='offset points',ha='right' if name=='west_face' else 'left',
                        fontsize=8,color='#202020',bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1))
        ax.set_title(title);ax.set_xlabel('Longitude (°E)');ax.set_ylabel('Latitude (°N)');ax.set_aspect('equal')
        fig.colorbar(mesh,ax=ax,label='Hs (m)',shrink=.75)
    fig.suptitle('Coastal exposure during the second cycle',fontsize=16)
    fig.supxlabel('Six-hourly maps with interpolated cycle endpoints · flooded-rock geography at the 28% water scenario',fontsize=9)
    path=output/'coastal_history_maps.png'
    fig.savefig(path,dpi=170);plt.close(fig)
    sidecar.update(image_sha256=sha256(path),exact_display_array_checks=2,reading_rule=record['maps']['reading_rule'])
    path.with_suffix('.json').write_text(json.dumps(sidecar,indent=2)+'\n')
    print(output/'coastal_history.png')


if __name__=='__main__':
    render()
