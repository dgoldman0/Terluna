"""Render the frozen state-growth product; no domain model is imported."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PRODUCT = ROOT/'atmosphere/snow/results/mcsnow_state_growth.json'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def render():
    d = json.loads(PRODUCT.read_text())
    if d['schema'] != 'terluna.atmosphere.mcsnow-state-growth/1':
        raise ValueError('Unsupported state-growth product')
    if not d.get('evidence') or not d.get('reading_rule'):
        raise ValueError('Missing evidence or reading rule')
    expected = [*d['producer'].items(), *d['inputs'].items(),
                ('atmosphere/snow/mcsnow_state_sources.json',d['source_record_sha256']),
                (d['previous_product_path'],d['previous_product_sha256'])]
    for name, sha in expected:
        if digest(ROOT/name)!=sha: raise ValueError('Input identity mismatch: '+name)
    selected=d['selected_family']
    fig, axes = plt.subplots(4,3,figsize=(11,11),sharex=True)
    labels = ['Mass (mg)','Width (mm)','Polar dimension (mm)']
    colors = {'original':'#999999','constant':'#d18a24',selected:'#087f8c'}
    names = {'original':'Original improved component','constant':'Best constant adjustment',
             selected:'Shared state-dependent adjustment'}
    plotted=[]
    def draw(ax, ages, values, scale, label, color, style, identity):
        y=np.asarray(values)*scale
        if np.any(~np.isfinite(y)) or np.any(y<=0): raise ValueError('Invalid plotted values')
        line,=ax.plot(ages,y,style,color=color,label=label)
        if not np.array_equal(line.get_xdata(),ages):
            raise ValueError('Plot ages changed')
        if not np.allclose(line.get_ydata()/scale,values,rtol=3e-15,atol=0):
            raise ValueError('Plot unit conversion failed')
        plotted.append(dict(identity=identity,ages_min=list(ages),native_values=list(values),scale=scale))
    for i,(t,init) in enumerate([(t,init) for t in d['training']['temperatures_C'] for init in ('sphere','axes')]):
        for j,scale in enumerate((1e6,1e3,1e3)):
            ax=axes[i,j]
            for name,color in colors.items():
                rows=[r for r in d['assessments'][name]['rows'] if r['temperature_C']==t and r['initialization']==init]
                draw(ax,[r['age_min'] for r in rows],[r['values'][j] for r in rows],scale,
                     names[name],color,'.-',[t,init,j,name])
            chosen_name=d['selected_temperature_searches'].get(str(t))
            if chosen_name:
                chosen=next(f for f in d['temperature_searches'] if f['name']==chosen_name)
                rows=[r for r in chosen['assessment']['rows'] if r['initialization']==init]
                draw(ax,[r['age_min'] for r in rows],[r['values'][j] for r in rows],scale,
                     'Temperature-specific adjustment','#7656a3','.-',[t,init,j,chosen_name])
            segments=d['laboratory_plot_segments'][str(t)][d['quantities'][j]]
            for k,segment in enumerate(segments):
                draw(ax,segment['age_min'],segment['values'],scale,
                     'Laboratory fitted history' if k==0 else None,'#222222','--',[t,init,j,'laboratory',k])
            ax.set_yscale('log');ax.grid(alpha=.18);ax.set_title(labels[j] if i==0 else '')
            if j==0: ax.set_ylabel(f'{t:g} °C\n'+('4 µm seed' if init=='sphere' else '3-minute restart'))
            if i==3: ax.set_xlabel('Age (minutes)')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,.97))
    fig.suptitle('Size-dependent growth: Earth diagnostic comparisons',y=.995)
    fig.text(.5,.012,'Fit ages: 5, 10, 20 min. Other measurements excluded from fitting. Laboratory curves are fitted histories.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.03,1,.90))
    out=HERE/'generated';out.mkdir(parents=True,exist_ok=True)
    identity=dict(product_sha256=digest(PRODUCT),renderer_sha256=digest(Path(__file__)),evidence=d['evidence'])
    png=out/'state_growth.png';pdf=out/'state_growth.pdf'
    fig.savefig(png,dpi=160,metadata={'Description':json.dumps(identity,sort_keys=True)})
    fig.savefig(pdf,metadata={'Title':'State-dependent Earth crystal growth','Subject':json.dumps(identity,sort_keys=True),
                             'CreationDate':None,'ModDate':None})
    plt.close(fig)
    manifest=dict(schema='terluna.visualization.snow-state-growth/1',**identity,
                  inputs=dict(expected),reading_rule=d['reading_rule'],plotted_values=plotted,
                  outputs={p.name:digest(p) for p in (png,pdf)})
    (out/'state_growth_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(curves=len(plotted),outputs=manifest['outputs']),indent=2))


if __name__=='__main__': render()
