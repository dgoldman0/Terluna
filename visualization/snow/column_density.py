"""Display the recorded warm-column error comparisons; no growth physics."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
PRODUCT=ROOT/'atmosphere/snow/results/mcsnow_column_density.json'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def render():
    data=json.loads(PRODUCT.read_text())
    if data['schema']!='terluna.atmosphere.mcsnow-column-density/1':raise ValueError('Unsupported product')
    expected={**data['producer'],**data['inputs'],data['previous_product_path']:data['previous_product_sha256']}
    for name,sha in expected.items():
        if digest(ROOT/name)!=sha:raise ValueError('Source identity changed: '+name)
    selected=[f for f in data['searches'] if f['name'] in data['selected_searches']]
    curves=[(data['previous_warm'],'Previous warm fit','#888888')]
    labels={
        'global_20261024':('Density refit, 40 µm stage','#c47b1b'),
        'compact_local_20261024':('Density and compact-stage refit','#8553a3'),
        'seed_state_density_floor':('Seed constraint, gentler density change','#087f8c'),
        data['preferred_search']:('Early constraints and balanced join weight','#b53d59')}
    for f in selected:
        if f['name'] in labels:curves.append((f,*labels[f['name']]))
    fig,axes=plt.subplots(2,3,figsize=(11,8),sharex=True)
    plotted=[]
    for i,start in enumerate(('sphere','axes')):
        for j,quantity in enumerate(('Mass','Width','Length')):
            ax=axes[i,j];ax.axhspan(-10,10,color='#dcecdc',alpha=.75,zorder=0)
            ax.axhline(0,color='#444444',lw=.7)
            for f,label,color in curves:
                source=f.get('interval_assessment',f['assessment'])
                rows=sorted([r for r in source['rows'] if r['initialization']==start],key=lambda r:r['age_min'])
                ages=np.array([r['age_min'] for r in rows]);ratios=np.array([r['ratio'][j] for r in rows])
                line,=ax.plot(ages,100*(ratios-1),'.-',label=label,color=color)
                if not np.array_equal(line.get_xdata(),ages):raise ValueError('Age coordinates changed')
                if not np.allclose(line.get_ydata()/100+1,ratios,rtol=3e-15,atol=0):raise ValueError('Error conversion changed ratios')
                plotted.append(dict(fit=f['name'],initialization=start,quantity=quantity,
                                    ages_min=ages.tolist(),native_ratios=ratios.tolist()))
            ax.grid(alpha=.18);ax.set_title(quantity);ax.set_ylim(-30,35)
            if j==0:ax.set_ylabel(('4 µm spherical seed' if start=='sphere' else 'Three-minute restart')+'\nError (%)')
            if i==1:ax.set_xlabel('Age (minutes)')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.945),frameon=False,ncol=2)
    fig.suptitle('Warm column at −5.3 °C: evolving deposition-density tests',y=.99)
    cold=max(f['cold_control']['maximum_relative_difference'] for f in selected)
    fig.text(.5,.035,'Shading: ±10%. Final fitting ages: 3, 4, 5, just before 6, 10 and 20 min. The length reference changes abruptly at 6 min.',ha='center',fontsize=8)
    fig.text(.5,.012,f'Cold benchmark: maximum change across the tested states = {cold:.1g}. Separate temperature-specific controls; sampled errors against fitted laboratory histories.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.075,1,.84))
    out=HERE/'generated';out.mkdir(parents=True,exist_ok=True)
    identity=dict(product_sha256=digest(PRODUCT),renderer_sha256=digest(Path(__file__)),evidence=data['evidence'])
    png=out/'column_density.png';pdf=out/'column_density.pdf'
    fig.savefig(png,dpi=160,metadata={'Description':json.dumps(identity,sort_keys=True)})
    fig.savefig(pdf,metadata={'Title':'Warm-column density and compact-stage tests','Subject':json.dumps(identity,sort_keys=True),'CreationDate':None,'ModDate':None})
    plt.close(fig)
    manifest=dict(schema='terluna.visualization.snow-column-density/1',**identity,inputs=expected,
        plotted_values=plotted,reading_rule='Error = 100 * (model / laboratory fitted history - 1); lines join saved ages. Cold parity is reported separately.',
        outputs={p.name:digest(p) for p in (png,pdf)})
    (out/'column_density_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(curves=len(plotted),outputs=manifest['outputs']),indent=2))


if __name__=='__main__':render()
