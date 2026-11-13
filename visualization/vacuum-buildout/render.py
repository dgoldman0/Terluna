"""Render the frozen year-200 enclosed-build and mining study."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'research/studies/vacuum_buildout/results/buildout.json'
OUT=ROOT/'research/runs/vacuum_buildout'
data=json.loads(SOURCE.read_text());OUT.mkdir(parents=True,exist_ok=True)
rows=[c for c in data['cases'] if c['scenario']=='base']
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
r=rows[0]['annual'];years=[x['year'] for x in r]
axes[0,0].plot(years,[x['population']/1e6 for x in r],color='#175f8a')
axes[0,0].set(ylabel='Occupied capacity (million people)',title='Planned enclosed capacity reaches 6.47 million')
axes[0,1].plot(years,[x['construction_people']/1000 for x in r],color='#175f8a')
axes[0,1].set(ylabel='Construction orders (thousand people/year)',title='Two-year commissioning; growth then taper')
for c in rows:
    r=c['annual'];label=c['strategy'].replace('_',' ')
    axes[1,0].plot(years,[x['cumulative_source_feed_kg']/1e12 for x in r],label=label)
    axes[1,1].plot(years,[x['max_capacity_use'] for x in r],label=label)
axes[1,0].axhline(data['inventory_comparators']['narrow_catalogue_mass_kg']/1e12,color='black',linestyle='--',label='Narrow catalogue mass (not reserves)')
axes[1,0].set(ylabel='Cumulative source feed (billion tonnes)',title='More source qualification is required')
axes[1,1].axhline(1,color='black',linestyle='--')
axes[1,1].set(ylabel='Largest engineering capacity fraction',title='Expanded machinery, base feed assumptions')
for ax in axes.flat:
    ax.set_xlabel('Elapsed programme year');ax.grid(alpha=.2)
    ax.axvline(100,color='grey',alpha=.4);ax.axvline(200,color='grey',alpha=.4)
handles,labels=axes[1,0].get_legend_handles_labels()
fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
fig.suptitle('Conditional year-200 vacuum build-out — no qualified reserve or financing claim',fontsize=13)
for ext in ('png','svg'):fig.savefig(OUT/f'buildout.{ext}',dpi=150)
(OUT/'figure-manifest.json').write_text(json.dumps({'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print(OUT/'buildout.png')
