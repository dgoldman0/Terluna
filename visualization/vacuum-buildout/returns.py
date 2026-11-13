"""Plot the frozen year-100/year-200 returns and completion-funding requirement."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'research/studies/vacuum_buildout/results/returns.json'
OUT=ROOT/'research/runs/vacuum_buildout';OUT.mkdir(parents=True,exist_ok=True)
data=json.loads(SOURCE.read_text())
rows=[r for r in data['horizons'] if r['mode']=='latest_U_completion' and r['market']=='breakthrough']
names=list(dict.fromkeys(r['strategy'] for r in rows))
labels=['Raw-price','Commodity','Specialty','Premium','Propellant','Computing','Services','Diversified','Prior mix']
lookup={(r['strategy'],r['horizon']):r for r in rows}
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained');y=np.arange(len(names))
for year,offset,color in [(100,-.18,'#2586b0'),(200,.18,'#b67027')]:
 axes[0,0].barh(y+offset,[lookup[n,year]['discounted']['0.03']['sponsor_cash_NPV']/1e12 for n in names],height=.35,label=f'Year {year}',color=color)
axes[0,0].axvline(0,color='black',linewidth=.8);axes[0,0].set_yticks(y,labels);axes[0,0].invert_yaxis();axes[0,0].legend()
axes[0,0].set(title='Early receipts can support positive discounted returns',xlabel='Sponsor cash NPV at 3% (trillion USD2041)')
axes[0,1].barh(y,[lookup[n,200]['additional_external_equity']/1e12 for n in names],color='#b67027');axes[0,1].set_yticks(y,labels);axes[0,1].invert_yaxis()
axes[0,1].set(title='Completing the full build needs substantial funding',xlabel='Additional equity through year 200 (trillion USD2041)')
axes[1,0].barh(y,[lookup[n,200]['annual_operating_available']/1e12 for n in names],color='#a94a43');axes[1,0].set_yticks(y,labels);axes[1,0].invert_yaxis()
axes[1,0].set(title='The completed system has operating deficits',xlabel='Year-200 operating cash (trillion USD2041/year)')
r=lookup['propellant_logistics',200];rates=[0.,.01,.03,.07]
axes[1,1].plot([v*100 for v in rates],[r['discounted'][str(v)]['sponsor_cash_NPV']/1e12 for v in rates],marker='o',color='#175f8a')
axes[1,1].axhline(0,color='black',linewidth=.8)
axes[1,1].set(title='Propellant/logistics: sensitivity to time valuation',xlabel='Real discount rate (%)',ylabel='Year-200 sponsor cash NPV (trillion USD2041)')
for ax in axes.flat:ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
fig.suptitle('Strong worldwide demand, reference mining, latest U financing — completion depends on outside equity',fontsize=12)
for ext in ('png','svg'):fig.savefig(OUT/f'returns.{ext}',dpi=150)
(OUT/'returns-figure-manifest.json').write_text(json.dumps({'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print(OUT/'returns.png')
