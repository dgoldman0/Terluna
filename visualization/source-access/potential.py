"""Render the resource and reclamation products without adding physical models."""
from pathlib import Path
import gzip,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
paths={'population':ROOT/'resources/source_access/results/population.json.gz','comparison':ROOT/'resources/source_access/results/potential.json.gz','reclamation':ROOT/'resources/source_access/inputs/reclamation.json.gz'}
p={k:json.loads(gzip.decompress(v.read_bytes())) for k,v in paths.items()}
out=ROOT/'research/runs/source_access/potential';out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(15,5.2));a,b,c=axes
pools=['narrow','wider','all_no_don_quixote'];labels=['Narrow orbital screen','Wider orbital screen','Broader NEAs\nwithout Don Quixote']
roundtrip={}
for i,pool in enumerate(pools):
 q=np.quantile(p['population']['scenarios']['central']['pools'][pool]['total']['water_equivalent_kg'],[.05,.5,.95])/1e12
 a.errorbar(q[1],i,xerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',capsize=6,color='#246b8e');roundtrip[pool]=q.tolist()
 text=' / '.join(f'{v:,.0f}' if v>=100 else f'{v:.2f}' if v<1 else f'{v:.1f}' for v in q)
 a.annotate(text+' Gt',(q[1],i),xytext=(0,12),textcoords='offset points',ha='center',fontsize=9)
a.set_xscale('log');a.set_yticks(range(3),labels);a.set_ylim(-.5,2.6);a.set_xlabel('Contained water equivalent (Gt)');a.set_title('Geological potential\nCentral-model p05 / median / p95');a.grid(axis='x',alpha=.2)
colors={'baseline':'#73655d','combined':'#178477','combined_without_residual_allowance':'#a35c8f'}
names={'baseline':'Existing recovery','combined':'Recovery + refurbishment','combined_without_residual_allowance':'Additional allowance removed'}
for variant in colors:
 case=next(x for x in p['reclamation']['cases'] if x['baseline_case']=='services:breakthrough:reference' and x['variant']==variant)
 rows=[r for r in case['annual'] if r['year']>=8]
 style='--' if variant=='combined_without_residual_allowance' else '-'
 b.plot([r['year'] for r in rows],[r['source_feed_kg']/1e12 for r in rows],label=names[variant],color=colors[variant],ls=style)
b.set_title('Fresh asteroid feed\nServices, 0.5% real market growth');b.set_xlabel('Year');b.set_ylabel('Gt per year');b.legend(fontsize=8,loc='upper left');b.grid(alpha=.2)
for pool,color,label in [('wider','#246b8e','Wider orbital screen'),('all_no_don_quixote','#cf7c35','Broader NEAs without Don Quixote')]:
 rows=[r for r in p['comparison']['comparisons'] if r['population_scenario']=='central' and r['pool']==pool and r['baseline_case']=='services:breakthrough:reference' and r['variant']=='combined' and r['horizon']==200]
 rows.sort(key=lambda r:r['available_share'])
 c.plot([r['available_share']*100 for r in rows],[r['joint_feed_and_water_probability']*100 for r in rows],marker='o',color=color,label=label)
c.set_xscale('log');c.set_xticks([1,5,10,25,50,100],[1,5,10,25,50,100]);c.set_xlabel('Available share of stock (%)');c.set_ylabel('Conditional simulation frequency (%)');c.set_ylim(-3,103);c.set_title('Stock covers year-200 demand\nServices with recovery + refurbishment');c.grid(alpha=.2);c.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18))
fig.suptitle('Material potential and fresh extraction in the enclosed phase',fontsize=16,y=1.02)
fig.text(.5,-.02,'Conditional models, not mine or delivery qualification. Water recovery: 70%. The 100% availability point is an inventory ceiling.\nThe additional-water-allowance removal is a maximum-credit sensitivity. 1 Gt = 1 trillion kg.',ha='center',fontsize=9)
fig.tight_layout();fig.savefig(out/'potential.png',dpi=170,bbox_inches='tight');fig.savefig(out/'potential.svg',bbox_inches='tight');plt.close(fig)
manifest={'inputs':{k:{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in paths.items()},'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'water_quantiles_Gt':roundtrip,'note':'All curves and quantiles render supplied domain products; no new source capacity or process model.'}
(out/'potential.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(out/'potential.png')
