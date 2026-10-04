"""Display saved local acquisition and orbital-retiming results."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
PRODUCTS=ROOT/'research/studies/solar_shield_array/results'
OUTPUT=Path(__file__).with_name('results')/'active.png'


def main():
    names=['active_acquisition.json','active_retime.json','active_validation.json']
    p={n:json.loads((PRODUCTS/n).read_text()) for n in names}
    acq=p[names[0]]; cases=p[names[1]]['cases']; validation=p[names[2]]
    teal,orange,red='#007f79','#ad7900','#ad4150'
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titleweight':'bold'})
    fig,axes=plt.subplots(2,2,figsize=(12.7,9.7))
    fig.subplots_adjust(top=.84,bottom=.16,left=.09,right=.96,hspace=.46,wspace=.30)
    fig.suptitle('Active assistance with separate 10 km tiles',x=.09,y=.975,ha='left',fontsize=21,fontweight='bold')
    fig.text(.09,.932,'50 g/m² · DE440 ephemeris · finite-Sun shadows · explicit electric thrust and propellant',fontsize=11)
    fig.text(.09,.900,'Local gap closure and individual orbital transfers; a complete collision-safe lunar fleet remains unresolved.',fontsize=10,color='#555555')

    ax=axes[0,0]
    coarse=acq['coverage'];t=[c['minutes'] for c in coarse]
    ax.plot(t,[100*c['ray_coverage'] for c in coarse],'-o',ms=3,color=teal,label='Ray interception (32 Sun directions)')
    ax.plot(t,[100*c['all_sampled_sun_coverage'] for c in coarse],'-o',ms=3,color=orange,label='All 32 sampled Sun directions')
    fine=validation['acquisition_128_sun_coverage']
    ax.scatter([c['minutes'] for c in fine],[100*c['all_sampled_sun_coverage'] for c in fine],
               color=red,marker='x',s=45,zorder=4,label='128-direction checks')
    ax.set(xlabel='Minutes after acquisition starts',ylabel='Coverage (%)',xlim=(0,180),ylim=(0,105),
           title='625 tiles close a local 100 km-wide window')
    ax.legend(fontsize=8,loc='lower right',frameon=False);ax.grid(axis='y',alpha=.2)

    ax=axes[0,1]
    ax.semilogy(acq['errors_times_minutes'],np.maximum(acq['maximum_tile_error_m'],1e-5),color=teal,lw=1.5)
    ax.axhline(50,color=red,ls='--',lw=1,label='50 m tracking allowance')
    ax.set(xlabel='Minutes after acquisition starts',ylabel='Largest tile position error (m)',xlim=(0,180),
           ylim=(1e-5,1000),title='Bounded feedback brings the formation together')
    ax.legend(fontsize=8,frameon=False);ax.grid(axis='y',alpha=.2)

    ordered=[cases[k] for k in sorted(cases,key=float)];labels=[f"{c['days']:g} days" for c in ordered]
    ax=axes[1,0]
    bars=ax.bar(labels,[c['delta_v_m_s'] for c in ordered],color=[teal,orange],width=.5)
    ax.set(ylabel='Total integrated correction (m/s)',ylim=(0,8.5),title='About 1,400 km of orbital phase repositioning')
    for bar,c in zip(bars,ordered):
        ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.15,f"{c['delta_v_m_s']:.2f} m/s",ha='center',fontsize=10)
        ax.text(bar.get_x()+bar.get_width()/2,.35,f"{c['propellant_kg']:,.0f} kg propellant",ha='center',color='white',fontsize=8)
    ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)

    ax=axes[1,1];x=np.arange(len(ordered));w=.32
    ax.bar(x-w/2,[c['mean_electric_power_over_maneuver_W']/1e6 for c in ordered],width=w,color=teal,label='Mean over maneuver')
    ax.bar(x+w/2,[c['peak_electric_power_W']/1e6 for c in ordered],width=w,color=orange,label='Peak during thrust arcs')
    ax.set(xticks=x,xticklabels=labels,ylabel='Electrical propulsion power per tile (MW)',
           title='Four smooth six-hour thrust arcs')
    ax.legend(fontsize=8,frameon=False);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)

    fig.text(.09,.093,'Transfers use one 5-million-kg tile each. Independent arrival errors: 0.34 m (3 days) and 0.02 m (7 days).',fontsize=9)
    fig.text(.09,.068,'The 7-day solution uses more total energy: 0.286 GWh versus 0.229 GWh. These are feasible transfers, with no energy-optimality claim.',fontsize=9)
    fig.text(.09,.043,'Assumed exhaust: 30 km/s, 70% efficiency, 45° cant. Collector/habitat loads, exhaust safety and hardware mass feedback remain open.',fontsize=9)
    OUTPUT.parent.mkdir(exist_ok=True,parents=True);fig.savefig(OUTPUT,dpi=170,facecolor='white');plt.close(fig)
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    OUTPUT.with_suffix('.json').write_text(json.dumps(dict(products={str((PRODUCTS/n).relative_to(ROOT)):digest(PRODUCTS/n) for n in names},
        renderer_sha256=digest(Path(__file__)),image_sha256=digest(OUTPUT)),indent=2)+'\n')
    print(OUTPUT)


if __name__=='__main__':main()
