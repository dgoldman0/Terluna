"""Draw the proposed section and calculated demand of the conditional buildout."""
from pathlib import Path
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    result=json.loads((HERE/'results/buildout.json').read_text())
    spec=json.loads((HERE/'buildout_inputs.json').read_text())
    c=result['cases']['reference'];g=spec['geometry'];r=g['pressure_hall_diameter_m']/2
    fig,(section,power)=plt.subplots(1,2,figsize=(14,5.7),gridspec_kw={'width_ratios':[1.4,1]})
    fig.patch.set_facecolor('#f9f8f3')
    for ax in (section,power):ax.set_facecolor('#f9f8f3')
    angle=np.linspace(0,math.pi,300)
    section.plot(g['cave_width_m']/2*np.cos(angle),g['cave_height_m']*np.sin(angle),color='#6c6254',lw=3,label='Assumed cave boundary')
    section.plot([-150,150],[0,0],color='#6c6254',lw=3)
    for lane in range(g['lanes']):
        x=(lane-(g['lanes']-1)/2)*(2*r+g['pressure_lane_clearance_m'])
        section.add_patch(Circle((x,r),r,fc='#d4e8eb',ec='#1d6371',lw=1.6))
        for z in g['deck_height_fractions']:
            half=r*math.sqrt(1-z*z)
            section.plot([x-half,x+half],[r*(1+z)]*2,color='#4c8066',lw=1.3)
    section.set(xlim=(-160,160),ylim=(-7,115),xlabel='Horizontal distance from cave centre (m)',ylabel='Height above cave floor (m)')
    section.set_aspect('equal')
    section.text(0,76,'300 m wide · 100 m high',ha='center',color='#554d43',fontsize=11)
    section.set_title('Four 60 m pressure halls with five decks',loc='left',fontweight='bold',pad=16)
    section.spines[['top','right']].set_visible(False)
    fig.text(.06,.13,'Proposed section; no surveyed city-sized cavity.\nEach city requires 11.2 km of equivalent corridor.\nThe shells close pressure loads; rock stability needs its own survey.',fontsize=10,linespacing=1.6)
    loads=c['utilities']['mean_load_components_W']
    food,basic,cooling,garden,industry=[loads[k]/1e9 for k in ('food','basic_services','cooling','civic_gardens','everyday_industry')]
    labels=['Food production','Basic urban services','Cooling electricity','Civic garden lighting','Everyday industry']
    values=[food,basic,cooling,garden,industry]
    bars=power.barh(labels,values,color=['#4c8066','#1d6371','#5c8d9e','#91a96e','#b58a56'],height=.63)
    power.invert_yaxis();power.set_xlim(0,4.3)
    power.bar_label(bars,labels=[f'{v:.2f}' for v in values],padding=6,fontsize=10)
    power.set_xlabel('Mean electrical demand (GW)')
    power.set_title('6.30 GW for all 250,000 residents',loc='left',fontweight='bold',pad=16)
    power.spines[['top','right','left']].set_visible(False)
    power.tick_params(axis='y',length=0)
    power.grid(axis='x',alpha=.15);power.set_axisbelow(True)
    fig.subplots_adjust(left=.06,right=.98,top=.82,bottom=.27,wspace=.48)
    fig.suptitle('Conditional lava-tube city engineering',x=.06,ha='left',fontsize=17,fontweight='bold')
    fig.text(.06,.03,'Reference design assumptions and calculated loads · Four cities plus 100 towns · November 2026',fontsize=10,color='#554d43')
    target=HERE/'results/buildout_engineering.png';fig.savefig(target,dpi=170,facecolor=fig.get_facecolor());plt.close(fig)
    print(target)


if __name__=='__main__':main()
