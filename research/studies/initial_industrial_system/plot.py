"""Reference capacity and capital figure from the committed calculation."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results/industrial_system.json').read_text());c=r['cases']['reference'];b=c['cost'];f=c['peak_flows']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,(a,z)=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[1.1,1]})
labels=['Asteroid mines and fleet','Lunar power constellation','Orbital depots and landers','Lunar factories and machinery','Earth–depot vehicles']
keys=['source','power','landing','lunar','earth_transfer'];values=[b['groups_USD2041'][k]/1e12 for k in keys]
a.barh(labels,values,color='#366d83');a.invert_yaxis();a.set_xlim(0,max(values)*1.26);a.set_xlabel('Trillion 2041 dollars, before programme allowances')
for i,v in enumerate(values):a.text(v+.10,i,f'{v:.2f}',va='center')
a.set_title('Industrial hardware and delivery from Earth',loc='left',pad=14)
flowkeys=['source_ore_kg_year','source_water_kg_year','nea_payload_kg_year','lunar_propellant_kg_year','LEO_kg_year','moon_payload_kg_year','moon_water_kg_year']
flowlabels=['Source material processed','Water extracted at sources','Cargo arriving at lunar depot','Lunar landing propellant','Earth payload and fuel to orbit','Cargo delivered to surface','Surface water delivery']
values=[f[k]/1e9 for k in flowkeys];z.barh(flowlabels,values,color='#60825d');z.set_xscale('log');z.set_xlim(.6,300);z.invert_yaxis();z.set_xlabel('Million tonnes per year — logarithmic scale')
for i,v in enumerate(values):z.text(v*1.08,i,f'{v:.2f}',va='center')
z.set_title('Rated production and transport flows',loc='left',pad=14)
fig.suptitle('Initial industrial supply for 250,000 lunar residents',fontsize=17,x=.02,ha='left',y=.99)
fig.text(.02,.03,'Conditional chemical-transport reference. Flow rows include connected streams and must not be added.\nCapital bars exclude research, preservation, starter water and contingency; total establishment is $16.2 trillion.',fontsize=10)
fig.subplots_adjust(left=.20,right=.96,top=.86,bottom=.20,wspace=.95)
fig.savefig(HERE/'results/industrial_system.png',dpi=160,facecolor='white')
