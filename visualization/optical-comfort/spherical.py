"""Plot the consistent spherical sky and its optical-comfort follow-up."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT/'research/runs/optical_comfort/figures'


def main():
    files = [ROOT/'illumination/sky/results/solved_sky.json',
             ROOT/'research/studies/optical_comfort/results/spherical_sky.json',
             ROOT/'research/studies/optical_comfort/results/native_sky.json',
             ROOT/'illumination/surface_light/results/surface_light.json']
    sky,new,old,surface = [json.loads(p.read_text()) for p in files]
    if sky['schema']!='terluna.illumination.solved-spherical-sky/1' or new['schema']!='terluna.research.optical-comfort-spherical-sky/1':
        raise ValueError('Unsupported spherical sky products')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(2,2,figsize=(12,8.5))
    colors = {'moon_1.2atm':'#4269a9','earth_control':'#309b86'}
    for name,label,case in [('moon_1.2atm','Shielded Moon','moon_1.2atm/design'),
                            ('earth_control','Earth control','earth_control/unfiltered')]:
        rows = sky['worlds'][name]['samples']
        day = [r for r in rows if r['sun_deg']>=0]
        dusk = [r for r in rows if -24<=r['sun_deg']<=5]
        axes[0,0].plot([r['sun_deg'] for r in day],[r['total_horizontal_lux']/1000 for r in day],
                       color=colors[name],label=label)
        earlier = sorted(surface['cases'][case]['summaries'].items(),key=lambda x:float(x[0]))
        axes[0,0].plot([float(k) for k,_ in earlier],[v['0.1']['illuminance_lux']/1000 for _,v in earlier],
                       color=colors[name],ls='--',alpha=.65)
        axes[0,1].semilogy([r['sun_deg'] for r in dusk],[r['total_horizontal_lux'] for r in dusk],
                           color=colors[name],label=label)
    axes[0,0].set(xlabel='Sun elevation (degrees)',ylabel='Horizontal illuminance (klux)',
                   title='Daylight · solid: spherical; dashed: two-stream')
    axes[0,0].legend(frameon=False)
    axes[0,1].set(xlabel='Sun elevation (degrees)',ylabel='Horizontal illuminance (lux)',
                   title='Clear twilight on the solved columns')
    axes[0,1].axvline(0,color='#777777',lw=.8)
    axes[0,1].grid(axis='y',alpha=.2)
    for sun,color in [(90.,'#4269a9'),(10.,'#c08b37')]:
        for p,style in [(new,'-'),(old,'--')]:
            rows = [r for r in p['angular_views'] if r['world']=='moon' and r['sun_deg']==sun and r['azimuth_from_sun_deg']==90]
            axes[1,0].plot([r['gaze_deg'] for r in rows],[r['ambient_lux']/1000 for r in rows],style,
                           color=color,label=f'Sun {sun:g}°' if style=='-' else None)
    axes[1,0].set(xlabel='Gaze elevation (degrees)',ylabel='Ambient eye illuminance (klux)',
                   title='Moon · solid: consistent sky; dashed: hybrid')
    axes[1,0].legend(frameon=False)
    comparisons = sky['additional_checks']['numerical_checks.json']['spectral_monte_carlo']
    for name,label in [('moon_1.2atm','Shielded Moon'),('earth_control','Earth control')]:
        rows = sorted([r for r in comparisons if r['world']==name and r['percent_comparison_resolved']],
                      key=lambda r:r['sun_deg'])
        axes[1,1].errorbar([r['sun_deg'] for r in rows],
                           [100*r['relative_difference'] for r in rows],
                           yerr=[200*r['relative_standard_error'] for r in rows],
                           color=colors[name],fmt='o',capsize=3,label=label)
    axes[1,1].axhline(0,color='#777777',lw=.8)
    axes[1,1].set(xlabel='Sun elevation (degrees)',ylabel='Difference from Monte Carlo (%)',
                   title='Independent Monte Carlo comparison')
    axes[1,1].legend(frameon=False)
    fig.suptitle('Solved spherical skies · titania shield and current atmospheric profiles',fontsize=15)
    fig.text(.5,.025,'Clear molecular atmosphere; ground albedo 0.1. Clouds, aerosols, refraction and polarization remain separate inputs.\n'
             'Bars: ±2 sampling standard errors; samples shown have relative standard error ≤2%. Earth below −6° needs further sampling.\n'
             'Scene geometry and surface reflectance match the earlier optical-comfort cases.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.1,1,.95))
    OUTPUT.mkdir(parents=True,exist_ok=True)
    for ext in ('png','svg'):
        fig.savefig(OUTPUT/f'spherical_sky.{ext}',dpi=180)
    manifest = dict(inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                    plotter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUTPUT/'spherical_sky.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(OUTPUT/'spherical_sky.png')


if __name__=='__main__':
    main()
