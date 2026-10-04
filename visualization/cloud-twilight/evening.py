"""Render regional maps, computed cloud views and a local interactive gallery.

MPLCONFIGDIR=/tmp/terluna-mpl python visualization/cloud-twilight/evening.py
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from scipy.ndimage import gaussian_filter

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from illumination.sky.colour_matching import M_XYZ_RGB
from shared.constants import SYNODIC_MONTH_DAYS

OUT=ROOT/'research/runs/optical_comfort/evening/figures'
RESULTS=ROOT/'research/studies/cloud_twilight/results'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def colour(xyz,scale=1500.):
    """Display-only hue-preserving compression and sRGB encoding."""
    rgb=np.maximum(xyz@M_XYZ_RGB.T,0)
    linear=rgb/(scale+rgb.max(axis=-1,keepdims=True))
    return np.where(linear<=.0031308,12.92*linear,1.055*linear**(1/2.4)-.055)


def load_scenes():
    all_scenes={};sources={}
    for mode in ('selected','history','sensitivity','deep'):
        path=RESULTS/f'evening_scenes_{mode}.json'
        product=json.loads(path.read_text());sources[str(path.relative_to(ROOT))]=digest(path)
        rows=[]
        for scene in product['scenes']:
            archive=ROOT/scene['archive']['path']
            if digest(archive)!=scene['archive']['sha256']:
                raise ValueError('Scene archive differs')
            with np.load(archive) as f:
                row=dict(scene=scene,xyz=f['xyz'],clear=f['clear_xyz'],se=f['standard_error_xyz'],
                         e=f['elevations_deg'],a=f['azimuths_deg'])
                np.testing.assert_allclose(row['xyz'],f['group_xyz'].mean(axis=2),rtol=1e-12,atol=1e-10)
            broad_colour=gaussian_filter(row['xyz'],sigma=(3.,3.,0))
            light=gaussian_filter(row['xyz'][:,:,1],sigma=1.)
            row['display_xyz']=broad_colour*light[:,:,None]/np.maximum(broad_colour[:,:,1,None],1e-30)
            np.testing.assert_allclose(row['display_xyz'][:,:,1],light,rtol=1e-12,atol=1e-10)
            rows.append(row)
        all_scenes[mode]=rows
    return all_scenes,sources


def coordinate(lat,lon):
    west=(lon+180)%360-180
    return f'{abs(lat):.1f}°{"N" if lat>=0 else "S"}, {abs(west):.1f}°{"E" if west>=0 else "W"}'


def regional_map(regional,deep=None):
    geography_path=ROOT/'research/runs/optical_comfort/evening/inputs/geography.json'
    admission=json.loads(geography_path.read_text());path=ROOT/admission['path']
    if digest(path)!=admission['sha256']:
        raise ValueError('Geography basemap differs')
    with np.load(path) as f:
        water=f['water_label']>0
    water=np.roll(water,water.shape[1]//2,axis=1)
    fig=plt.figure(figsize=(14,10.4),layout='constrained')
    fig.get_layout_engine().set(rect=(0,.055,1,.94))
    grid=fig.add_gridspec(2,2,height_ratios=[1.2,1])
    ax=fig.add_subplot(grid[0,:])
    ax.imshow(water,extent=(-180,180,-90,90),cmap=ListedColormap(['#eee9dd','#c8dae1']),interpolation='nearest')
    rows=regional['regional_bins']
    scatter=ax.scatter([(r['longitude_deg']+180)%360-180 for r in rows],
                       [r['latitude_deg'] for r in rows],c=[100*r['cloud_fraction'] for r in rows],
                       cmap='magma',vmin=0,vmax=50,s=75,edgecolor='white',linewidth=.35,zorder=3)
    for i,s in enumerate(regional['selected_scenes']+([] if deep is None else deep['selected_scenes'])):
        lon=(s['observer_longitude_deg']+180)%360-180;lat=s['observer_latitude_deg']
        ax.scatter(lon,lat,marker='*',s=170,edgecolor='#123742',facecolor='#72e2ca',zorder=4)
        ax.annotate(str(i+1),(lon,lat),xytext=(3,8),textcoords='offset points',fontsize=10,weight='bold',zorder=5,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.6))
    ax.scatter(-113.9,-44.7,marker='s',s=65,facecolor='white',edgecolor='#123742')
    ax.annotate('Highland box:\nclear evenings',(-113.9,-44.7),xytext=(-169,-32),fontsize=9,
                arrowprops=dict(arrowstyle='-',color='#455760',lw=.8))
    ax.set(xlabel='Longitude (degrees east)',ylabel='Latitude',xlim=(-180,180),ylim=(-80,80),
           title='Where evening clouds occur in the sampled tracks')
    fig.colorbar(scatter,ax=ax,fraction=.023,pad=.02,label='Cloudy columns during first 88.6 evening hours (%)')
    ax.grid(alpha=.2)
    b=fig.add_subplot(grid[1,0]);c=fig.add_subplot(grid[1,1])
    colours=['#26748a','#b15e35','#9270a7']
    labels=['Equatorial ring','Tilted ring through 45°E','Tilted ring through 135°E']
    for case,col,label in zip(regional['case_summaries'][:3],colours,labels):
        rows=case['latitude_bands']
        b.plot([0.] if case['case']=='ring_equator' else [sum(r['abs_latitude_deg'])/2 for r in rows],[100*r['cloud_fraction'] for r in rows],
               'o-',color=col,label=label)
    b.set(xlabel='Absolute latitude (degrees)',ylabel='Cloudy evening columns (%)',ylim=(0,25),xlim=(-2,75),
          title='Evening clouds concentrate near the equator')
    b.legend(fontsize=9,frameon=False);b.grid(alpha=.2)
    context=regional['darkness_context']
    hours=-np.array(context['sun_deg'])*SYNODIC_MONTH_DAYS*24/360
    mask=hours>=0
    c.semilogy(hours[mask],np.array(context['surface_lux'])[mask],color='#354a60',lw=2.5,label='Clear open surface')
    for level in (100,10,1):
        c.axhline(level,color='#8b96a1',ls=':',lw=.8)
    c.set(xlabel='Hours after equatorial geometric sunset',ylabel='Horizontal illuminance (lux)',
          title='A long, luminous twilight',xlim=(0,120),ylim=(.05,10000))
    c.grid(alpha=.2)
    fig.suptitle('Evening cloud regions on the Open Moon',fontsize=22,weight='bold')
    fig.text(.5,.003,'Second lunar cycle only • Three 2-D rings and one terrain box • Sampled columns, one cycle • Stars identify rendered views\n'
             'Basemap: hydrostatic 28% water scenario. Light curve: solved mean molecular column, design shield, ground albedo 0.1; solar light only.',
             ha='center',fontsize=9)
    fig.savefig(OUT/'regions.png',dpi=170,bbox_inches='tight');plt.close(fig)
    return {str(geography_path.relative_to(ROOT)):digest(geography_path),str(path.relative_to(ROOT)):digest(path)}


def image_panel(ax,row,scale=1500.):
    ax.imshow(colour(row['display_xyz'],scale),origin='lower',extent=(-50,50,row['e'][0],row['e'][-1]),
              interpolation='bilinear',aspect='equal')
    ax.set(xlabel='Azimuth offset from selected direction (degrees)',ylabel='Elevation (degrees)')
    s=row['scene']
    ax.text(.03,.05,f"Open ground {s['surface_lux']:.0f} lx\nBlack recess {s['recess_lux']:.2g} lx",transform=ax.transAxes,
            color='white',fontsize=9,bbox=dict(facecolor='#071c27',alpha=.7,edgecolor='none',pad=4))


def selected_views(rows):
    fig,axes=plt.subplots(3,2,figsize=(14,12.5),layout='constrained')
    fig.get_layout_engine().set(rect=(0,.07,1,.925))
    for i,(ax,row) in enumerate(zip(axes.flat,rows)):
        s=row['scene'];image_panel(ax,row)
        ax.set_title(f"{i+1}   {coordinate(s['observer_latitude_deg'],s['observer_longitude_deg'])}\n"
                     f"Sun {s['observer_sun_elevation_deg']:.1f}°; {s['initial_pick']['hours_after_sunset']:.1f} h after sunset",loc='left',fontsize=11)
    fig.suptitle('Evening clouds: broad colour and brightness',fontsize=22,weight='bold')
    fig.text(.5,.003,'Common screen exposure scale: 1,500 cd/m². Colours derive from spectral XYZ; display compression preserves hue.\n'
             'Display averages: colour 3 pixels, luminance 1 pixel; sampling noise remains. Raw values and errors are preserved. 2-D cloud sections extruded across 200 km.\n'
             'The recess admits a 20° wide opening at elevations 5–25°; its lux value describes that opening. Ground/terrain imagery is outside these views.',
             ha='center',fontsize=9)
    fig.savefig(OUT/'cloud_views.png',dpi=160,bbox_inches='tight');plt.close(fig)


def history(rows,analysis):
    fig=plt.figure(figsize=(15,11),layout='constrained')
    fig.get_layout_engine().set(rect=(0,.055,1,.94))
    grid=fig.add_gridspec(3,4,height_ratios=[1.05,1,1])
    a=fig.add_subplot(grid[0,:2]);b=fig.add_subplot(grid[0,2:])
    times=[r['hours_after_sunset'] for r in analysis]
    for key,error,label,col in [('surface_lux','surface_standard_error_lux','Open ground','#315d83'),
                                ('recess_lux','recess_standard_error_lux','Black viewing recess','#b86d38')]:
        a.errorbar(times,[r[key] for r in analysis],yerr=[2*r[error] for r in analysis],fmt='o-',color=col,label=label,capsize=3)
    a.set(yscale='log',xlabel='Hours after sunset',ylabel='Horizontal illuminance (lux)',title='Light at one fixed observer')
    a.legend(frameon=False);a.grid(alpha=.2)
    b.errorbar(times,[r['aperture']['contrast'] for r in analysis],yerr=[2*r['aperture']['contrast_standard_error'] for r in analysis],
               fmt='o-',color='#2b8273',capsize=3)
    b.axhline(0,color='#999999',lw=.8)
    b.set(xlabel='Hours after sunset',ylabel='Contrast relative to the clear sightline',title='The same viewing aperture through changing clouds')
    b.grid(alpha=.2)
    for i,row in enumerate(rows):
        ax=fig.add_subplot(grid[1+i//4,i%4]);image_panel(ax,row,300.)
        ax.set_title(f"{times[i]:.1f} h · Sun {row['scene']['observer_sun_elevation_deg']:.1f}°",fontsize=11)
        ax.set_xlabel('Azimuth offset (degrees)')
    empty=fig.add_subplot(grid[2,3]);empty.axis('off')
    empty.text(.02,.9,'Fixed site: 0°, 19.7°W\nFixed viewing direction\nEvolving saved CM1 fields\n\nCommon exposure: 300 cd/m²\nError bars: ±2 sampling SE\n\nThe interactive gallery allows\nexposure adjustment for\nthe late, faint views.\n\nFaint colour fields retain\nsubstantial sampling noise.',va='top',fontsize=11,linespacing=1.5)
    fig.suptitle('An evening cloud system through 72 hours of change',fontsize=22,weight='bold')
    fig.text(.5,.003,'Each frame uses its saved weather snapshot. The initial cloud target defines the fixed viewing direction.\n'
             'Solar light only; Earthlight is an additional source at this near-side location. Optical and across-ring geometry assumptions remain scenarios.',
             ha='center',fontsize=9)
    fig.savefig(OUT/'evening_history.png',dpi=165,bbox_inches='tight');plt.close(fig)


def deep_view(row,stats):
    fig,axes=plt.subplots(1,2,figsize=(14,4.5),layout='constrained')
    fig.get_layout_engine().set(rect=(0,.12,1,.87))
    extent=(-30,30,row['e'][0],row['e'][-1])
    for ax,xyz,title in zip(axes,[row['display_xyz'],row['clear']],['Saved cloud field + atmosphere','Corresponding clear atmosphere']):
        ax.imshow(colour(xyz,100.),origin='lower',extent=extent,aspect='equal',interpolation='bilinear')
        ax.add_patch(plt.Rectangle((-10,.5),20,4.5,fill=False,color='white',linewidth=1.3,linestyle='--'))
        ax.set(xlabel='Azimuth offset (degrees)',ylabel='Elevation (degrees)',title=title)
    fig.suptitle('A distant cloud in the deep evening',fontsize=22,weight='bold')
    fig.text(.5,.03,f"0°, 113.1°W · {stats['hours_after_sunset']:.1f} hours after sunset · Cloud top 73.8 km, about 399 km away\n"
             f"Open surface: {stats['surface_lux']:.1f} lux. Low viewing recess: {stats['recess_lux']:.3f} ± {2*stats['recess_standard_error_lux']:.3f} lux (2 SE).\n"
             'Dashed opening: 20° wide, elevations 0.5–5°. Shared exposure scale 100 cd/m²; colour averaged across 3 pixels, luminance across 1; solar source only.',
             ha='center',fontsize=10)
    fig.savefig(OUT/'deep_evening.png',dpi=170,bbox_inches='tight');plt.close(fig)


def gallery(scenes,analysis):
    data=[]
    for mode,rows in scenes.items():
        for i,row in enumerate(rows):
            s=row['scene'];stats=analysis[mode][i]
            data.append(dict(mode=mode,label=f"{mode.title()} {i+1}: {coordinate(s['observer_latitude_deg'],s['observer_longitude_deg'])}, "
                        f"Sun {s['observer_sun_elevation_deg']:.1f}°",shape=list(row['xyz'].shape[:2]),
                        xyz=row['display_xyz'].round(7).tolist(),raw_xyz=row['xyz'].round(7).tolist(),
                        clear=row['clear'].round(7).tolist(),se=row['se'].round(7).tolist(),e=row['e'].tolist(),
                        a=(row['a']-s['camera_azimuth_in_ring_basis_deg']).tolist(),stats=stats))
    template=Path(__file__).with_name('evening.html').read_text()
    text=template.replace('/*SCENES*/[]',json.dumps(data,separators=(',',':'),allow_nan=False))
    (OUT/'index.html').write_text(text)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    regional_path=RESULTS/'regional_evenings.json';appearance_path=RESULTS/'evening_appearance.json'
    regional=json.loads(regional_path.read_text());appearance=json.loads(appearance_path.read_text())
    deep_path=RESULTS/'deep_evenings.json';deep=json.loads(deep_path.read_text())
    scenes,sources=load_scenes()
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in (regional_path,appearance_path,deep_path,Path(__file__),Path(__file__).with_name('evening.html'))})
    sources.update(regional_map(regional,deep));selected_views(scenes['selected']);history(scenes['history'],appearance['history'])
    deep_view(scenes['deep'][0],appearance['deep'][0])
    gallery(scenes,appearance)
    (OUT/'manifest.json').write_text(json.dumps(dict(schema='terluna.visualization.evening-clouds/1',inputs=sources,
        outputs={p.name:digest(p) for p in OUT.iterdir() if p.suffix in ('.png','.html')},
        evidence='Computed cloud radiance and regional sampled occurrence. Display smoothing and exposure act on saved XYZ; '
                 'the physical values retain their sampling errors. Rendering checks recover per-pixel means from photon blocks.'),indent=2)+'\n')
    print(OUT)


if __name__=='__main__':
    main()
