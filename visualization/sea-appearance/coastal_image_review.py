"""Image-space checks for the selected cloud illustrations.

These compare composition and display values; they do not recover radiance,
validate human colour perception, or measure wave height from a photograph.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
from scipy.stats import spearmanr


def luminance(path):
    rgb=np.asarray(Image.open(path).convert('RGB'),float)/255
    linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
    return linear@np.array([.2126729,.7151522,.072175])


def measure(path,display):
    y=luminance(path);height,width=y.shape;camera=display['camera']
    rgb=np.asarray(Image.open(path).convert('RGB'),float)/255
    linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
    focal=width/(2*np.tan(np.radians(camera['horizontal_fov_deg']/2)))
    horizon=height/2+focal*np.tan(np.radians(camera['pitch_deg']))
    lo=max(0,int(horizon-.035*height));hi=min(height,int(horizon+.035*height))
    # A coherent colour edge avoids mistaking individual wave crests for the horizon.
    mean_rgb=gaussian_filter(linear.mean(axis=1),(.7,0))
    edge=np.linalg.norm(np.diff(mean_rgb,axis=0),axis=1)
    horizon_row=lo+int(np.argmax(edge[lo:hi]))
    skyline=np.full(width,horizon_row/height)
    regions={'upper_sky':y[:max(1,int(.12*height))],
             'middle_cloud_field':y[int(.27*height):int(.55*height),int(.3*width):int(.7*width)],
             'water':y[int(.91*height):]}
    rec=dict(frame=[width,height],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        expected_geometric_horizon_fraction=float(horizon/height),
        fitted_horizon_coherent_fraction=float(np.median(skyline)),
        median_linear_display_Y={k:float(np.median(v)) for k,v in regions.items()},
        median_linear_display_RGB=dict(upper_sky=np.median(linear[:max(1,int(.12*height))].reshape(-1,3),axis=0).tolist(),
            middle_cloud_field=np.median(linear[int(.27*height):int(.55*height),int(.3*width):int(.7*width)].reshape(-1,3),axis=0).tolist(),
            water=np.median(linear[int(.91*height):].reshape(-1,3),axis=0).tolist()))
    texture=None
    if not display.get('earth_below_horizon',False):
        proj=display['earth_projection'];rw,rh=camera['frame']
        cx,cy=np.array(proj['centre'])*[width/rw,height/rh]
        rx,ry=np.array(proj['disk_dimensions_px'])/2*[width/rw,height/rh]
        phi=np.linspace(0,2*np.pi,160,endpoint=False);rad=np.linspace(.65,1.45,241)
        xx=cx+rx*np.cos(phi)[:,None]*rad;yy=cy+ry*np.sin(phi)[:,None]*rad
        sample=map_coordinates(gaussian_filter(y,.5),[yy,xx],order=1)
        edge=rad[np.argmax(-np.gradient(sample,rad,axis=1),axis=1)]
        ex=cx+rx*np.cos(phi)*edge;ey=cy+ry*np.sin(phi)*edge
        fit=least_squares(lambda q:np.hypot((ex-cx-q[0])/rx,(ey-cy-q[1])/ry)-q[2],
            [0,0,1],loss='soft_l1',f_scale=.015,bounds=([-rx*.35,-ry*.35,.6],[rx*.35,ry*.35,1.4]))
        u,v=np.meshgrid(np.linspace(-.8,.8,81),np.linspace(-.8,.8,81));inside=u*u+v*v<.64
        texture=map_coordinates(y,[cy+fit.x[1]+ry*fit.x[2]*v[inside],cx+fit.x[0]+rx*fit.x[2]*u[inside]],order=1)
        rec['earth_limb_fit']=dict(centre_shift_px=fit.x[:2].tolist(),radius_scale=float(fit.x[2]),
            radial_rms_fraction=float(np.sqrt(np.mean(fit.fun**2))))
    return rec,skyline,texture


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--display',type=Path,required=True)
    p.add_argument('--image',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();display=json.loads(a.display.read_text());ref,line,texture=measure(a.reference,display);records={}
    for path in a.image:
        rec,row,tex=measure(path,display)
        target=np.interp(np.linspace(0,1,len(row)),np.linspace(0,1,len(line)),line)
        rec['horizon_rms_fraction_of_height']=float(np.sqrt(np.mean((row-target)**2)))
        rec['region_display_luminance_ratios']={k:rec['median_linear_display_Y'][k]/max(v,1e-12) for k,v in ref['median_linear_display_Y'].items()}
        if tex is not None:rec['registered_earth_interior_rank_correlation']=float(spearmanr(texture,tex).statistic)
        records[path.name]=rec
    result=dict(schema='terluna.visualization.coastal-image-checks/1',reference=ref,images=records,
        method='Coherent row-mean colour edge in the known horizon band and robust joint Earth-limb centre/radius fit; registered interior ranks; fixed-region linear display luminance.',
        limits=['Image-space screens only; no absolute photometry, human-vision validation or recovered wave height.',
                'Cloud microtexture and cross-ring morphology are illustrative; the cloud experiment does not resolve them.',
                'Limb fits can respond to glare and low pixel count; interpret together with visual inspection.'],
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
