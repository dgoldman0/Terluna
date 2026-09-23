"""Generate coordinated albedo/AO and normal/roughness/height surface atlases.

Run `python bake/surfaces.py` (or with --if-stale to skip matching textures).

All inputs are deterministic authored morphology. Heights and tile extents are
in metres. Albedo is encoded sRGB; other channels are linear. The source and
manifest make the maps editable and reproducible without external assets.
"""
from pathlib import Path
import hashlib, json
import numpy as np
from scipy.ndimage import gaussian_filter
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'assets/surfaces'; OUT.mkdir(parents=True,exist_ok=True)
if '--if-stale' in __import__('sys').argv:
    # Skip regeneration when every texture already matches the committed manifest.
    files=json.loads((OUT/'manifest.json').read_text())['files']
    if all((OUT/n).is_file() and hashlib.sha256((OUT/n).read_bytes()).hexdigest()==v['sha256'] for n,v in files.items()):
        print('Surface textures match their manifest.');raise SystemExit(0)
N=512; PAD=16; TILE=N+2*PAD
Y,X=np.mgrid[:N,:N].astype(float)
rng=np.random.default_rng(21891)
def smooth_noise(sigma):
    a=gaussian_filter(rng.standard_normal((N,N)),sigma,mode='wrap')
    return (a-a.mean())/(a.std()+1e-12)
def srgb(a):
    a=np.clip(a,0,1)
    return np.where(a<=.0031308,a*12.92,1.055*a**(1/2.4)-.055)
def ellipsoid(h,colour,cx,cy,rx,ry,angle,height,tint):
    margin=int(max(rx,ry)+2)
    ix=np.arange(int(cx)-margin,int(cx)+margin+1);iy=np.arange(int(cy)-margin,int(cy)+margin+1)
    xx,yy=np.meshgrid(ix-cx,iy-cy);u=(xx*np.cos(angle)+yy*np.sin(angle))/rx;v=(-xx*np.sin(angle)+yy*np.cos(angle))/ry
    r=u*u+v*v; mask=r<1; dome=np.sqrt(np.clip(1-r,0,1))*height
    sy=iy%N;sx=ix%N;old=h[np.ix_(sy,sx)];local=colour[np.ix_(sy,sx)]
    edge=np.clip((1-r)*9,0,1);old=np.maximum(old,dome)
    local=local*(1-edge[...,None])+(np.array(tint)[None,None,:]*(.83+.17*dome[...,None]/max(height,1e-9)))*edge[...,None]
    h[np.ix_(sy,sx)]=old;colour[np.ix_(sy,sx)]=local
specs=[
 ('sand',.75,(.36,.31,.23),.88,.006),
 ('gravel',1.5,(.25,.22,.17),.86,.060),
 ('basalt',2.8,(.105,.105,.093),.84,.090),
 ('soil',1.1,(.135,.092,.049),.94,.024),
 ('litter',1.1,(.090,.061,.029),.96,.045),
 ('bark',1.1,(.115,.076,.040),.91,.024)
]
albedo=np.zeros((TILE*2,TILE*3,4),dtype=np.uint8); packed=np.zeros_like(albedo);records=[]
for layer,(name,metres,base,rough,height_range) in enumerate(specs):
    broad=smooth_noise(22);fine=smooth_noise(1);meso=smooth_noise(5)
    colour=np.array(base)[None,None,:]*(1+.075*broad[...,None]+.045*fine[...,None])
    h=(.00012*fine+.00020*meso)
    if name=='sand':
        # Sorting laminations, tiny mineral grains and irregular ripple crests.
        warp=.55*smooth_noise(35);phase=Y/29+X/112+warp
        h+=.0015*(.5+.5*np.cos(phase*2*np.pi))**5
        grains=rng.uniform(size=(N,N));colour*=1+((grains-.5)*.23)[...,None]
        dark=grains>.995;colour[dark]*=.38
        for _ in range(210):
            ellipsoid(h,colour,*rng.uniform(0,N,2),rng.uniform(1,3.6),rng.uniform(1,2.7),rng.uniform(0,6.28),rng.uniform(.0008,.0022),rng.uniform(.14,.28,3))
    elif name=='gravel':
        colour*=.76;h+=.0006*meso
        for _ in range(1450):
            r=rng.uniform(2,11);v=rng.uniform(.07,.32);tint=(v,v*.94,v*.80)
            ellipsoid(h,colour,*rng.uniform(0,N,2),r,r*rng.uniform(.5,.9),rng.uniform(0,6.28),r/N*metres*rng.uniform(.3,.65),tint)
    elif name=='basalt':
        h+=.004*smooth_noise(4)+.008*smooth_noise(29)
        fissure=np.ones((N,N))
        for _ in range(14):
            a=rng.uniform(0,np.pi);off=rng.uniform(-N,N);dist=np.abs((X-N/2)*np.cos(a)+(Y-N/2)*np.sin(a)-off+2.0*meso)
            fissure=np.minimum(fissure,np.clip(dist/2,0,1))
        h-=.020*(1-fissure);colour*=.48+.52*fissure[...,None]
        mineral=np.clip(smooth_noise(.4)-1.4,0,1)*.10;colour+=mineral[...,None]*np.array([.65,.63,.51])
        colour*=1+.09*meso[...,None]
    elif name in ('soil','litter'):
        h+=.0015*meso+.0012*fine
        for _ in range(750):
            r=rng.uniform(1.0,8);v=rng.uniform(.055,.19)
            ellipsoid(h,colour,*rng.uniform(0,N,2),r,r*.72,rng.uniform(0,6.28),rng.uniform(.0005,.008),(v,v*.67,v*.35))
        if name=='litter':
            for _ in range(155):
                cx,cy=rng.uniform(0,N,2);ang=rng.uniform(0,6.28);rx=rng.uniform(9,27);ry=rx*rng.uniform(.25,.50);margin=int(rx+3)
                ix=np.arange(int(cx)-margin,int(cx)+margin+1);iy=np.arange(int(cy)-margin,int(cy)+margin+1);xx,yy=np.meshgrid(ix-cx,iy-cy)
                u=(xx*np.cos(ang)+yy*np.sin(ang))/rx;v=(-xx*np.sin(ang)+yy*np.cos(ang))/ry
                outline=(u*u+v*v<1)&(np.abs(v)<np.maximum(0,1-np.abs(u))**.60)
                veins=np.exp(-np.abs(v)*24)+.34*np.exp(-np.abs(np.sin((u*5+np.abs(v)*1.7)*np.pi))*24)
                sy=iy%N;sx=ix%N;local=colour[np.ix_(sy,sx)];value=rng.uniform(.075,.23);tint=np.array([value,value*rng.uniform(.53,.72),value*.25])
                leaf=tint[None,None,:]*(.80+.20*np.maximum(0,1-u*u-v*v)[...,None])*(1-.21*veins[...,None])
                local[outline]=leaf[outline];colour[np.ix_(sy,sx)]=local
                old=h[np.ix_(sy,sx)];old[outline]=np.maximum(old[outline],.010+.006*np.abs(v[outline]));h[np.ix_(sy,sx)]=old
    else:
        warp=2.2*smooth_noise(22);ridge=(.5+.5*np.cos((X/17+warp*.07)*2*np.pi))**.8
        pores=smooth_noise(2);h+=.013*ridge+.0009*pores
        cracks=np.clip(smooth_noise(6)-1.2,0,1);h-=.008*cracks;colour*=.46+.75*ridge[...,None];colour*=1-.30*cracks[...,None]
    # Positive and negative micro relief share one height source across channels.
    dx=(np.roll(h,-1,1)-np.roll(h,1,1))/(2*metres/N)
    dy=(np.roll(h,-1,0)-np.roll(h,1,0))/(2*metres/N)
    norm=np.sqrt(dx*dx+dy*dy+1);nx=-dx/norm;ny=-dy/norm
    cavity=np.clip((gaussian_filter(h,2,mode='wrap')-h)/max(.001,height_range*.10),0,1)
    ao=1-.44*cavity; roughness=np.clip(rough+.035*meso-.045*np.clip(h/height_range,0,1),.38,1)
    ar=np.empty((N,N,4),np.uint8);ar[...,:3]=(srgb(colour)*255+.5).astype(np.uint8);ar[...,3]=(ao*255+.5).astype(np.uint8)
    nr=np.stack([nx*.5+.5,ny*.5+.5,roughness,np.clip(.5+h/height_range,0,1)],axis=-1)
    nr=(nr*255+.5).astype(np.uint8)
    ar=np.pad(ar,((PAD,PAD),(PAD,PAD),(0,0)),mode='wrap');nr=np.pad(nr,((PAD,PAD),(PAD,PAD),(0,0)),mode='wrap')
    row,col=divmod(layer,3);albedo[row*TILE:(row+1)*TILE,col*TILE:(col+1)*TILE]=ar;packed[row*TILE:(row+1)*TILE,col*TILE:(col+1)*TILE]=nr
    records.append({'name':name,'layer':layer,'tile_metres':metres,'height_range_metres':height_range,'mean_albedo_linear':np.clip(colour,0,1).mean(axis=(0,1)).tolist(),'mean_roughness':float(roughness.mean()),'height_min_max_metres':[float(h.min()),float(h.max())]})
Image.fromarray(albedo).save(OUT/'albedo-ao.png',optimize=True)
Image.fromarray(packed).save(OUT/'normal-roughness-height.png',optimize=True)
manifest={'schema':'open-moon-surfaces/1','generator':'bake/surfaces.py','seed':21891,'provenance':'Original deterministic morphology; no external image assets','license':"Project author's licensing choice",'tile_pixels':N,'padding_pixels':PAD,'layout':[3,2],'size_pixels':[TILE*3,TILE*2],'albedo_encoding':'sRGB RGB, linear AO alpha','normal_encoding':'tangent XY in RG, linear roughness B, height/range + 0.5 in A','layers':records,'files':{}}
for name in ['albedo-ao.png','normal-roughness-height.png']:
    b=(OUT/name).read_bytes();manifest['files'][name]={'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest['files'],indent=2))
