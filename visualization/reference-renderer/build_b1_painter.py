"""Build the offline B1 viewer from actually computed frames; no fabricated fill."""
from pathlib import Path
import base64,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'tools/b1'))
from render_scene import SCENARIOS,raw_identity
import scene_transport

def enc(x):return base64.b64encode(np.asarray(x,dtype='<f4').tobytes()).decode()

def main():
    gen=ROOT/'data/b1/generated';frames=[];domes=[];maps=[];controls=[]
    for p in sorted((gen/'final').glob('*.json')):
        m=json.loads(p.read_text());d=np.load(p.with_suffix('.npz'))
        if m['sources']!=raw_identity():raise RuntimeError('Stale frame source identity: '+p.name)
        item=dict(meta=m,raw=enc(d['linear_srgb']),filtered=enc(d['filtered_linear_srgb']),summary=dict(median_Y=float(np.median(d['XYZ'][:,:,1])),surface_fraction=float(d['surface_path_XYZ'][:,:,1].sum()/max(1e-30,d['XYZ'][:,:,1].sum()))))
        if m['camera']['projection']=='perspective':frames.append(item)
        else:domes.append(dict(scenario=m['scenario'],sun=m['sun_elevation_deg'],width=m['width'],height=m['height'],raw=item['raw'],filtered=item['filtered']))
    for p in sorted((gen/'controls').glob('*.json')):
        m=json.loads(p.read_text());d=np.load(p.with_suffix('.npz'))
        if m['sources']!=raw_identity():raise RuntimeError('Stale control identity: '+p.name)
        controls.append(dict(meta=m,raw=enc(d['linear_srgb']),filtered=enc(d['filtered_linear_srgb']),summary=dict(median_Y=float(np.median(d['XYZ'][:,:,1])),surface_fraction=float(d['surface_path_XYZ'][:,:,1].sum()/max(1e-30,d['XYZ'][:,:,1].sum())))))
    spec=json.loads(scene_transport.locate('data/a1/spectral-inputs.json').read_text());M=np.array(spec['xyz_to_linear_srgb'])
    for p in sorted((gen/'maps').glob('*.json')):
        m=json.loads(p.read_text());d=np.load(p.with_suffix('.npz'));maps.append(dict(scenario=m['scenario'],sun=m['sun_elevation_deg'],size=m['size'],direct=enc(d['direct_XYZ']@M.T),diffuse=enc(d['diffuse_XYZ']@M.T),total=enc(d['total_XYZ']@M.T),shadow=enc(d['direct_beam_visibility'])))
    tests=json.loads((ROOT/'validation/b1/boundary-tests.json').read_text())
    payload=dict(schema='open-moon-b1-painter/1',scenarios=SCENARIOS,frames=frames,domes=domes,maps=maps,controls=controls,validation=dict(tests_passed=tests['passed'],tests_failed=tests['failed']))
    if not frames or any(f['meta']['unresolved_paths'] for f in frames):raise RuntimeError('Missing frames or unfinished path histories')
    raw=json.dumps(payload,separators=(',',':'),allow_nan=False).replace('</','<\\/')
    html=(ROOT/'index.b1-scene-painter.html').read_text().replace('__B1_DATA__',raw).replace('__B1_SCRIPT__',(ROOT/'src/b1-scene-painter.js').read_text())
    out=ROOT/'Open_Moon_B1_Scene_Painter.html';out.write_text(html);print(out,len(html),hashlib.sha256(out.read_bytes()).hexdigest())
if __name__=='__main__':main()
