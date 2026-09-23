"""Build the self-contained A1 diagnostic lab from verified local reference data."""
from __future__ import annotations
import argparse,hashlib,json,html as html_module
from pathlib import Path
ROOT=Path(__file__).resolve().parent
# These labs display the reference computations of illumination/references.
REFERENCES=ROOT.parents[1]/'illumination'/'references'

def build(data_dir:Path,output:Path)->None:
    manifest=json.loads((data_dir/'profiles.json').read_text())
    data={'manifest':manifest,'skies':{},'profileRaw':{}}
    for r in manifest['profiles']:
        raw=(data_dir/r['file']).read_text();sha=hashlib.sha256(raw.encode()).hexdigest()
        if sha!=r['profile_sha256']:raise ValueError('Profile identity mismatch: '+r['id'])
        sky=json.loads((data_dir/(r['id']+'.sky.json')).read_text())
        if sky['profile_sha256']!=sha:raise ValueError('Sky/profile identity mismatch')
        for key,path in [('solver_sha256',REFERENCES/'tools/a1/reference_optics.py'),('spectral_sha256',REFERENCES/'data/a1/spectral-inputs.json')]:
            if sky[key]!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('Stale clear-sky reference: '+key)
        data['skies'][r['id']]=sky;data['profileRaw'][r['id']]=raw
    for key,suffix in [('cloudMC','mc'),('cloudGPU','gpu'),('cloudField','cloud')]:data[key]=json.loads((data_dir/('moon-fair.'+suffix+'.json')).read_text())
    cloud_sha=hashlib.sha256((data_dir/'moon-fair.cloud.json').read_bytes()).hexdigest()
    for key in ('cloudGPU','cloudMC'):
        if not data[key]['passed'] or data[key]['field_sha256']!=cloud_sha:raise ValueError('Cloud evidence mismatch')
    if data['cloudMC']['solver_sha256']!=hashlib.sha256((REFERENCES/'tools/a1/cloud_reference.py').read_bytes()).hexdigest():raise ValueError('Stale cloud reference')
    payload=json.dumps(data,separators=(',',':'),ensure_ascii=True).replace('</','<\\/')
    html=(ROOT/'index.accuracy-lab.html').read_text().replace('__A1_DATA__',payload).replace('__A1_SCRIPT__',(ROOT/'src/accuracy-lab.js').read_text()).replace('__A1_LICENSE__',html_module.escape((REFERENCES/'data/a1/BRUNETON_LICENSE.txt').read_text()))
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(html)
    print(json.dumps({'output':str(output),'bytes':output.stat().st_size,'sha256':hashlib.sha256(output.read_bytes()).hexdigest()}))
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data',type=Path,default=REFERENCES/'data/a1/generated');ap.add_argument('--output',type=Path,default=ROOT/'Open_Moon_A1_Accuracy_Lab.html');a=ap.parse_args();build(a.data,a.output)
