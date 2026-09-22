"""Assemble the self-contained landscape viewer from editable source and assets."""
from pathlib import Path
import base64, hashlib, json
ROOT=Path(__file__).resolve().parent
html=(ROOT/'index.landscape.html').read_text()
html=html.replace('__SKY_DATA__',(ROOT/'data/atmosphere.json').read_text())
vendor=ROOT/'vendor/three.cjs'
if vendor.exists():
 b=vendor.read_bytes();h=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
 if h!='ca4833532c363b72477b2e8a6f47cc0e2fc7b09a':raise ValueError('Three.js library hash mismatch')
 encoded=base64.b64encode(b).decode()
else:encoded=''
html=html.replace('__THREE_VENDOR__',encoded)
assets=ROOT/'assets/surfaces';manifest=json.loads((assets/'manifest.json').read_text())
images={}
for key,name in [('albedo','albedo-ao.png'),('packed','normal-roughness-height.png')]:
 b=(assets/name).read_bytes()
 if hashlib.sha256(b).hexdigest()!=manifest['files'][name]['sha256']:raise ValueError('Surface asset hash mismatch: '+name)
 images[key]='data:image/png;base64,'+base64.b64encode(b).decode()
images['manifest']=manifest
html=html.replace('__SURFACE_ASSETS__',json.dumps(images,separators=(',',':')))
for name in ['landscape','core','atmosphere','materials','terrain','surface-water','ecology','scene','water','audio','app','loader']:
 html=html.replace('__'+name.upper()+'__',(ROOT/'src'/f'{name}.js').read_text().replace('</script','<\\/script'))
if any(t in html for t in ['__SURFACE_ASSETS__','__LANDSCAPE__','__TERRAIN__']):raise ValueError('Unexpanded build token')
out=ROOT/'Open_Moon_Shoreline.html';out.write_text(html)
print(f'{out.name}: {out.stat().st_size:,} bytes; dependency embedded: {vendor.exists()}')
