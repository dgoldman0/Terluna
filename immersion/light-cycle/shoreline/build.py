"""Build the complete viewer from checked-in scene code and packed sky data.
A downloaded dependency is optional: without it the HTML fetches the pinned
Three.js library at first opening and offers a self-contained offline export.
"""
from pathlib import Path
import base64,hashlib
ROOT=Path(__file__).resolve().parent
html=(ROOT/'index.template.html').read_text()
html=html.replace('__SKY_DATA__',(ROOT/'data/atmosphere.json').read_text())
vendor=ROOT/'vendor/three.cjs'
if vendor.exists():
 b=vendor.read_bytes(); h=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
 if h!='ca4833532c363b72477b2e8a6f47cc0e2fc7b09a':raise ValueError('Three.js library hash mismatch')
 encoded=base64.b64encode(b).decode()
else: encoded=''
html=html.replace('__THREE_VENDOR__',encoded)
for name in ['core','atmosphere','scene','water','audio','app','loader']:
 html=html.replace('__'+name.upper()+'__',(ROOT/'src'/f'{name}.js').read_text().replace('</script','<\\/script'))
out=ROOT/'Open_Moon_Shoreline.html';out.write_text(html)
print(f'{out.name}: {out.stat().st_size:,} bytes; dependency embedded: {vendor.exists()}')
