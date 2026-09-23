"""Build the offline r186 viewer. Verify original dependency/assets before linking.

No dependency download is performed by the builder. tools/import_three.py imports
just the required release files from the supplied ZIP, retaining their hashes.
"""
from pathlib import Path
import base64, hashlib, json, re
ROOT=Path(__file__).resolve().parent

def dependency(include_webgpu=False):
    folder=ROOT/'vendor/three-r186'
    manifest=json.loads((folder/'manifest.json').read_text())
    if manifest['version']!='0.186.0':
        raise ValueError('Expected Three.js 0.186.0')
    names=['three.core.js','three.module.js']
    if include_webgpu:names+=['three.webgpu.js','three.tsl.js']
    modules={}
    for name in names:
        data=(folder/name).read_bytes();spec=manifest['files'][name]
        if len(data)!=spec['bytes'] or hashlib.sha256(data).hexdigest()!=spec['sha256']:
            raise ValueError('Dependency integrity mismatch: '+name)
        modules[name]=base64.b64encode(data).decode()
    return json.dumps({'manifest':manifest,'modules':modules},separators=(',',':'))

def build(template='index.landscape.html',output='Open_Moon_Shoreline.html',webgpu=False):
    html=(ROOT/template).read_text()
    html=html.replace('__SKY_DATA__',(ROOT/'data/atmosphere.json').read_text())
    html=html.replace('__THREE_VENDOR__',dependency(webgpu))
    assets=ROOT/'assets/surfaces';manifest=json.loads((assets/'manifest.json').read_text())
    images={'manifest':manifest}
    for key,name in [('albedo','albedo-ao.png'),('packed','normal-roughness-height.png')]:
        data=(assets/name).read_bytes()
        if hashlib.sha256(data).hexdigest()!=manifest['files'][name]['sha256']:
            raise ValueError('Surface asset integrity mismatch: '+name)
        images[key]='data:image/png;base64,'+base64.b64encode(data).decode()
    html=html.replace('__SURFACE_ASSETS__',json.dumps(images,separators=(',',':')))
    for name in ['landscape','core','weather-column','cloud-optics','cloud-renderer','atmosphere','materials','terrain','surface-water','ponds','ecology','scene','water','audio','app','renderer-lab','node-materials','loader']:
        token='__'+name.upper()+'__'
        # The column thermodynamics are the atmosphere domain's; everything else is the viewer's own.
        source=ROOT.parents[2]/'atmosphere'/'column'/'weather-column.js' if name=='weather-column' else ROOT/'src'/f'{name}.js'
        if token in html:html=html.replace(token,source.read_text().replace('</script','<\\/script'))
    unexpanded=re.findall(r'__[A-Z][A-Z_-]+__',html)
    if unexpanded:raise ValueError('Unexpanded tokens: '+str(unexpanded))
    path=ROOT/output;path.write_text(html)
    print(f'{output}: {path.stat().st_size:,} bytes; Three.js r186 ESM embedded')
    return path
if __name__=='__main__':
    build()
    if (ROOT/'index.renderer-lab.html').exists():build('index.renderer-lab.html','Open_Moon_Renderer_Lab.html',True)
