"""Capture A2 cloud grids from the exact revision-04 density field.

The suite separates voxel-resolution convergence from finite-crop convergence.
No production renderer source is modified.
"""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, time
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[2]
# The experience whose real-time cloud rendering these references measure.
VIEWER=ROOT.parents[1]/'immersion'/'dist'
ENTRY='experiences/shoreline/index.html'

def serve_viewer():
    """Serve the built experience (run `npm run build` in immersion/) on localhost."""
    import functools, threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
    class Quiet(SimpleHTTPRequestHandler):
        def log_message(self, *args): pass
    server=ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Quiet,directory=str(VIEWER)))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    return server,f'http://127.0.0.1:{server.server_port}/'
CONFIGS=[
    ('grid-coarse',[25,17,25],1.0),
    ('grid-a1',[49,33,49],1.0),
    ('crop-075',[37,33,37],0.75),
    ('crop-150',[73,33,73],1.5),
    ('crop-200',[97,33,97],2.0),
    ('grid-fine',[97,65,97],1.0),
]
def canonical_field(field): return (json.dumps(field,sort_keys=True,separators=(',',':'))+'\n').encode()
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'data/a2/generated');ap.add_argument('--probe-count',type=int,default=96);ap.add_argument('--only',nargs='*');a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    html=(VIEWER/ENTRY).read_bytes();server,base=serve_viewer();html_sha=hashlib.sha256(html).hexdigest();profiles=json.loads((ROOT/'data/a1/generated/profiles.json').read_text());profile=next(x for x in profiles['profiles'] if x['id']=='moon-fair')
    report={'schema':'open-moon-a2-cloud-capture/1','html_sha256':html_sha,'profile_sha256':profile['profile_sha256'],'configs':[],'errors':[],'warnings':[],'external_requests':[],'passed':False}
    xvfb=None;start=time.monotonic()
    try:
      if not Path('/tmp/.X11-unix/X98').exists():xvfb=subprocess.Popen(['Xvfb',':98','-screen','0','1100x850x24','-nolisten','tcp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);time.sleep(.5)
      with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=False,env={**os.environ,'DISPLAY':':98'},args=['--no-sandbox','--disable-dev-shm-usage','--disable-gpu-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':800,'height':600});page.set_default_timeout(180000);page.on('pageerror',lambda e:report['errors'].append(str(e)))
        def console(m):
          if m.type=='error' or any(x in m.text.lower() for x in ['context lost','context_lost','shader error']):report['errors'].append(m.text)
          elif m.type=='warning':report['warnings'].append(m.text)
        page.on('console',console);page.on('request',lambda r:report['external_requests'].append(r.url) if not r.url.startswith(('blob:','data:',base)) else None)
        page.expose_function('a2Digest',lambda name,values:list(hashlib.new({'SHA-256':'sha256','SHA-1':'sha1'}[name],bytes(values)).digest()))
        page.evaluate("""()=>{if(!crypto.subtle)Object.defineProperty(window,'crypto',{value:{subtle:{digest:async(name,b)=>new Uint8Array(await a2Digest(name,Array.from(new Uint8Array(b.buffer||b,b.byteOffset||0,b.byteLength)))).buffer}}});window.requestAnimationFrame=cb=>{window.__a2Frame=cb;return 1;};}""")
        page.goto(base+ENTRY,wait_until='domcontentloaded');page.wait_for_function('window.openMoonShoreline?.ready || window.openMoonShorelineError',polling=100);err=page.evaluate('window.openMoonShorelineError||null');
        if err:raise RuntimeError(err)
        page.add_script_tag(content=(ROOT/'src/frozen-cloud-field.js').read_text());probe=(ROOT/'tests/a2_frozen_probe.js').read_text()
        for name,dims,mult in CONFIGS:
          if a.only and name not in a.only: continue
          before=len(report['errors']);result=page.evaluate(probe,{'regime':'fair','world':'moon','dimensions':dims,'extentMultiplier':mult,'probeExtentMultiplier':min(1.0,mult),'probeCount':a.probe_count});page.wait_for_timeout(150)
          if len(report['errors'])!=before:raise RuntimeError(f'Browser errors during {name}')
          field=result.pop('field');field['provenance']['html_sha256']=html_sha;field['provenance']['shared_atmospheric_profile_sha256']=profile['profile_sha256'];raw=canonical_field(field);fn=f'moon-fair-{name}.cloud.json';(a.output/fn).write_bytes(raw)
          entry={'id':name,'dimensions':dims,'extent_multiplier':mult,'field_file':fn,'field_sha256':hashlib.sha256(raw).hexdigest(),**result};report['configs'].append(entry);(a.output/f'capture-{name}.json').write_text(json.dumps({'schema':'open-moon-a2-cloud-capture-item/1','html_sha256':html_sha,'profile_sha256':profile['profile_sha256'],'entry':entry},indent=2)+'\n');print(json.dumps({'id':name,'voxel':entry['voxelization']['sum_absolute_difference_over_sum_continuous'],'elapsed_s':time.monotonic()-start}),flush=True)
        if report['external_requests']:raise RuntimeError('External requests recorded')
        if any(x['field_parity']['maximum_absolute_extinction_error_m1']>1e-7 for x in report['configs']):raise RuntimeError('Frozen CPU/GPU field sampling mismatch')
        report['passed']=True;browser.close()
    except Exception as e:report['errors'].append(type(e).__name__+': '+str(e))
    finally:
      if xvfb:xvfb.terminate()
      report['elapsed_s']=time.monotonic()-start;(a.output/'cloud-capture.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':report['passed'],'errors':report['errors'],'elapsed_s':report['elapsed_s']},indent=2))
    return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
