"""Capture a frozen production cloud and exercise its solver-neutral GPU sampler.

Every console context-loss message is a failure. Browser policy and production
source remain unchanged; the harness provides SHA-256 for its opaque origin.
"""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, time
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[2]
# The experience whose real-time cloud rendering these references measure.
VIEWER=ROOT.parents[1]/'immersion'/'light-cycle'/'shoreline'

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--regime',choices=['fair','convection','fog'],default='fair')
    ap.add_argument('--world',choices=['moon','earth','moon_no_ozone'],default='moon')
    ap.add_argument('--output',type=Path,default=ROOT/'data/a1/generated')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    html=(VIEWER/'Open_Moon_Shoreline.html').read_bytes()
    report={'schema':'open-moon-a1-gpu-field-capture/1','html_sha256':hashlib.sha256(html).hexdigest(),'errors':[],'warnings':[],'external_requests':[],'passed':False}
    xvfb=None;start=time.monotonic()
    try:
        if not Path('/tmp/.X11-unix/X97').exists():
            xvfb=subprocess.Popen(['Xvfb',':97','-screen','0','1000x800x24','-nolisten','tcp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);time.sleep(.5)
        with sync_playwright() as p:
            browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=False,env={**os.environ,'DISPLAY':':97'},args=['--no-sandbox','--disable-dev-shm-usage','--disable-gpu-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
            page=browser.new_page(viewport={'width':800,'height':600});page.set_default_timeout(120000)
            page.on('pageerror',lambda e:report['errors'].append(str(e)))
            def console(m):
                if m.type=='error' or any(x in m.text.lower() for x in ['context lost','context_lost','shader error']):report['errors'].append(m.text)
                elif m.type=='warning':report['warnings'].append(m.text)
            page.on('console',console)
            page.on('request',lambda r:report['external_requests'].append(r.url) if not r.url.startswith(('blob:','data:')) else None)
            page.expose_function('a1Digest',lambda name,values:list(hashlib.new({'SHA-256':'sha256','SHA-1':'sha1'}[name],bytes(values)).digest()))
            page.evaluate("""()=>{if(!crypto.subtle)Object.defineProperty(window,'crypto',{value:{subtle:{digest:async(name,b)=>new Uint8Array(await a1Digest(name,Array.from(new Uint8Array(b.buffer||b,b.byteOffset||0,b.byteLength)))).buffer}}});window.requestAnimationFrame=cb=>{window.__a1Frame=cb;return 1;};}""")
            page.set_content(html.decode(),wait_until='domcontentloaded')
            page.wait_for_function('window.openMoonShoreline?.ready || window.openMoonShorelineError',polling=100)
            error=page.evaluate('window.openMoonShorelineError||null')
            if error:raise RuntimeError(error)
            page.add_script_tag(content=(ROOT/'src/frozen-cloud-field.js').read_text())
            result=page.evaluate((ROOT/'tests/a1_frozen_probe.js').read_text(),{'regime':a.regime,'world':a.world})
            page.wait_for_timeout(300)
            if report['errors'] or report['external_requests']:raise RuntimeError('Browser errors or external dependencies recorded')
            field=result.pop('field');field['provenance']['html_sha256']=report['html_sha256']
            profiles=json.loads((a.output/'profiles.json').read_text());pr=next(x for x in profiles['profiles'] if x['id']==a.world+'-'+a.regime)
            field['provenance']['shared_atmospheric_profile_sha256']=pr['profile_sha256']
            raw=(json.dumps(field,sort_keys=True,separators=(',',':'))+'\n').encode();field_name=a.world+'-'+a.regime+'.cloud.json'
            (a.output/field_name).write_bytes(raw)
            report.update(result);report['field_file']=field_name;report['field_sha256']=hashlib.sha256(raw).hexdigest()
            if result['field_parity']['maximum_absolute_extinction_error_m1']>1e-7:raise RuntimeError('GPU/CPU frozen sampling mismatch')
            report['passed']=True;browser.close()
    except Exception as e:
        report['errors'].append(type(e).__name__+': '+str(e))
    finally:
        if xvfb:xvfb.terminate()
        report['elapsed_s']=time.monotonic()-start
        (a.output/(a.world+'-'+a.regime+'.gpu.json')).write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:report[k] for k in ['passed','errors','elapsed_s']},indent=2))
    return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
