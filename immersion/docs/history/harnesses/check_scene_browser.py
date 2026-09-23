"""Pending end-to-end check, for an environment able to load genuine Three.js.

Serve the complete HTML over localhost. Use actual Chromium WebGL; record failure
rather than substituting a renderer. Install Playwright/Chromium separately when
running this optional check. This script has not been run successfully here.
"""
import functools,json,threading
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
server=ThreadingHTTPServer(('127.0.0.1',0),functools.partial(SimpleHTTPRequestHandler,directory=str(ROOT)))
threading.Thread(target=server.serve_forever,daemon=True).start()
record={'scope':__doc__,'success':False,'scenes':[],'errors':[]}
try:
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True)
        page=browser.new_page(viewport={'width':1280,'height':800})
        page.on('pageerror',lambda e:record['errors'].append(str(e)))
        page.goto(f'http://127.0.0.1:{server.server_port}/Open_Moon_Shoreline.html',wait_until='domcontentloaded')
        page.wait_for_function('window.openMoonShoreline?.ready || window.openMoonShorelineError',timeout=90000)
        error=page.evaluate('window.openMoonShorelineError || null')
        if error:raise RuntimeError(error)
        page.locator('#help-close').click()
        for name,phase,t,location in [('noon',0,0,'shore'),('rain',0,7800,'shelter'),('after-rain',0,14400,'path'),('sunset',.25,0,'overlook')]:
            page.evaluate('([p,t,l])=>{openMoonShoreline.setPhase(p);openMoonShoreline.setWeather(t);openMoonShoreline.setLocation(l);}',[phase,t,location])
            page.wait_for_timeout(1200)
            assert not page.evaluate('window.openMoonShorelineError || null')
            result=page.evaluate('()=>({frame:openMoonShoreline.frameCount,calls:openMoonShoreline.renderer.info.render.calls,triangles:openMoonShoreline.renderer.info.render.triangles,position:openMoonShoreline.camera.position.toArray(),water:openMoonShoreline.state.ledger})')
            assert result['frame']>0 and result['calls']>0
            page.screenshot(path=str(ROOT/'validation'/f'browser_{name}.png'))
            record['scenes'].append({'name':name,**result})
        page.evaluate('openMoonShoreline.setLocation("shore")')
        before=page.evaluate('openMoonShoreline.camera.position.toArray()')
        page.locator('#world').click(position={'x':630,'y':420});page.keyboard.down('KeyW');page.wait_for_timeout(800);page.keyboard.up('KeyW')
        after=page.evaluate('openMoonShoreline.camera.position.toArray()');assert before!=after
        record['movement']={'before':before,'after':after};record['success']=True;browser.close()
except Exception as e:
    record['errors'].append(str(e))
finally:
    server.shutdown();(ROOT/'validation/full-browser-check.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
assert record['success']
