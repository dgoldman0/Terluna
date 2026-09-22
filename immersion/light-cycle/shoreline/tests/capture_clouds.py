"""Exact-build cloud regression captures; every context-loss message is a failure.
In-memory HTML + harness-only Web Crypto digest bridge where restricted origins
lack it. Browser policy and renderer source remain unchanged. No GPU-speed claim.
"""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, time
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--regime',default='fair',choices=['reference','fog','fair','convection']);ap.add_argument('--world',default='moon');ap.add_argument('--sun',type=float,default=45);ap.add_argument('--width',type=int,default=800);ap.add_argument('--height',type=int,default=500);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--reset',action='store_true');ap.add_argument('--ui',action='store_true');ap.add_argument('--probe',action='store_true');ap.add_argument('--html',type=Path);ap.add_argument('--yaw',type=float,default=0);ap.add_argument('--pitch',type=float,default=.46);ap.add_argument('--cloud-quality',choices=['economy','balanced','high']);ap.add_argument('--ev',type=float,default=0);ap.add_argument('--advance',type=float,default=0);ap.add_argument('--quality',choices=['economy','balanced','high'],default='economy');args=ap.parse_args()
    output=args.output;output.parent.mkdir(parents=True,exist_ok=True)
    html=(args.html or ROOT/'Open_Moon_Shoreline.html').read_bytes();rec={'schema':'open-moon-cloud-capture/1','html_sha256':hashlib.sha256(html).hexdigest(),'regime':args.regime,'world':args.world,'sunElevation_deg':args.sun,'viewport':[args.width,args.height],'errors':[],'warnings':[],'external_requests':[],'checks':{},'harness':'Exact in-memory HTML; SHA-256 bridge; queued application RAF; headed Chromium SwiftShader; synchronous readback. No native WebGPU test.'};start=time.monotonic()
    xvfb=None
    if not Path('/tmp/.X11-unix/X97').exists():
        xvfb=subprocess.Popen(['Xvfb',':97','-screen','0','1600x1000x24','-nolisten','tcp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);time.sleep(.6)
    try:
      with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=False,env={**os.environ,'DISPLAY':':97'},args=['--no-sandbox','--disable-dev-shm-usage','--disable-gpu-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':args.width,'height':args.height},device_scale_factor=1);page.set_default_timeout(180000)
        page.on('pageerror',lambda e:rec['errors'].append(str(e)))
        def console(m):
            if m.type=='error' or any(x in m.text.lower() for x in ['context lost','context_lost','shader error']):rec['errors'].append(m.text)
            elif m.type=='warning':rec['warnings'].append(m.text)
        page.on('console',console)
        page.on('request',lambda r:rec['external_requests'].append(r.url) if not r.url.startswith(('blob:','data:')) else None)
        page.expose_function('cloudDigest',lambda name,values:list(hashlib.new({'SHA-256':'sha256','SHA-1':'sha1'}[name],bytes(values)).digest()))
        page.evaluate("""()=>{if(!crypto.subtle)Object.defineProperty(window,'crypto',{value:{subtle:{digest:async(name,b)=>new Uint8Array(await cloudDigest(name,Array.from(new Uint8Array(b.buffer||b,b.byteOffset||0,b.byteLength)))).buffer}}});window.requestAnimationFrame=cb=>{window.__cloudFrame=cb;return 1;};}""")
        page.set_content(html.decode(),wait_until='domcontentloaded');page.wait_for_function('window.openMoonShoreline?.ready || window.openMoonShorelineError',polling=100)
        error=page.evaluate('window.openMoonShorelineError||null');assert not error,error
        rec['checks']['boot']=True
        page.evaluate("()=>{document.getElementById('world').addEventListener('webglcontextlost',()=>window.__cloudContextLost=true);document.getElementById('unhide').style.display='none';document.getElementById('reticle').style.display='none';}")
        if args.reset:
            before=page.screenshot()
        page.evaluate("""({regime,world,sun,yaw,pitch,quality,cloudQuality,ev})=>{const a=openMoonShoreline;a.setWorld(world);a.quality(quality);if(cloudQuality)a.atmosphere.columnClouds.setQuality(cloudQuality);a.state.ev=ev;a.setPhase(OpenMoonCore.phaseForElevation(sun));if(a.setColumn)a.setColumn(regime);else if(regime!=='reference')throw Error('Source has no column model');a.state.pitch=pitch;a.state.yaw=yaw;a.renderOnce();const g=a.renderer.getContext(),pixel=new Uint8Array(4);g.readPixels(1,1,1,1,g.RGBA,g.UNSIGNED_BYTE,pixel);g.finish();} """,{'regime':args.regime,'world':args.world,'sun':args.sun,'yaw':args.yaw,'pitch':args.pitch,'quality':args.quality,'cloudQuality':args.cloud_quality,'ev':args.ev})
        if args.advance:
            rec['advection']=page.evaluate('''(seconds)=>{const a=openMoonShoreline,c=a.atmosphere.columnClouds,before=c.bakeCount,water=JSON.stringify(a.geography.surfaceWater.ledger);a.state.motionTime+=seconds;a.renderOnce();const g=a.renderer.getContext(),pixel=new Uint8Array(4);g.readPixels(1,1,1,1,g.RGBA,g.UNSIGNED_BYTE,pixel);g.finish();return {seconds,before,after:c.bakeCount,waterUnchanged:water===JSON.stringify(a.geography.surfaceWater.ledger),cacheKey:c.renderKey};}''',args.advance)
            assert rec['advection']['waterUnchanged'];assert rec['advection']['after']>rec['advection']['before'];rec['checks']['wind_advects_without_invented_surface_rain']=True
        rec['state']=page.evaluate('openMoonShoreline.snapshotState()')
        rec['gpu']=page.evaluate("()=>{const g=openMoonShoreline.renderer.getContext(),e=g.getExtension('WEBGL_debug_renderer_info');return {lost:g.isContextLost(),renderer:e?g.getParameter(e.UNMASKED_RENDERER_WEBGL):g.getParameter(g.RENDERER)}}")
        assert not rec['gpu']['lost'];assert not page.evaluate('window.__cloudContextLost||false');assert not rec['errors'],rec['errors']
        page.screenshot(path=str(output));rec['image_sha256']=hashlib.sha256(output.read_bytes()).hexdigest();rec['checks']['scenario_render']=True
        if args.probe:
            rec['fieldProbe']=page.evaluate((ROOT/'tests/cloud_field_probe.js').read_text());assert rec['fieldProbe']['maxAbsoluteError']<.00001;rec['checks']['gpu_column_and_light_fields']=True
        if args.ui:
            page.evaluate("()=>{document.body.classList.remove('minimal');document.getElementById('panel').classList.add('open');}")
            rec['controls']=page.evaluate('''()=>Array.from(document.querySelectorAll('[data-column],#look-sky,#column-export')).map(e=>{const r=e.getBoundingClientRect();return {label:e.textContent,left:r.left,right:r.right,top:r.top,bottom:r.bottom,visible:r.width>0&&r.height>0};})''')
            assert all(v['visible'] and v['left']>=0 and v['right']<=args.width and v['top']>=0 and v['bottom']<=args.height for v in rec['controls']);rec['checks']['column_controls_fit_viewport']=True
            page.screenshot(path=str(output.with_name(output.stem+'-controls.png')))
        if args.reset:
            page.evaluate("()=>{document.body.classList.add('minimal');document.getElementById('panel').classList.remove('open');openMoonShoreline.benchmark();const g=openMoonShoreline.renderer.getContext(),pixel=new Uint8Array(4);g.readPixels(1,1,1,1,g.RGBA,g.UNSIGNED_BYTE,pixel);g.finish();}")
            after=page.screenshot();rec['checks']['pixel_identical_reset']=before==after;assert before==after,'Reset changed reference pixels'
        # Let queued console / context-loss events arrive before deciding success.
        page.wait_for_timeout(300);assert not page.evaluate('window.__cloudContextLost||false');assert not rec['errors'],rec['errors'];assert not rec['external_requests'],rec['external_requests']
        rec['passed']=True;browser.close()
    except Exception as e:rec['passed']=False;rec['errors'].append(type(e).__name__+': '+str(e))
    finally:
        if xvfb:xvfb.terminate()
        rec['elapsed_s']=time.monotonic()-start;output.with_suffix('.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps({'passed':rec['passed'],'errors':rec['errors'],'seconds':rec['elapsed_s'],'state':rec.get('state',{}).get('cloudColumn')}))
    return 0 if rec['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
