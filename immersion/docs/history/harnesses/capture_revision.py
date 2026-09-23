"""Exact-build browser checks. Restricted runner uses in-memory HTML and a
SHA-256 bridge; release code is unchanged. TSL fallback is recorded as WebGL2.
"""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,time
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--lab',action='store_true');ap.add_argument('--case',default='reference',choices=['reference','ponds','reset','walk','benchmark']);ap.add_argument('--width',type=int,default=640);ap.add_argument('--height',type=int,default=400);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    path=ROOT/('Open_Moon_Renderer_Lab.html' if args.lab else 'Open_Moon_Shoreline.html');html=path.read_text();args.output.parent.mkdir(parents=True,exist_ok=True)
    rec={'html_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'lab':args.lab,'case':args.case,'viewport':[args.width,args.height],'errors':[],'console':[],'checks':{},'external_requests':[],'harness':'Headed Chromium / SwiftShader, about:blank in-memory exact source; SHA-256 bridge; queued main-viewer RAF only; normal RAF for node async readback; explicitly forced TSL WebGL2 when --lab'}
    if not Path('/tmp/.X11-unix/X97').exists():
        subprocess.Popen(['Xvfb',':97','-screen','0','1600x1000x24','-nolisten','tcp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);time.sleep(.8)
    start=time.monotonic()
    try:
      with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=False,env={**os.environ,'DISPLAY':':97'},args=['--no-sandbox','--disable-dev-shm-usage','--disable-gpu-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':args.width,'height':args.height},device_scale_factor=1);page.set_default_timeout(180000)
        page.on('pageerror',lambda e:(rec['errors'].append(str(e)),print('ERROR',e,flush=True)))
        page.on('console',lambda m:(rec['console'].append(m.type+': '+m.text),print(m.type,m.text[:1200],flush=True)) if m.type in ['warning','error'] else None)
        page.on('request',lambda r:rec['external_requests'].append(r.url) if not r.url.startswith(('data:','blob:')) else None)
        page.expose_function('testDigest',lambda name,values:list(hashlib.new({'SHA-256':'sha256','SHA-1':'sha1'}[name],bytes(values)).digest()))
        page.evaluate("""()=>{if(!crypto.subtle)Object.defineProperty(window,'crypto',{value:{subtle:{digest:async(name,bytes)=>new Uint8Array(await testDigest(name,Array.from(new Uint8Array(bytes.buffer||bytes,bytes.byteOffset||0,bytes.byteLength)))).buffer}}});window.OM_LAB_BACKEND='webgl';window.__nativeRAF=window.requestAnimationFrame;window.requestAnimationFrame=cb=>{window.__nextFrame=cb;return 1;};}""")
        if args.lab:page.evaluate('()=>{window.requestAnimationFrame=window.__nativeRAF;}')
        page.set_content(html,wait_until='domcontentloaded');print('injected',flush=True)
        page.wait_for_function('window.rendererLab?.ready || window.openMoonShoreline?.ready || window.openMoonShorelineError',polling=100)
        error=page.evaluate('window.openMoonShorelineError||null')
        if error:raise RuntimeError(error)
        name='rendererLab' if args.lab else 'openMoonShoreline';rec['checks']['boot']=True
        page.evaluate("document.body.classList.add('hide-ui');document.getElementById('unhide')?.setAttribute('style','display:none');document.getElementById('reticle')?.setAttribute('style','display:none');")
        if args.lab:
            rec['backend']=page.evaluate('rendererLab.report');rec['state']=page.evaluate('rendererLab.snapshot()')
        else:
            rec['state']=page.evaluate('openMoonShoreline.snapshotState()')
        if args.case=='ponds':
            if args.lab:page.evaluate("async()=>{await rendererLab.setWater(14400);document.getElementById('pools').click();await rendererLab.complete();}")
            else:page.evaluate("()=>{const a=openMoonShoreline;a.setWeather(14400,'episode');a.inspectPonds();a.renderOnce();a.renderer.getContext().finish();}")
            rec['ponds']=page.evaluate(name+'.geography.ponds.model.summary');assert rec['ponds']['reconstructedVolume_m3']>0;assert abs(rec['ponds']['accountingResidual_m3'])<1e-7;rec['checks']['pond_geometry_volume']=True
        if args.case=='reset' and not args.lab:
            ref=page.screenshot();page.evaluate("()=>{const a=openMoonShoreline;a.setWeather(14400,'episode');a.inspectPonds();a.renderOnce();a.benchmark();a.renderer.getContext().finish();}");after=page.screenshot();rec['resetPixelIdentical']=ref==after;assert ref==after;rec['checks']['reset_identical']=True
        if args.case=='benchmark':
            if args.lab:rec['route']=page.evaluate('rendererLab.benchmark(12)')
            else:rec['route']=page.evaluate("""()=>{const a=openMoonShoreline,times=[],cpu=[],draws=[],pixel=new Uint8Array(4),g=a.renderer.getContext();const complete=()=>{g.readPixels(1,1,1,1,g.RGBA,g.UNSIGNED_BYTE,pixel);g.finish();};a.renderer.info.autoReset=false;a.benchmark();for(let i=0;i<3;i++){a.renderer.info.reset();a.renderOnce();complete();}const state=a.snapshotState();for(let i=0;i<12;i++){const x=i*.15,z=9+i*.025;a.camera.position.set(x,a.geography.height(x,z)+1.7,z);a.renderer.info.reset();const start=performance.now();a.renderOnce();cpu.push(performance.now()-start);complete();times.push(performance.now()-start);draws.push({...a.renderer.info.render});if(a.renderer.getContext().isContextLost())throw Error('Context lost on benchmark frame '+i);}const quant=(a,q)=>[...a].sort((a,b)=>a-b)[Math.min(a.length-1,Math.floor(a.length*q))];return {frames:12,state,median_ms:quant(times,.5),p95_ms:quant(times,.95),max_ms:Math.max(...times),submissionMedian_ms:quant(cpu,.5),times_ms:times,draws,memory:{...a.renderer.info.memory},completionPixel:Array.from(pixel),method:'Synchronous one-pixel readback plus finish; identical prescribed 12-frame route'};}""")
            rec['checks']['benchmark_completed']=True
        if args.case=='walk':
            if args.lab:rec['route']=page.evaluate('rendererLab.benchmark(12)')
            else:
                rec['route']=page.evaluate("""()=>{const a=openMoonShoreline,times=[],cpu=[];for(let i=0;i<24;i++){a.state.keys={KeyD:true};const start=performance.now();window.__nextFrame(a.state.last+60);cpu.push(performance.now()-start);a.renderer.getContext().finish();times.push(performance.now()-start);if(a.renderer.getContext().isContextLost())throw Error('Context lost on frame '+i);}a.state.keys={};return {frames:24,position:a.camera.position.toArray(),times_ms:times,submission_ms:cpu,method:'real walking callbacks, synchronized GPU completion'};}""")
            rec['checks']['walk_completed']=True
        if args.lab and args.case=='ponds':
            rec['fieldReadback']=page.evaluate((ROOT/'tests/node_field_probe.js').read_text());assert rec['fieldReadback']['maxAbsoluteError']<.001;rec['checks']['node_fields_reach_gpu']=True
        if not args.lab:
            rec['gpu']=page.evaluate("()=>{const g=openMoonShoreline.renderer.getContext(),e=g.getExtension('WEBGL_debug_renderer_info');return {lost:g.isContextLost(),renderer:g.getParameter(e.UNMASKED_RENDERER_WEBGL)}}")
            assert not rec['gpu']['lost']
        await_done=page.evaluate(name+'.complete()') if args.lab else page.evaluate(name+'.renderer.getContext().finish()')
        rec['finalState']=page.evaluate(name+('.snapshot()' if args.lab else '.snapshotState()'));
        page.screenshot(path=str(args.output));rec['checks']['capture']=True;rec['passed']=not rec['errors'];browser.close()
    except Exception as e:
        rec['errors'].append(type(e).__name__+': '+str(e));rec['passed']=False
    rec['elapsed_s']=time.monotonic()-start;args.output.with_suffix('.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2),flush=True);return 0 if rec['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
