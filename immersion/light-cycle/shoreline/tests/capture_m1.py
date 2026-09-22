"""Render the actual bundled Three.js application and capture repeatable evidence.

Normal mode serves the build on localhost and uses browser Web Crypto. A restricted
runner may use --in-memory to inject local HTML into about:blank, where the missing
Web Crypto SHA-1 operation is bridged to Python hashlib. The vendor bytes and the
application's expected Git-blob hash are identical in both modes. The bridge is a
harness facility, never included in the distributed viewer.

For deterministic stills the harness queues requestAnimationFrame callbacks and
advances the real application callback explicitly. Rendering/shaders stay intact.
"""
from __future__ import annotations
import argparse
import functools
import hashlib
import json
import os
from pathlib import Path
import threading
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from playwright.sync_api import sync_playwright


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--width', type=int, default=1440)
    parser.add_argument('--height', type=int, default=900)
    parser.add_argument('--browser', default='/usr/bin/chromium')
    parser.add_argument('--headed', action='store_true')
    parser.add_argument('--software', action='store_true')
    parser.add_argument('--display', default=os.environ.get('DISPLAY', ':97'))
    parser.add_argument('--in-memory', action='store_true')
    parser.add_argument('--interactions', action='store_true')
    parser.add_argument('--probe', action='store_true')
    parser.add_argument('--ui-only', action='store_true')
    args = parser.parse_args()
    if min(args.width, args.height) < 128:
        parser.error('Viewport dimensions must be at least 128 pixels.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    html_path = args.root / 'Open_Moon_Shoreline.html'
    html = html_path.read_text()
    record = {
        'html_sha256': hashlib.sha256(html_path.read_bytes()).hexdigest(),
        'viewport': [args.width, args.height], 'errors': [], 'console': [],
        'external_requests': [], 'checks': {},
        'harness': {'in_memory': args.in_memory, 'sha1_bridge': args.in_memory,
                    'deterministic_animation_callbacks': True, 'software': args.software},
    }
    server = None
    started = time.monotonic()
    try:
        with sync_playwright() as p:
            flags = ['--no-sandbox', '--disable-dev-shm-usage']
            if args.software:
                flags += ['--use-angle=swiftshader', '--enable-unsafe-swiftshader']
            env = {**os.environ, 'DISPLAY': args.display}
            if args.software:
                env['LIBGL_ALWAYS_SOFTWARE'] = '1'
            browser = p.chromium.launch(executable_path=args.browser,
                headless=not args.headed, env=env, args=flags)
            page = browser.new_page(viewport={'width':args.width,'height':args.height}, device_scale_factor=1)
            page.set_default_timeout(180000)
            page.on('pageerror', lambda e: record['errors'].append(str(e)))
            page.on('console', lambda m: record['console'].append(m.type+': '+m.text) if m.type in ('error','warning') else None)
            page.on('request', lambda r: record['external_requests'].append(r.url) if not r.url.startswith(('http://127.0.0.1:', 'data:', 'blob:')) else None)
            queue = "window.requestAnimationFrame=cb=>{window.__m1NextFrame=cb;return 1;};"
            if args.in_memory:
                page.expose_function('testSHA1', lambda values: list(hashlib.sha1(bytes(values)).digest()))
                page.evaluate("""() => { if (!window.crypto?.subtle) Object.defineProperty(window,'crypto',{value:{subtle:{digest:async(name,bytes)=>{if(name!=='SHA-1')throw Error('Only SHA-1 is bridged');return new Uint8Array(await window.testSHA1(Array.from(new Uint8Array(bytes.buffer||bytes,bytes.byteOffset||0,bytes.byteLength)))).buffer;}}}}); }""")
                page.evaluate(queue)
                page.set_content(html, wait_until='domcontentloaded', timeout=180000)
            else:
                class QuietHandler(SimpleHTTPRequestHandler):
                    def log_message(self, *unused): pass
                server=ThreadingHTTPServer(('127.0.0.1',0),functools.partial(QuietHandler,directory=str(args.root)))
                threading.Thread(target=server.serve_forever,daemon=True).start()
                page.add_init_script(queue)
                page.goto(f'http://127.0.0.1:{server.server_port}/Open_Moon_Shoreline.html',wait_until='domcontentloaded')
            page.wait_for_function('window.openMoonShoreline?.ready || window.openMoonShorelineError', polling=100, timeout=180000)
            error=page.evaluate('window.openMoonShorelineError || null')
            if error: raise RuntimeError(error)
            record['state']=page.evaluate('openMoonShoreline.snapshotState()')
            record['gpu']=page.evaluate("""()=>{const gl=openMoonShoreline.renderer.getContext(),e=gl.getExtension('WEBGL_debug_renderer_info');return {version:gl.getParameter(gl.VERSION),renderer:e?gl.getParameter(e.UNMASKED_RENDERER_WEBGL):gl.getParameter(gl.RENDERER),maxTextureSize:gl.getParameter(gl.MAX_TEXTURE_SIZE)};}""")
            record['counts']=page.evaluate('openMoonShoreline.geography.counts')
            record['initial_environment_updates']=page.evaluate('openMoonShoreline.atmosphere.environmentUpdates || 0')
            record['checks']['boot']=True
            page.evaluate("document.getElementById('unhide').style.visibility='hidden';document.getElementById('reticle').style.visibility='hidden';")
            page.screenshot(path=str(args.output), timeout=180000)
            record['checks']['reference_screenshot']=True
            record['screenshot_sha256']=hashlib.sha256(args.output.read_bytes()).hexdigest()
            print('Reference screenshot captured.', flush=True)
            if args.ui_only:
                page.evaluate('window.__m1NextFrame(performance.now()+300)')
                page.evaluate("document.getElementById('unhide').style.visibility='visible';")
                page.locator('#unhide').click()
                assert not page.evaluate("document.body.classList.contains('minimal')")
                record['ui_bounds']=page.evaluate("""()=>{return ['m1-benchmark','help-toggle','panel-toggle'].map(id=>{const r=document.getElementById(id).getBoundingClientRect();return {id,x:r.x,y:r.y,right:r.right,bottom:r.bottom,visible:r.width>0&&r.height>0};});}""")
                assert all(b['visible'] and b['x']>=0 and b['right']<=args.width for b in record['ui_bounds']),record['ui_bounds']
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_controls.png')))
                record['checks']['controls_fit_viewport']=True
            if args.probe:
                record['probes']=page.evaluate((args.root/'tests/m1_probe.js').read_text())
                card=record['probes']['grayCard']
                assert abs(card['referenceDisplayLuminance']-.18)<.003,card
                assert all(d['geometricDepth_m']>0 and d['estimatedOpticalPath_m']>0 for d in record['probes']['depth'])
                record['checks']['gpu_gray_card_and_depth']=True
                print('GPU gray-card and depth probes passed.', flush=True)
            if args.interactions:
                frozen=page.evaluate("""()=>{const a=openMoonShoreline,before=a.water.passes;for(let i=0;i<5;i++)window.__m1NextFrame(performance.now()+20*i);return {before,after:a.water.passes,time:a.state.motionTime};}""")
                assert frozen['before']==frozen['after'] and frozen['time']==0,frozen
                record['checks']['paused_callbacks_skip_render']=True
                page.evaluate("document.getElementById('unhide').style.visibility='visible';")
                page.locator('#unhide').click()
                assert not page.evaluate("document.body.classList.contains('minimal')")
                record['checks']['show_controls']=True
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_controls.png')))
                page.locator('#shore-edge').click()
                page.evaluate('openMoonShoreline.renderOnce()')
                assert not page.evaluate('openMoonShoreline.renderer.getContext().isContextLost()')
                edge=page.evaluate('openMoonShoreline.snapshotState()')
                assert edge['position'][2]==-1,edge
                record['edge_state']=edge
                record['checks']['second_view_with_depth_and_reflection']=True
                page.evaluate("document.body.classList.add('minimal');document.getElementById('unhide').style.visibility='hidden';")
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_water_edge.png')))
                page.evaluate("openMoonShoreline.atmosphere.whiteBalanceStrength=0;openMoonShoreline.state.ev=2;openMoonShoreline.benchmark();")
                restored=page.evaluate('openMoonShoreline.snapshotState()')
                assert restored['whiteBalanceStrength']==1 and restored['ev']==0 and restored['phase']==0 and restored['freeze']
                assert restored['world']=='moon' and restored['weatherKind']=='clear' and restored['position'][2]==9
                record['checks']['benchmark_reset']=True
                assert not page.evaluate('openMoonShoreline.renderer.getContext().isContextLost()')
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_reset.png')))
                assert not page.evaluate('openMoonShoreline.renderer.getContext().isContextLost()'), 'Context lost after reset capture'
                record['checks']['multiple_full_render_passes']=True
                print('Repeated rendering and control checks passed.', flush=True)
            record['final_context_lost']=page.evaluate('openMoonShoreline.renderer.getContext().isContextLost()')
            assert not record['final_context_lost']
            record['final_environment_updates']=page.evaluate('openMoonShoreline.atmosphere.environmentUpdates || 0')
            record['final_render_info']=page.evaluate('({...openMoonShoreline.renderer.info.render})')
            browser.close()
        assert not record['errors'] and not record['external_requests'],record
        record['passed']=True
    except Exception as exc:
        record['errors'].append(type(exc).__name__+': '+str(exc)); record['passed']=False
    finally:
        if server: server.shutdown()
        record['elapsed_s']=time.monotonic()-started
        args.output.with_suffix('.json').write_text(json.dumps(record,indent=2))
        print(json.dumps(record,indent=2),flush=True)
    return 0 if record['passed'] else 1

if __name__=='__main__':
    raise SystemExit(main())
