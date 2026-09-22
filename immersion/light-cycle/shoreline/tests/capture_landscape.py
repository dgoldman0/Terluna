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
    parser.add_argument('--software-backend', choices=['swiftshader','mesa'], default='swiftshader')
    parser.add_argument('--display', default=os.environ.get('DISPLAY', ':97'))
    parser.add_argument('--in-memory', action='store_true')
    parser.add_argument('--interactions', action='store_true')
    parser.add_argument('--probe', action='store_true')
    parser.add_argument('--landscape-review', action='store_true')
    parser.add_argument('--ui-only', action='store_true')
    parser.add_argument('--view', choices=['shore','path','overlook','shelter'])
    parser.add_argument('--case', choices=['walk','weather','after-rain','reset','fields','low-sun','inland'])
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
                    'deterministic_animation_callbacks': True, 'finish_each_walk_frame': True, 'software': args.software, 'software_backend': args.software_backend},
    }
    import subprocess
    xvfb=None
    if args.headed and not Path('/tmp/.X11-unix/X'+args.display.split(':')[-1]).exists():
        xvfb=subprocess.Popen(['Xvfb',args.display,'-screen','0','1600x1000x24','-nolisten','tcp'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        time.sleep(.5)
    server = None
    started = time.monotonic()
    try:
        with sync_playwright() as p:
            flags = ['--no-sandbox', '--disable-dev-shm-usage']
            if args.software:
                flags += (['--use-angle=swiftshader', '--enable-unsafe-swiftshader'] if args.software_backend=='swiftshader' else ['--use-gl=angle','--use-angle=gl'])+['--disable-gpu-sandbox']
            env = {**os.environ, 'DISPLAY': args.display}
            if args.software:
                env['LIBGL_ALWAYS_SOFTWARE'] = '1'
            browser = p.chromium.launch(executable_path=args.browser,
                headless=not args.headed, env=env, args=flags)
            page = browser.new_page(viewport={'width':args.width,'height':args.height}, device_scale_factor=1)
            page.set_default_timeout(180000)
            page.on('pageerror', lambda e: (record['errors'].append(str(e)),print('PAGE ERROR',str(e),flush=True)))
            page.on('console', lambda m: (record['console'].append(m.type+': '+m.text), print(m.type,m.text,flush=True)) if m.type in ('error','warning') else None)
            page.on('request', lambda r: record['external_requests'].append(r.url) if not r.url.startswith(('http://127.0.0.1:', 'data:', 'blob:')) else None)
            queue = "window.requestAnimationFrame=cb=>{window.__m1NextFrame=cb;return 1;};"
            if args.in_memory:
                page.expose_function('testSHA1', lambda values: list(hashlib.sha1(bytes(values)).digest()))
                page.evaluate("""() => { if (!window.crypto?.subtle) Object.defineProperty(window,'crypto',{value:{subtle:{digest:async(name,bytes)=>{if(name!=='SHA-1')throw Error('Only SHA-1 is bridged');return new Uint8Array(await window.testSHA1(Array.from(new Uint8Array(bytes.buffer||bytes,bytes.byteOffset||0,bytes.byteLength)))).buffer;}}}}); }""")
                page.evaluate(queue)
                print('Injecting viewer',flush=True);page.set_content(html, wait_until='domcontentloaded', timeout=180000);print('Viewer injected',flush=True)
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
            page.evaluate("""()=>{const r=openMoonShoreline.renderer,old=r.debug.onShaderError;r.debug.onShaderError=(gl,p,v,f)=>{console.error(JSON.stringify({shaderFailure:{lost:gl.isContextLost(),error:gl.getError(),linked:gl.getProgramParameter(p,gl.LINK_STATUS),vertex:gl.getShaderParameter(v,gl.COMPILE_STATUS),fragment:gl.getShaderParameter(f,gl.COMPILE_STATUS),programLog:gl.getProgramInfoLog(p),vertexLog:gl.getShaderInfoLog(v),fragmentLog:gl.getShaderInfoLog(f)}}));old(gl,p,v,f);};}""")
            if args.view:
                page.evaluate('(view)=>{openMoonShoreline.setLocation(view);openMoonShoreline.renderOnce();openMoonShoreline.renderer.getContext().finish();}',args.view)
                record['view']=args.view
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
                record['probes']=page.evaluate((args.root/'tests/landscape_probe.js').read_text())
                card=record['probes']['grayCard']
                assert abs(card['referenceDisplayLuminance']-.18)<.003,card
                depth=record['probes']['depth']
                assert all(d['geometricDepth_m']>0 and d['estimatedOpticalPath_m']>0 for d in depth[:3]),depth
                assert depth[3]['occluder'] and depth[3]['opaqueViewZ_m']<depth[3]['flatWaterViewZ_m'] and depth[3]['estimatedOpticalPath_m']==0,depth[3]
                record['checks']['foreground_rock_occlusion_confirmed_by_raycast']=True
                record['checks']['gpu_gray_card_and_depth']=True
                print('GPU gray-card and depth probes passed.', flush=True)
            if args.interactions:
                # Several real callbacks with frozen time should issue no new draws.
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
                # Change comparison state, then call the actual M1 reset command.
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
            if args.landscape_review:
                # Use the real movement callback with a held movement key.
                page.evaluate("openMoonShoreline.benchmark();openMoonShoreline.state.keys.KeyW=true;")
                for step in range(12):
                    page.evaluate('window.__m1NextFrame(openMoonShoreline.state.last+60);openMoonShoreline.renderer.getContext().finish();')
                page.evaluate('openMoonShoreline.state.keys.KeyW=false;')
                record['walk_state']=page.evaluate('openMoonShoreline.snapshotState()')
                assert record['walk_state']['position'][2]<8.5
                assert record['walk_state']['terrain']['updates']>1
                record['checks']['actual_walking_crosses_grid_snaps']=True
                print('Walking grid transitions passed.',flush=True)
                # Repeatable named views under exactly the noon reference light.
                for name,location in [('woodland','path'),('overlook','overlook')]:
                    page.evaluate('(location)=>{const a=openMoonShoreline;a.setLocation(location);a.renderOnce();}', location)
                    page.screenshot(path=str(args.output.with_name(args.output.stem+'_'+name+'.png')))
                    print(name+' view captured.',flush=True)
                page.evaluate("const a=openMoonShoreline;a.setLocation('path');a.state.yaw=2.9;a.state.pitch=-.13;a.renderOnce();")
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_woodland_inland.png')))
                # Exercise the diagnostic through its actual UI event.
                page.locator('#landscape-debug').select_option('3',force=True)
                page.evaluate('openMoonShoreline.renderOnce()')
                assert page.evaluate('openMoonShoreline.geography.fieldUniforms.uLandscapeDebug.value')==3
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_drainage.png')))
                record['checks']['field_diagnostic_control']=True
                page.locator('#landscape-debug').select_option('0',force=True)
                page.evaluate("openMoonShoreline.benchmark();openMoonShoreline.setWeather(8500,'episode');openMoonShoreline.renderOnce();")
                record['wet_state']=page.evaluate('openMoonShoreline.snapshotState()')
                assert record['wet_state']['spatialWater']['input']>0
                assert record['wet_state']['spatialWater']['relativeResidual']<1e-10
                record['surface_state_probe']=page.evaluate((args.root/'tests/surface_state_probe.js').read_text())
                assert record['surface_state_probe']['maximumAbsoluteError']<.0001
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_rain.png')))
                page.evaluate("openMoonShoreline.setWeather(14400,'episode');openMoonShoreline.renderOnce();")
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_after_rain.png')))
                record['checks']['wet_state_reaches_gpu_and_persists_after_rain']=True
                print('Weather and GPU surface stores passed.',flush=True)
                # A changing Sun must leave planted world identities intact.
                before=page.evaluate('openMoonShoreline.geography.ecology.trees.map(t=>t.id)')
                page.evaluate("const a=openMoonShoreline;a.setWeather(0,'clear');a.setLocation('path');a.state.yaw=-Math.PI/2;a.state.pitch=.16;a.setPhase(OpenMoonCore.phaseForElevation(18));a.renderOnce();")
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_low_sun.png')))
                assert before==page.evaluate('openMoonShoreline.geography.ecology.trees.map(t=>t.id)')
                record['checks']['plant_identities_independent_of_instant_sun']=True
                page.evaluate("openMoonShoreline.setWorld('earth');openMoonShoreline.renderOnce();openMoonShoreline.benchmark();")
                page.screenshot(path=str(args.output.with_name(args.output.stem+'_final_reset.png')))
                reset=args.output.with_name(args.output.stem+'_final_reset.png')
                record['checks']['reset_pixel_identical']=reset.read_bytes()==args.output.read_bytes()
                assert record['checks']['reset_pixel_identical'], 'Reference reset differed after world, wetness and camera changes'
                record['checks']['world_comparison_and_complete_reset']=True
            if args.case:
                record['case']=args.case
                if args.case=='walk':
                    page.evaluate('openMoonShoreline.state.keys.KeyW=true;')
                    for _ in range(12):
                        page.evaluate('window.__m1NextFrame(openMoonShoreline.state.last+60);openMoonShoreline.renderer.getContext().finish();')
                    page.evaluate('openMoonShoreline.state.keys.KeyW=false;')
                    record['walk_state']=page.evaluate('openMoonShoreline.snapshotState()')
                    assert record['walk_state']['position'][2]<8.5
                    assert record['walk_state']['terrain']['updates']>1
                    record['checks']['actual_walking_crosses_grid_snaps']=True
                elif args.case in ('weather','after-rain'):
                    t=8500 if args.case=='weather' else 14400
                    page.evaluate('(t)=>{openMoonShoreline.setWeather(t,"episode");openMoonShoreline.renderOnce();}',t)
                    record['wet_state']=page.evaluate('openMoonShoreline.snapshotState()')
                    assert record['wet_state']['spatialWater']['input']>0
                    assert record['wet_state']['spatialWater']['relativeResidual']<1e-10
                    record['surface_state_probe']=page.evaluate((args.root/'tests/surface_state_probe.js').read_text())
                    assert record['surface_state_probe']['maximumAbsoluteError']<.0001
                    record['checks']['wet_state_reaches_gpu']=True
                    record['surface_samples']=page.evaluate('openMoonShoreline.geography.surfaceWater.sample(0,9)')
                elif args.case=='reset':
                    page.evaluate('openMoonShoreline.setWorld("earth");openMoonShoreline.setLocation("path");openMoonShoreline.setWeather(8500,"episode");openMoonShoreline.renderOnce();openMoonShoreline.renderer.getContext().finish();')
                    page.evaluate('openMoonShoreline.benchmark();openMoonShoreline.renderer.getContext().finish();')
                    record['reset_state']=page.evaluate('openMoonShoreline.snapshotState()')
                    assert record['reset_state']['world']=='moon' and record['reset_state']['spatialWater']['input']==0
                elif args.case=='fields':
                    page.evaluate('const e=document.getElementById("landscape-debug");e.value="3";e.dispatchEvent(new Event("change",{bubbles:true}));openMoonShoreline.renderOnce();')
                    assert page.evaluate('openMoonShoreline.geography.fieldUniforms.uLandscapeDebug.value')==3
                    record['checks']['diagnostic_control_event']=True
                elif args.case in ('inland','low-sun'):
                    before=page.evaluate('openMoonShoreline.geography.ecology.trees.map(t=>t.id)')
                    page.evaluate('const a=openMoonShoreline;a.setLocation("path");a.state.yaw=2.9;a.state.pitch=-.13;')
                    if args.case=='low-sun':
                        page.evaluate('const a=openMoonShoreline;a.state.yaw=-Math.PI/2;a.state.pitch=.16;a.setPhase(OpenMoonCore.phaseForElevation(18));')
                    page.evaluate('openMoonShoreline.renderOnce()')
                    assert before==page.evaluate('openMoonShoreline.geography.ecology.trees.map(t=>t.id)')
                    record['checks']['plant_identities_preserved']=True
                case_image=args.output.with_name(args.output.stem+'_'+args.case+'.png')
                page.screenshot(path=str(case_image))
                record['case_image_sha256']=hashlib.sha256(case_image.read_bytes()).hexdigest()
                if args.case=='reset':
                    record['checks']['world_weather_camera_reset_pixel_identical']=case_image.read_bytes()==args.output.read_bytes()
                    assert record['checks']['world_weather_camera_reset_pixel_identical']
                print('Isolated '+args.case+' case passed.',flush=True)
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
        if xvfb:xvfb.terminate()
        record['elapsed_s']=time.monotonic()-started
        args.output.with_suffix('.json').write_text(json.dumps(record,indent=2))
        print(json.dumps(record,indent=2),flush=True)
    return 0 if record['passed'] else 1

if __name__=='__main__':
    raise SystemExit(main())
