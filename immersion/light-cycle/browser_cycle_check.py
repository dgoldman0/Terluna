"""Functional browser checks against the actual offline deliverable.
This container's Chromium disables WebGL; fallback is exercised here.
Exact GPU shaders are independently compiled, linked and rendered with Mesa EGL.
"""
from pathlib import Path
import json,time
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parent
with sync_playwright() as p:
 b=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--disable-dev-shm-usage'])
 page=b.new_page(viewport={'width':1600,'height':1000});errs=[];requests=[]
 page.on('pageerror',lambda e:errs.append(str(e)));page.on('request',lambda r:requests.append(r.url))
 page.set_content((ROOT/'Open_Moon_Full_Cycle.html').read_text(),wait_until='load',timeout=60000)
 page.wait_for_function('window.openMoon?.ready||window.openMoonError',timeout=60000)
 assert page.evaluate('!!window.openMoon?.ready'),page.evaluate('window.openMoonError')
 page.evaluate('openMoon.state.animate=false;document.querySelector("#animate").checked=false;openMoon.draw()')
 ready=page.evaluate('({software:openMoon.software,phase:openMoon.state.phase,angle:openMoon.state.angle,moon:openMoon.illuminance("moon"),earth:openMoon.illuminance("earth"),draws:openMoon.drawCount})');assert ready['angle']==90 and ready['moon']>80000
 page.screenshot(path=str(ROOT/'validation'/'browser_midday.png'))
 page.click('[data-weather="scattered"]');page.evaluate('openMoon.draw()');page.screenshot(path=str(ROOT/'validation'/'browser_clouds.png'))
 page.click('[data-weather="rain"]');page.evaluate('openMoon.draw()');page.screenshot(path=str(ROOT/'validation'/'browser_rain.png'))
 page.click('[data-weather="clear"]');page.click('[data-phase="0.25"]');page.evaluate('openMoon.draw()');assert abs(page.evaluate('openMoon.state.angle'))<1e-10
 page.screenshot(path=str(ROOT/'validation'/'browser_sunset.png'))
 page.click('#compare');page.evaluate('openMoon.draw()');assert page.evaluate('openMoon.state.world')=='compare';page.screenshot(path=str(ROOT/'validation'/'browser_compare.png'))
 page.click('#moonOnly');page.click('[data-phase="0.5"]');page.evaluate('openMoon.draw()');assert page.evaluate('openMoon.state.angle')==-90
 page.select_option('#exposureMode','adaptive');page.wait_for_timeout(2000);page.evaluate('openMoon.draw()');page.screenshot(path=str(ROOT/'validation'/'browser_night_auto.png'))
 page.click('[data-phase="0.75"]');page.evaluate('openMoon.draw()');assert page.evaluate('openMoon.state.side')==-1
 page.click('#east');assert abs(page.evaluate('openMoon.state.yaw')-3.141592653589793)<1e-8
 page.click('#conditions');page.locator('#cover').evaluate('(e)=>{e.value=75;e.dispatchEvent(new Event("input",{bubbles:true}));}');assert page.evaluate('openMoon.state.conditions.cover')==.75
 page.locator('#visibility').evaluate('(e)=>{e.value=-.5;e.dispatchEvent(new Event("input",{bubbles:true}));}');assert page.evaluate('openMoon.state.conditions.visibility')<1000
 page.screenshot(path=str(ROOT/'validation'/'browser_conditions.png'));page.click('#closeWeather')
 page.click('#info');page.select_option('#ozone','moon_no_ozone');assert page.evaluate('openMoon.state.moon')=='moon_no_ozone';page.screenshot(path=str(ROOT/'validation'/'browser_notes.png'));page.click('#closeInfo')
 page.click('#west');page.click('[data-weather="clear"]');page.select_option('#exposureMode','fixed');page.click('[data-phase="0"]');page.select_option('#clock','earth');page.select_option('#speed','60');page.click('#play');t0=time.monotonic();page.wait_for_timeout(1100);page.click('#play');t1=time.monotonic();after=page.evaluate('openMoon.state.phase');assert after>0
 # Loop crossing and exact endpoints via API; playback uses wall-clock increments.
 page.evaluate('openMoon.state.loop=true;openMoon.state.phase=.9999;openMoon.state.speed=259200;openMoon.state.playing=true');page.wait_for_timeout(600);page.evaluate('openMoon.state.playing=false');loop=page.evaluate('openMoon.state.phase');assert 0<=loop<1
 page.select_option('#clock','moon');page.evaluate('openMoon.setPhase(0);openMoon.setView("moon");openMoon.setWeather("clear");openMoon.state.animate=false;openMoon.state.auto=false;openMoon.state.moon="moon";openMoon.draw()')
 page.set_viewport_size({'width':390,'height':844});page.wait_for_timeout(300);page.evaluate('openMoon.draw()');page.screenshot(path=str(ROOT/'validation'/'browser_mobile.png'))
 rects=page.evaluate('Array.from(document.querySelectorAll(".header button,.dock button,.dock select,.viewtools button")).filter(e=>e.getClientRects().length).map(e=>{let r=e.getBoundingClientRect();return{id:e.id,text:e.textContent.trim(),x:r.x,right:r.right,y:r.y,bottom:r.bottom}})')
 outside=[r for r in rects if r['x']<-.5 or r['right']>390.5 or r['y']<-.5 or r['bottom']>844.5]
 record=dict(ready=ready,page_errors=errs,external_requests=[r for r in requests if r.startswith('http')],controls_outside_mobile=outside,earth_clock_playback=dict(wall_seconds=t1-t0,end_phase=after),loop_phase=loop,checked=['offline loading','actual noon at +90 degrees','all cycle quarters','east/west reversal','weather presets','custom cloud and visibility settings','Earth/Moon comparison','ozone setting','camera exposure','clock switching','play/pause','looping','mobile layout'],browser='Chromium page.set_content: software fallback; file navigation disabled by environment policy; GPU renders separately in EGL')
 (ROOT/'validation'/'browser_cycle_checks.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2));assert not errs;assert not outside;b.close()
