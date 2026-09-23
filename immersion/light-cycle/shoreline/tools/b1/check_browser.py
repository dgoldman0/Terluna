"""Exercise the exact B1 offline build and inspect the actual canvas pixels.

No browser security policy is modified. If file navigation is disallowed by
this environment, render the identical document with set_content and record it.
"""
from __future__ import annotations
import hashlib,json,time
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'validation/b1'
PIXELS=r'''id=>{const c=document.getElementById(id),x=c.getContext('2d'),a=x.getImageData(0,0,c.width,c.height).data;let lo=255,hi=0,s=0,ss=0,n=0,white=0,opaque=0,hash=2166136261;const step=Math.max(1,Math.floor(Math.sqrt(c.width*c.height/8192)));for(let y=0;y<c.height;y+=step)for(let z=0;z<c.width;z+=step){const i=4*(y*c.width+z),v=.2126*a[i]+.7152*a[i+1]+.0722*a[i+2];lo=Math.min(lo,v);hi=Math.max(hi,v);s+=v;ss+=v*v;n++;if(a[i+3]===255)opaque++;if(a[i]>=254&&a[i+1]>=254&&a[i+2]>=254)white++;for(let k=0;k<4;k++){hash^=a[i+k];hash=Math.imul(hash,16777619);}}return{width:c.width,height:c.height,css_width:c.getBoundingClientRect().width,sampled_pixels:n,opaque_fraction:opaque/n,luminance_min:lo,luminance_max:hi,luminance_range:hi-lo,standard_deviation:Math.sqrt(Math.max(0,ss/n-(s/n)**2)),white_fraction:white/n,pixel_hash:hash>>>0};}'''

def main():
 OUT.mkdir(parents=True,exist_ok=True);t0=time.monotonic();f=ROOT/'Open_Moon_B1_Scene_Painter.html';html=f.read_text()
 report=dict(schema='open-moon-b1-browser/1',html_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),errors=[],console_errors=[],external_requests=[],navigation_attempts=[],layouts=[],state_checks=0,passed=False)
 try:
  with sync_playwright() as pw:
   browser=pw.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--disable-dev-shm-usage'])
   for name,w,h,dpr in [('desktop',1280,1000,1),('desktop-hidpi',1280,1000,2),('mobile',390,844,2)]:
    context=browser.new_context(viewport={'width':w,'height':h},device_scale_factor=dpr,accept_downloads=True)
    if name=='desktop':
     probe=context.new_page()
     try:
      probe.goto(f.as_uri(),wait_until='load',timeout=10000)
      probe.wait_for_function('window.B1Data?.frames?.length>0',timeout=3000)
      report['navigation_attempts'].append(dict(method='file',success=True))
     except Exception as e:report['navigation_attempts'].append(dict(method='file',success=False,error=str(e)[:1500]))
     probe.close()
    page=context.new_page()
    page.on('pageerror',lambda e:report['errors'].append(str(e)))
    page.on('console',lambda m:report['console_errors'].append(m.text) if m.type=='error' else None)
    page.on('request',lambda r:report['external_requests'].append(r.url) if not r.url.startswith(('data:','blob:','file:','about:')) else None)
    page.set_content(html,wait_until='load');page.wait_for_function('window.B1Data?.frames?.length===18')
    layout=dict(name=name,device_pixel_ratio=dpr,cases=[],checks=0)
    cases=page.evaluate('B1Data.frames.map(f=>({scenario:f.meta.scenario,sun:f.meta.sun_elevation_deg}))')
    for case in cases:
     page.select_option('#scenario',case['scenario']);page.select_option('#sun',str(case['sun']))
     assert page.evaluate('B1CurrentMetadata.scenario')==case['scenario']
     assert page.evaluate('B1CurrentMetadata.sun_elevation_deg')==case['sun']
     checks=[]
     for filtered in [False,True]:
      page.locator('#filtered').set_checked(filtered)
      im=page.evaluate(PIXELS,'image');sky=page.evaluate(PIXELS,'dome')
      assert im['opaque_fraction']==1 and sky['opaque_fraction']==1
      assert im['luminance_range']>3 and im['standard_deviation']>1
      assert im['white_fraction']<.1 and sky['white_fraction']<.1
      for pix in [im,sky]:assert abs(pix['width']-round(pix['css_width']*dpr))<=1
      checks.append(dict(filtered=filtered,image=im,sky=sky));report['state_checks']+=1;layout['checks']+=1
     grounds=[]
     for mode in ['direct','diffuse','total','shadow']:
      page.select_option('#mapMode',mode);px=page.evaluate(PIXELS,'ground')
      assert px['opaque_fraction']==1 and abs(px['width']-round(px['css_width']*dpr))<=1
      # A visibility map may legitimately be all zero after sunset. Never
      # demand arbitrary contrast in a physically dark direct-light product.
      grounds.append(dict(mode=mode,**px));report['state_checks']+=1;layout['checks']+=1
     assert page.evaluate('document.documentElement.scrollWidth-innerWidth')<=0
     layout['cases'].append(dict(**case,views=checks,ground=grounds))
    page.select_option('#scenario','fair');page.select_option('#sun','12');page.select_option('#mapMode','direct');page.locator('#filtered').check()
    before=page.evaluate('JSON.stringify(B1CurrentMetadata)+B1Data.frames.map(f=>f.raw).join("")')
    px0=page.evaluate(PIXELS,'image')
    page.locator('#exposure').evaluate('e=>{e.value="2";e.dispatchEvent(new Event("input",{bubbles:true}))}')
    px1=page.evaluate(PIXELS,'image')
    assert px0['pixel_hash']!=px1['pixel_hash']
    assert before==page.evaluate('JSON.stringify(B1CurrentMetadata)+B1Data.frames.map(f=>f.raw).join("")')
    page.click('#reset');assert page.locator('#exposure').input_value()=='0'
    assert not page.locator('#blackControl').is_disabled()
    page.locator('#blackControl').check();black=page.evaluate(PIXELS,'image')
    assert black['pixel_hash']!=px0['pixel_hash']
    assert page.evaluate('B1CurrentMetadata.surface_boundary')=='black control'
    if name=='desktop':
     with page.expect_download() as event:page.click('#download')
     dest=OUT/'metadata-download.json';event.value.save_as(dest)
     expected=page.evaluate('JSON.stringify(B1CurrentMetadata,null,2)+"\\n"')
     assert dest.read_text()==expected
     layout['download_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();dest.unlink()
    page.locator('#blackControl').uncheck()
    page.click('#play');sun0=page.locator('#sun').input_value();page.wait_for_timeout(1450)
    assert page.locator('#sun').input_value()!=sun0 and page.locator('#exposure').input_value()=='0'
    page.click('#play');page.select_option('#sun','12');page.select_option('#mapMode','direct')
    page.evaluate('scrollTo(0,0)');page.screenshot(path=str(OUT/('lab-'+name+'.png')),full_page=True)
    layout['exposure_mutates_data']=False;layout['black_control']=black;layout['overflow_px']=page.evaluate('document.documentElement.scrollWidth-innerWidth');layout['default_pixels']=page.evaluate(PIXELS,'image')
    report['state_checks']+=5;layout['checks']+=5;report['layouts'].append(layout)
    context.close()
   browser.close()
  report['passed']=not any(report[k] for k in ['errors','console_errors','external_requests'])
 except Exception as e:report['errors'].append(type(e).__name__+': '+str(e))
 report['load_method']='Exact standalone HTML inserted using set_content; local navigation result separately recorded; no policy changes.'
 report['elapsed_s']=time.monotonic()-t0
 (OUT/'browser.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='layouts'},indent=2))
 return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
