"""Exercise exact offline A3 and corrected A2 builds, including real canvas pixels."""
from __future__ import annotations
import hashlib, json, time
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'validation/a3'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
PIXELS=r'''(id)=>{const c=document.getElementById(id),d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;
 let min=255,max=0,white=0,opaque=0,colour=0;for(let i=0;i<d.length;i+=4){if(d[i+3]<128)continue;opaque++;const y=.2126*d[i]+.7152*d[i+1]+.0722*d[i+2];min=Math.min(min,y);max=Math.max(max,y);if(d[i]>=254&&d[i+1]>=254&&d[i+2]>=254)white++;if(Math.max(d[i],d[i+1],d[i+2])-Math.min(d[i],d[i+1],d[i+2])>25)colour++;}
 return {width:c.width,height:c.height,css_width:c.getBoundingClientRect().width,opaque_pixels:opaque,luminance_min:min,luminance_max:max,luminance_range:max-min,white_fraction:white/Math.max(1,opaque),chromatic_pixels:colour};}'''
def main():
 OUT.mkdir(parents=True,exist_ok=True);start=time.monotonic()
 report=dict(schema='open-moon-a3-browser/1',html_sha256=sha(ROOT/'Open_Moon_A3_Accuracy_Lab.html'),a2_html_sha256=sha(ROOT/'Open_Moon_A2_Accuracy_Lab.html'),errors=[],console_errors=[],external_requests=[],states_checked=0,layouts=[],a2_pixel_checks=[],load_method='Exact HTML bytes inserted with set_content; file navigation is blocked by browser policy. No policy settings changed.',passed=False)
 try:
  with sync_playwright() as p:
   browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--disable-dev-shm-usage'])
   for label,w,h,dpr in [('desktop',1280,980,1),('desktop-hidpi',1280,980,2),('mobile',390,844,2)]:
    context=browser.new_context(viewport={'width':w,'height':h},device_scale_factor=dpr,accept_downloads=True)
    page=context.new_page()
    page.on('pageerror',lambda e:report['errors'].append(str(e)))
    page.on('console',lambda m:report['console_errors'].append(m.text) if m.type=='error' else None)
    page.on('request',lambda r:report['external_requests'].append(r.url) if not r.url.startswith(('file:','data:','blob:')) else None)
    page.set_content((ROOT/'Open_Moon_A3_Accuracy_Lab.html').read_text(),wait_until='load')
    page.wait_for_function('window.a3AccuracyLab?.ready === true')
    assert page.locator('#ready').inner_text()=='Offline · matched A3 benchmark'
    assert 'FAIL' in page.locator('#gateRows').inner_text(), 'Inherited A2 failure must remain visible'
    layout=dict(name=label,device_pixel_ratio=dpr,pixel_checks=[])
    for ci in range(8):
     page.select_option('#case',str(ci))
     for comp in ['coupling','closure','orders']:
      page.select_option('#comparison',comp)
      assert 'NaN' not in page.locator('#treatmentRows').inner_text()
      assert page.locator('#treatmentRows tr').count()==7
      spectrum=page.evaluate(PIXELS,'spectrum');density=page.evaluate(PIXELS,'density')
      assert spectrum['opaque_pixels']>1000 and spectrum['chromatic_pixels']>100
      assert spectrum['luminance_range']>30 and spectrum['white_fraction']<.01
      assert density['luminance_range']>50 and density['white_fraction']<.01
      assert abs(spectrum['width']-round(spectrum['css_width']*dpr))<=1
      layout['pixel_checks'].append(dict(case=ci,comparison=comp,spectrum=spectrum,density=density));report['states_checked']+=1
    # Color encoding and exposure operate on a display copy; raw bytes and values remain stable.
    page.select_option('#case','0');page.select_option('#comparison','coupling')
    before=page.evaluate('JSON.stringify(a3AccuracyLab.data)')
    chip0=page.locator('.chip').evaluate_all('(es)=>es.map(e=>getComputedStyle(e).backgroundColor)')
    vals=page.evaluate('a3AccuracyLab.data.cases[0].treatments.coupled_mc.color.linear_srgb.map(a3AccuracyLab.encode)')
    assert all(0<=v<=255 for v in vals) and max(vals)<254 and min(vals)>20
    page.locator('#exposure').evaluate('(e)=>{e.value="2";e.dispatchEvent(new Event("input",{bubbles:true}))}')
    chip1=page.locator('.chip').evaluate_all('(es)=>es.map(e=>getComputedStyle(e).backgroundColor)')
    assert chip0!=chip1 and before==page.evaluate('JSON.stringify(a3AccuracyLab.data)')
    page.locator('#exposure').evaluate('(e)=>{e.value="0";e.dispatchEvent(new Event("input",{bubbles:true}))}')
    layout['zero_EV_coupled_srgb_codes']=vals;layout['exposure_changes_data']=False
    # Export actual downloaded files; compare exact bytes, not just parsed objects.
    if label=='desktop':
     for button,source in [('exportEvidence',ROOT/'data/a3/generated/benchmark.json'),('exportValidation',ROOT/'A3_VALIDATION.json')]:
      with page.expect_download() as event:page.click('#'+button)
      dest=OUT/(button+'-download.json');event.value.save_as(dest)
      assert dest.read_bytes()==source.read_bytes()
      report[button+'_sha256']=sha(dest);dest.unlink()
    page.evaluate('scrollTo(0,0)')
    layout['overflow_px']=page.evaluate('document.documentElement.scrollWidth-innerWidth');assert layout['overflow_px']==0
    page.screenshot(path=str(OUT/('lab-'+label+'.png')),full_page=True)
    report['layouts'].append(layout)
    # The source-reconciled A2 display is checked independently, without rewriting historic A2 validation.
    if label in ['desktop','mobile']:
     page.set_content((ROOT/'Open_Moon_A2_Accuracy_Lab.html').read_text(),wait_until='load');page.wait_for_function('window.a2AccuracyLab?.ready===true')
     for si in range(8):
      page.select_option('#sun',str(si))
      for mode in ['single','multiple']:
       page.select_option('#skyMode',mode);report['states_checked']+=1
     page.select_option('#sun','6');page.select_option('#skyMode','multiple')
     page.locator('#exposure').evaluate('(e)=>{e.value="0";e.dispatchEvent(new Event("input",{bubbles:true}))}')
     ids=page.locator('canvas').evaluate_all('(es)=>es.map(e=>e.id)')
     sky_id=next(i for i in ids if 'sky' in i.lower());px=page.evaluate(PIXELS,sky_id)
     assert px['luminance_range']>20 and px['white_fraction']<.01
     assert abs(px['width']-round(px['css_width']*dpr))<=2
     report['a2_pixel_checks'].append(dict(layout=label,canvas_id=sky_id,**px))
     page.locator('#'+sky_id).screenshot(path=str(OUT/('a2-corrected-sky-'+label+'.png')))
    context.close()
   browser.close()
  report['passed']=not any(report[k] for k in ['errors','console_errors','external_requests'])
 except Exception as e:report['errors'].append(type(e).__name__+': '+str(e))
 report['elapsed_s']=time.monotonic()-start
 (OUT/'lab-browser.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='layouts'},indent=2))
 return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
