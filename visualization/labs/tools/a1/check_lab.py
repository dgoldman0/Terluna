"""Browser smoke, control, export and layout checks for the generated A1 lab."""
from pathlib import Path
import argparse,hashlib,json,time
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[2]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'validation/a1');a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    html=ROOT/'Open_Moon_A1_Accuracy_Lab.html';raw=html.read_bytes();report={'schema':'open-moon-a1-lab-browser/1','html_sha256':hashlib.sha256(raw).hexdigest(),'errors':[],'external_requests':[],'states_checked':0,'passed':False};start=time.monotonic()
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--disable-dev-shm-usage'])
            page=browser.new_page(viewport={'width':1280,'height':980},accept_downloads=True)
            page.on('pageerror',lambda e:report['errors'].append(str(e)))
            page.on('request',lambda r:report['external_requests'].append(r.url) if not r.url.startswith(('blob:','data:')) else None)
            page.set_content(raw.decode(),wait_until='domcontentloaded');page.wait_for_function('window.a1AccuracyLab?.ready')
            for world in ['moon','earth','moon_no_ozone']:
                page.select_option('#world',world)
                for regime in ['fair','convection','fog']:
                    page.select_option('#regime',regime)
                    for sun in ['-6','45']:
                        page.select_option('#sun',sun)
                        assert page.evaluate('a1AccuracyLab.currentId')==world+'-'+regime
                        assert page.locator('#profileHash').inner_text()==page.evaluate("a1AccuracyLab.data.manifest.profiles.find(x=>x.id===a1AccuracyLab.currentId).profile_sha256")
                        report['states_checked']+=1
            page.select_option('#world','moon');page.select_option('#regime','fair');page.select_option('#sun','45')
            with page.expect_download() as event:page.click('#exportProfile')
            download=event.value;dest=a.output/'profile-export.json';download.save_as(dest)
            report['export_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();report['export_expected_sha256']=page.locator('#profileHash').inner_text();assert report['export_sha256']==report['export_expected_sha256'];dest.unlink()
            report['desktop_overflow_px']=page.evaluate('document.documentElement.scrollWidth-innerWidth');assert report['desktop_overflow_px']==0
            page.screenshot(path=str(a.output/'lab-desktop.png'),full_page=True)
            page.set_viewport_size({'width':390,'height':844});page.wait_for_timeout(250)
            report['mobile_overflow_px']=page.evaluate('document.documentElement.scrollWidth-innerWidth');assert report['mobile_overflow_px']==0
            page.screenshot(path=str(a.output/'lab-mobile.png'),full_page=True)
            report['passed']=not report['errors'] and not report['external_requests'];browser.close()
    except Exception as e:report['errors'].append(type(e).__name__+': '+str(e))
    report['elapsed_s']=time.monotonic()-start;(a.output/'lab-browser.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
