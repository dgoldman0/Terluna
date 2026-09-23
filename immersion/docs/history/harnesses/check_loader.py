"""Exercise the real HTML's dependency-failure UI in Chromium.

The sole external dependency request is deliberately aborted. No vendor code is
substituted. This tests the error/recovery UI, not Three.js or scene rendering.
"""
from pathlib import Path
import json,re
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
html=(ROOT/'Open_Moon_Shoreline.html').read_text()
record={'scope':__doc__,'integrated_scene_tested':False,'views':[]}
with sync_playwright() as p:
    browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
    for width,height in [(1440,900),(390,844)]:
        page=browser.new_page(viewport={'width':width,'height':height})
        errors=[];requests=[]
        page.on('pageerror',lambda e:errors.append(str(e)))
        def abort(route):
            requests.append(route.request.url)
            route.abort()
        page.route('https://cdn.jsdelivr.net/**',abort)
        page.set_content(html,wait_until='domcontentloaded',timeout=30000)
        page.locator('#boot-actions').wait_for(state='visible',timeout=25000)
        entry={'width':width,'height':height,'failure_status':page.locator('#boot-status').inner_text(),'recovery_visible':page.locator('#retry').is_visible(),'page_errors':errors}
        page.locator('#boot-about').click();entry['help_opens']=page.locator('#help').evaluate('(e)=>e.open');page.locator('#help-close').click();entry['help_closes']=not page.locator('#help').evaluate('(e)=>e.open')
        entry['horizontal_overflow']=page.evaluate('document.documentElement.scrollWidth > innerWidth')
        ids=re.findall(r'\bid="([^"]+)"',html)
        needed=set(re.findall(r"\$\('([^']+)'\)",(ROOT/'src/app.js').read_text()))
        entry['missing_control_ids']=sorted(needed-set(ids));entry['duplicate_ids']=sorted(k for k in set(ids) if ids.count(k)>1)
        entry['requests']=requests
        assert entry['recovery_visible'] and entry['help_opens'] and entry['help_closes'] and not entry['horizontal_overflow'] and not errors and not entry['missing_control_ids'] and not entry['duplicate_ids'],entry
        record['views'].append(entry);page.close()
    browser.close()
record['passed']=True
(ROOT/'validation/loader-ui-checks.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
