"""Browser smoke, control, export and responsive-layout checks for the generated A2 lab."""
from pathlib import Path
import argparse, hashlib, json, time
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = ROOT.parents[1]/'illumination'/'references'

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=ROOT/'validation/a2')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    html = ROOT/'Open_Moon_A2_Accuracy_Lab.html'
    raw = html.read_bytes()
    report = {
        'schema': 'open-moon-a2-lab-browser/1',
        'html_sha256': hashlib.sha256(raw).hexdigest(),
        'errors': [], 'console_errors': [], 'external_requests': [],
        'states_checked': 0, 'passed': False,
    }
    start = time.monotonic()
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(
                executable_path='/usr/bin/chromium', headless=True,
                args=['--no-sandbox','--disable-dev-shm-usage'])
            page = browser.new_page(viewport={'width': 1280, 'height': 980}, accept_downloads=True)
            page.on('pageerror', lambda e: report['errors'].append(str(e)))
            page.on('console', lambda m: report['console_errors'].append(m.text) if m.type == 'error' else None)
            page.on('request', lambda r: report['external_requests'].append(r.url)
                    if not r.url.startswith(('blob:','data:')) else None)
            page.set_content(raw.decode(), wait_until='domcontentloaded')
            page.wait_for_function('window.a2AccuracyLab?.ready === true')
            assert page.locator('#ready').inner_text() == 'Offline · A2 checkpoint passed'
            assert page.locator('#statusBadge').inner_text() == 'reference checkpoint passed'

            # Exercise every solar elevation and both sky fields. Changing the render
            # mode must leave the stored evidence object untouched while redrawing.
            suns = page.locator('#sun option').count()
            for si in range(suns):
                page.select_option('#sun', str(si))
                for mode in ('single','multiple'):
                    page.select_option('#skyMode', mode)
                    page.wait_for_timeout(25)
                    assert page.locator('#ratioStat').inner_text().endswith('×')
                    report['states_checked'] += 1

            # Every fine-grid ray must render its spectrum and convergence annotation.
            rays = page.locator('#ray option').count()
            for ri in range(rays):
                page.select_option('#ray', str(ri))
                note = page.locator('#rayNote').inner_text()
                assert 'multiple/single photopic luminance ratio' in note
                assert 'standard/fine spectral L1' in note
                report['states_checked'] += 1

            # Exposure changes display encoding only.
            before = page.evaluate('JSON.stringify(a2AccuracyLab.data.sky.multiple_XYZ[6][12][16])')
            page.locator('#exposure').evaluate('(e)=>{e.value="2.5";e.dispatchEvent(new Event("input",{bubbles:true}))}')
            after = page.evaluate('JSON.stringify(a2AccuracyLab.data.sky.multiple_XYZ[6][12][16])')
            assert before == after
            assert page.locator('#evLabel').inner_text() == '+2.5 EV'
            report['states_checked'] += 1

            # Export validation and compare bytes semantically and by the embedded hash.
            with page.expect_download() as event:
                page.click('#exportValidation')
            dest = args.output/'validation-export.json'
            event.value.save_as(dest)
            exported = json.loads(dest.read_text())
            expected = json.loads((REFERENCES/'A2_VALIDATION.json').read_text())
            assert exported == expected
            report['validation_export_semantic_match'] = True
            report['validation_embedded_sha256'] = page.evaluate('a2AccuracyLab.data.validationSha')
            report['validation_source_sha256'] = hashlib.sha256((REFERENCES/'A2_VALIDATION.json').read_bytes()).hexdigest()
            assert report['validation_embedded_sha256'] == report['validation_source_sha256']
            dest.unlink()

            # Restore a representative review state before screenshots.
            page.select_option('#sun', '6')  # +45 deg
            page.select_option('#skyMode', 'multiple')
            page.select_option('#ray', '1')
            page.locator('#exposure').evaluate('(e)=>{e.value="-1.5";e.dispatchEvent(new Event("input",{bubbles:true}))}')
            report['desktop_overflow_px'] = page.evaluate('document.documentElement.scrollWidth-innerWidth')
            assert report['desktop_overflow_px'] == 0
            page.screenshot(path=str(args.output/'lab-desktop.png'), full_page=True)
            page.set_viewport_size({'width':390,'height':844})
            page.wait_for_timeout(250)
            report['mobile_overflow_px'] = page.evaluate('document.documentElement.scrollWidth-innerWidth')
            assert report['mobile_overflow_px'] == 0
            page.screenshot(path=str(args.output/'lab-mobile.png'), full_page=True)
            report['passed'] = not report['errors'] and not report['console_errors'] and not report['external_requests']
            browser.close()
    except Exception as exc:
        report['errors'].append(type(exc).__name__ + ': ' + str(exc))
    report['elapsed_s'] = time.monotonic() - start
    (args.output/'lab-browser.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 0 if report['passed'] else 1

if __name__ == '__main__':
    raise SystemExit(main())
