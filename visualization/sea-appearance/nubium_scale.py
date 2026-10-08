"""Make a standalone angular-size comparison from the Nubium scene product."""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root))
    from shared.constants import EARTH_RADIUS, MOON_RADIUS

    source = root / 'research/studies/sea_appearance/results/nubium-midnight.json'
    product = json.loads(source.read_text())
    data = dict(earth_deg=product['scene']['earth']['diameter_deg'],
                earth_radius=EARTH_RADIUS, moon_radius=MOON_RADIUS,
                source='research/studies/sea_appearance/results/nubium-midnight.json',
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    html = '''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Earth and Moon at the same viewing scale</title>
<style>
body{margin:2rem auto;padding:0 1rem;max-width:960px;font:17px/1.5 system-ui;background:#151819;color:#eee}
h1{font-size:1.6rem}p{max-width:70ch}label{display:block;margin:1rem 0}input{width:min(100%,450px)}
.views{display:grid;grid-template-columns:1fr 1fr;gap:1rem}figure{margin:0}canvas{display:block;width:100%;background:#303437;border:1px solid #555}figcaption{margin:.5rem 0}small{color:#bec5c8}button{margin:.4rem .5rem .4rem 0;padding:.5rem;font:inherit}
@media(max-width:560px){body{margin:1rem auto}.views{gap:.5rem}h1{font-size:1.3rem}}
</style>
<h1>Same viewing scale</h1>
<p>Both panels contain exactly the same angular width of sky. The circles show physical disk size only; brightness, glare and surface detail are deliberately absent.</p>
<label for="fov">Width of sky in each panel: <strong id="degrees"></strong></label>
<input id="fov" type="range" min="5" max="90" step="0.5" value="10">
<div><button id="tight">10° crop</button><button id="wide">65° sea view</button></div>
<div class="views">
<figure><canvas id="moon" width="800" height="800"></canvas><figcaption>Moon from Earth<br><small id="mooninfo"></small></figcaption></figure>
<figure><canvas id="earth" width="800" height="800"></canvas><figcaption>Earth from the Moon<br><small id="earthinfo"></small></figcaption></figure>
</div>
<p id="ratio"></p>
<p>In the Nubium midnight scene, Earth is 53° above the horizon. A 10° crop centred on Earth cannot also contain the sea. The earlier portrait includes 65° horizontally and 87.4° vertically; displaying all of that sky in a small image shrinks every angular feature.</p>
<p>The supplied Moon photograph has unknown lens, crop and glare. Its luminous blob cannot establish its field of view. This comparison instead uses equal projection and equal Earth–Moon distance.</p>
<small id="provenance"></small>
<script>
const data = __DATA__;
const rad = Math.PI / 180;
const earth = data.earth_deg;
const moon = 2 * Math.asin(Math.sin(earth * rad / 2) * data.moon_radius / data.earth_radius) / rad;
const fov = document.getElementById('fov');
function draw(){
 const field=Number(fov.value); document.getElementById('degrees').textContent=field.toFixed(1)+'°';
 for(const [name,angle] of [['moon',moon],['earth',earth]]){
  const canvas=document.getElementById(name),ctx=canvas.getContext('2d');
  const fraction=Math.tan(angle*rad/2)/Math.tan(field*rad/2);
  ctx.clearRect(0,0,800,800);ctx.fillStyle='#e4e4df';ctx.beginPath();ctx.arc(400,400,400*fraction,0,Math.PI*2);ctx.fill();
  document.getElementById(name+'info').textContent=angle.toFixed(3)+'° across · '+(100*fraction).toFixed(2)+'% of panel width';
 }
 document.getElementById('ratio').textContent='Earth is '+(Math.tan(earth*rad/2)/Math.tan(moon*rad/2)).toFixed(2)+' times the diameter and about '+Math.round((earth/moon)**2)+' times the apparent disk area at equal distance. Zooming both panels changes their screen sizes together.';
}
fov.addEventListener('input',draw);
for(const [id,value] of [['tight',10],['wide',65]])document.getElementById(id).onclick=()=>{fov.value=value;draw()};
document.getElementById('provenance').textContent='Calculated from '+data.source+'; SHA-256 '+data.source_sha256+'. Rectilinear, centred disks; no atmospheric refraction.';
draw();
</script></html>'''.replace('__DATA__', json.dumps(data))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    print(args.out)


if __name__ == '__main__':
    main()
