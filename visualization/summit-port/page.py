"""The summit port page: the architectural model's views with numbered markers, and the stacking plan.

    python3 visualization/summit-port/page.py   # reads build/renders, writes build/summit-port-blueprints.html
"""
import json
from pathlib import Path

from PIL import Image

import drawings

HERE = Path(__file__).resolve().parent
ARCH = HERE / 'build' / 'renders'
IMG = HERE / 'build' / 'images'
IMG.mkdir(parents=True, exist_ok=True)
LABELS = json.loads((ARCH / 'labels.json').read_text())

VIEWS = [
    ('tower', 'The whole port', 'The frame model\'s full height, 24 km, from 60 km away. At this distance the building reads as a '
     'lattice with shelves: the occupied wings fill about one part in 2,000 of the space inside the frame. The views below '
     'step down in scale to show it working.',
     [('z3', 'Top of the open floors, 3 km'), ('z10', '10 km: enclosed floors below, regional docks above'),
      ('z15', '15 km: the sealed zone begins'), ('crown', 'The crown: long-haul terminal and docking arms in the flight band'),
      ('segment', 'The segment shown in view 3'), ('base', 'The base shown in view 5'),
      ('burj', 'The Burj Khalifa, 828 m, for scale')]),
    ('inside', 'Inside the frame, looking up', 'From the ground at the centre of the frame. The interior is open air, 7 km '
     'across at the base. Each of the six legs is a tower of its own, carrying a wing every 600 m; the wings step one leg '
     'round per band, so they climb the tower as a triple helix.',
     [('z3', 'Wings on the legs, stepping upward')]),
    ('segment', 'A segment at 6–7 km', 'One leg tower in the hotel and commerce zone, about 1.5 km away. Air taxis and '
     'regional ships move between the wings; the frame\'s face bracing crosses the view.',
     [('wing_6600', 'A wing: four storeys, 40 m deep, 1.25 km long in two arms'),
      ('lobby', 'Transfer lobby where the wing meets the leg'), ('leg', 'The leg: a lattice column 117 m across'),
      ('core', 'The lift-and-stair core inside the leg'), ('berths', 'Berth fingers every 60 m, with air taxis'),
      ('wing_7200', 'The next wing on this leg, 600 m up')]),
    ('wing', 'On a wing, at human scale', 'Standing off the outer rim of a wing at 6.6 km, looking back towards its leg. '
     'Four glazed storeys under a roof garden with a walk between two planted beds; berth fingers at the third storey; '
     'people on the roof walk are the dark marks.',
     [('roof', 'Roof garden and walk'), ('berth', 'Berth finger with an air taxi'), ('core', 'The leg\'s lattice and core')]),
    ('base', 'The base at a leg foot', 'One of the six leg feet from 2 km away. The Burj Khalifa stands beside it for scale: '
     'this leg\'s first wing is at 600 m and its second at 1.2 km.',
     [('gate', 'Gate interchange and footing at the leg foot'), ('wing_600', 'First wing on this leg, 600 m up'),
      ('wing_1200', 'Next wing, 1.2 km up'), ('rail', 'Rail round the base, joining the six gates')]),
    ('crown', 'The crown in the flight band', 'The top of the tower at 35 km above sea level. Long-haul ships 230 m long berth '
     'along eight docking arms on the side sheltered from the westerly; the sealed terminal wraps the top of the frame.',
     [('drum', 'The sealed long-haul terminal'), ('arm0', 'A docking arm with long-haul ships'),
      ('arm4', 'Arms step up 65 m apart, turned 137.5° from one another')]),
]


def jpeg(name):
    src = Image.open(ARCH / f'{name}.png').convert('RGB')
    src.save(IMG / f'{name}.jpg', quality=86, optimize=True, progressive=True)
    return f'images/{name}.jpg'


def figure(i, name, title, caption, marks):
    lab = LABELS[name]
    w, h = lab['size']
    src = jpeg(name)
    pins, items = [], []
    for n, (key, text) in enumerate(marks, 1):
        x, y, _ = lab['anchors'][key]
        x, y = min(max(x, 14), w - 14), min(max(y, 14), h - 14)
        pins.append(f'<span class="pin" style="left:{100 * x / w:.2f}%;top:{100 * y / h:.2f}%">{n}</span>')
        items.append(f'<li><span class="key">{n}</span>{text}</li>')
    return f'''
    <section class="view" id="{name}">
      <header class="view-head"><span class="view-no">View {i}</span><h2>{title}</h2></header>
      <div class="plate" style="aspect-ratio:{w}/{h}">
        <img src="{src}" alt="{title}: {caption}" width="{w}" height="{h}" loading="{'eager' if i < 3 else 'lazy'}">
        {''.join(pins)}
      </div>
      <div class="view-text">
        <p class="caption">{caption}</p>
        <ol class="keys">{''.join(items)}</ol>
      </div>
    </section>'''


def page():
    views = ''.join(figure(i, *v) for i, v in enumerate(VIEWS, 1))
    stacking = drawings.stacking()
    return f'''<title>Summit Port Blueprints</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Barlow+Condensed:wght@500;600&family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
  :root {{
    --paper: #eef2f5; --ink: #13222f; --muted: #4d6173; --rule: #c9d4dd; --accent: #1f5f8b; --pin: #c8553d;
    --pin-ink: #ffffff; --sheet: #0f3a5c;
    --display: 'Barlow Condensed', 'Arial Narrow', sans-serif;
    --body: 'IBM Plex Sans', system-ui, -apple-system, 'Segoe UI', sans-serif;
    --mono: 'IBM Plex Mono', ui-monospace, SFMono-Regular, Menlo, monospace;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{ --paper: #0b151e; --ink: #e0e9f0; --muted: #95a9ba; --rule: #22384a; --accent: #7fc4ec; color-scheme: dark; }}
  }}
  :root[data-theme="dark"] {{ --paper: #0b151e; --ink: #e0e9f0; --muted: #95a9ba; --rule: #22384a; --accent: #7fc4ec; color-scheme: dark; }}
  body {{ background: var(--paper); color: var(--ink); font-family: var(--body); font-size: 15px; line-height: 1.55; }}
  .wrap {{ max-width: 1180px; margin: 0 auto; padding-inline: 20px; padding-block: 40px 72px; }}
  .masthead {{ display: grid; gap: 14px; border-bottom: 1px solid var(--rule); padding-bottom: 26px; }}
  .kicker {{ font-family: var(--mono); font-size: 12px; letter-spacing: 0.12em; text-transform: uppercase; color: var(--accent); }}
  h1 {{ font-family: var(--display); font-weight: 600; font-size: clamp(34px, 5vw, 54px); line-height: 1.02; letter-spacing: 0.02em;
        text-transform: uppercase; margin: 0; text-wrap: balance; }}
  .lede {{ max-width: 70ch; margin: 0; }}
  .proposal {{ margin: 0; padding-left: 1.1em; max-width: 76ch; display: grid; gap: 4px; }}
  .facts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 1px; background: var(--rule);
            border: 1px solid var(--rule); margin-top: 6px; }}
  .fact {{ background: var(--paper); padding: 12px 14px; display: grid; gap: 2px; }}
  .fact b {{ font-family: var(--mono); font-weight: 600; font-size: 16px; font-variant-numeric: tabular-nums; }}
  .fact span {{ font-size: 12.5px; color: var(--muted); }}
  .view {{ margin-top: 48px; display: grid; gap: 12px; }}
  .view-head {{ display: flex; align-items: baseline; gap: 14px; flex-wrap: wrap; }}
  .view-no {{ font-family: var(--mono); font-size: 13px; color: var(--accent); letter-spacing: 0.08em; }}
  h2 {{ font-family: var(--display); font-weight: 600; font-size: 28px; letter-spacing: 0.04em; text-transform: uppercase; margin: 0; }}
  .plate {{ position: relative; width: 100%; max-width: 100%; border: 1px solid var(--rule); background: #dfe5ea; }}
  .plate img {{ display: block; width: 100%; height: 100%; object-fit: cover; }}
  .pin {{ position: absolute; transform: translate(-50%, -50%); width: 24px; height: 24px; border-radius: 50%;
          background: var(--pin); color: var(--pin-ink); font: 600 12px/24px var(--mono); text-align: center;
          box-shadow: 0 0 0 2px rgba(255,255,255,0.85); }}
  .view-text {{ display: grid; grid-template-columns: minmax(0, 1.1fr) minmax(0, 1fr); gap: 20px 36px; }}
  @media (max-width: 760px) {{ .view-text {{ grid-template-columns: 1fr; }} }}
  .caption {{ margin: 0; color: var(--muted); max-width: 64ch; }}
  .keys {{ list-style: none; margin: 0; padding: 0; display: grid; gap: 6px; }}
  .keys li {{ display: grid; grid-template-columns: 24px 1fr; gap: 10px; align-items: start; }}
  .key {{ width: 22px; height: 22px; border-radius: 50%; background: var(--pin); color: var(--pin-ink);
          font: 600 11.5px/22px var(--mono); text-align: center; }}
  .reading {{ margin-top: 56px; border-top: 1px solid var(--rule); padding-top: 22px; display: grid; gap: 10px; max-width: 76ch; }}
  .reading h3, .stack h3 {{ font-family: var(--display); font-size: 22px; text-transform: uppercase; letter-spacing: 0.05em; margin: 0; }}
  .reading p {{ margin: 0; }}
  .stack {{ margin-top: 48px; display: grid; gap: 12px; }}
  .drawing {{ background: var(--sheet); border: 1px solid var(--rule); overflow-x: auto; }}
  .drawing svg {{ display: block; width: 100%; min-width: 640px; height: auto; }}
  .notes {{ margin-top: 40px; color: var(--muted); max-width: 76ch; font-size: 14px; }}
  @media (max-width: 520px) {{ .wrap {{ padding-inline: 16px; }} .pin {{ width: 20px; height: 20px; line-height: 20px; font-size: 11px; }} }}
</style>
<div class="wrap">
  <header class="masthead">
    <div class="kicker">Open Moon · central port · architectural model</div>
    <h1>Summit Port Blueprints</h1>
    <p class="lede">The port modelled as a building. Its height, width at every level, floor bands, floor area and zones
    come from the tower model. The architecture on the frame is a proposal, drawn in so the port can be seen working:</p>
    <ul class="proposal">
      <li><b>Six leg towers.</b> Each leg of the hexagonal frame is a lattice column, 170 m across at the base and 117 m at
      6.6 km, with a lift-and-stair core inside.</li>
      <li><b>Wings.</b> Every 200 m a four-storey wing, 40 m deep, wraps two opposite legs, with arms along both faces.
      The pair steps one leg round per band, so each leg carries a wing every 600 m, where its shafts break and transfer.</li>
      <li><b>Docks and gardens.</b> Berth fingers every 60 m on each wing's outer rim; roof gardens on top.</li>
      <li><b>Base and crown.</b> A gate interchange at each leg foot, joined by rail; a sealed long-haul terminal at the
      crown with eight docking arms in the flight band.</li>
    </ul>
    <div class="facts">
      <div class="fact"><b>24 km</b><span>tall; 8.2 km across at the base</span></div>
      <div class="fact"><b>1.25 km</b><span>of wing per storey at 6.6 km, in two arms</span></div>
      <div class="fact"><b>117 m</b><span>across a leg at 6.6 km</span></div>
      <div class="fact"><b>1 in 2,000</b><span>of the frame's inside is building</span></div>
      <div class="fact"><b>1.5 million</b><span>people on a regular basis</span></div>
    </div>
  </header>
{views}
  <section class="reading">
    <h3>What the model shows about the design</h3>
    <p>The building occupies about one part in 2,000 of the space the frame encloses: 0.14 km³ of storeys inside 275 km³.
    That is why the whole tower reads as a lattice with shelves, and why the port only becomes legible at the scale of a
    leg and its wings. The architecture works at that scale: every wing reaches a core, every rim has berths, and ships
    have open air around every floor.</p>
    <p>If the port should read as one building from the highland, the floor would need gathering into fewer, deeper
    masses, such as districts of twenty storeys every 2 km in place of four storeys every 200 m, or the frame would need to
    shrink towards the floors. Both are choices about the design, and both change the frame's sizing.</p>
  </section>
  <section class="stack">
    <h3>Where the uses go</h3>
    <div class="drawing">{stacking}</div>
  </section>
  <p class="notes">Rendered from a Blender model built from the summit tower study's results. The frame, heights, widths,
  bands, floor area and zones are the model's; the leg columns, cores, wings, berths, gates, rail, terminal, arms, ships,
  trees and people are a proposal drawn at plausible sizes. The hexagonal frame is drawn to the square frame's widths and
  has not been sized. Nothing here is structural design.</p>
</div>
'''


if __name__ == '__main__':
    (HERE / 'build' / 'summit-port-blueprints.html').write_text(page())
    print('ok')
