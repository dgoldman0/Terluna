"""The long-haul decision sheet: each ship class at the crown in perspective (from crown.py's renders) beside its berth
plan, drawn to scale from the sky-fleet study's product, with the numbers that separate the classes, and the side
elevation of all classes beside the Hindenburg.

    python3 visualization/sky-fleet/sheet.py      # build/plans/*.png, build/long-haul-classes.png and .html

Run crown.py first for the renders. Every dimension and number is read from
research/studies/sky_fleet/results/sky_fleet.json; the plans draw the hull outline the renders use (s^0.45 (1-s)^0.8,
which holds the flight models' volume) and the fins at the model's area.
"""
from __future__ import annotations
import base64
import hashlib
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT = ROOT / 'research' / 'studies' / 'sky_fleet' / 'results' / 'sky_fleet.json'
BUILD = HERE / 'build'
RENDERS = BUILD / 'renders'
FONT = 'Inter'
INK, MUTED, RULE = '#1f2a33', '#5d6b76', '#b8c2ca'
HULL, FIN, TOWER, ARM, WIND = '#dfe4e8', '#c3cad0', '#8a949c', '#a9a39a', '#2f6f9f'
PEAK = 0.45 / 1.25
NORM = PEAK ** 0.45 * (1 - PEAK) ** 0.8


def profile(s):
    s = np.clip(s, 0.0, 1.0)
    return s ** 0.45 * (1.0 - s) ** 0.8 / NORM


def fonts():
    names = {f.name for f in font_manager.fontManager.ttflist}
    family = FONT if FONT in names else 'DejaVu Sans'
    plt.rcParams.update({'font.family': family, 'font.size': 10, 'axes.edgecolor': RULE, 'text.color': INK,
                         'axes.labelcolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED})
    return family


def hull_outline(length, diameter, x0, y0):
    """The hull's plan outline, nose at (x0, y0), pointing west (the ship lies east of its nose)."""
    s = np.linspace(0.0, 1.0, 200)
    r = diameter / 2 * profile(s)
    top = np.c_[x0 + s * length, y0 + r]
    bottom = np.c_[x0 + s[::-1] * length, y0 - r[::-1]]
    return np.vstack([top, bottom])


def fins_outline(length, diameter, x0, y0):
    """The two horizontal fins in plan, each of the model's area, as crown.py draws them."""
    vol = 0.475 * length * diameter ** 2
    area = 0.25 * vol ** (2 / 3) / 4
    span = np.sqrt(0.5 * area)
    mean = area / span
    root, tip = 4 * mean / 3, 2 * mean / 3
    x_te = 0.965 * length
    x_le = x_te - root
    r_mid = diameter / 2 * profile((x_le + x_te) / 2 / length)
    shapes = []
    for side in (1, -1):
        xs = np.linspace(x_le, x_te, 12)
        root_pts = [(x0 + x, y0 + side * 0.97 * diameter / 2 * profile(x / length)) for x in xs]
        shapes.append(root_pts + [(x0 + x_te, y0 + side * (r_mid + span)), (x0 + x_te - tip, y0 + side * (r_mid + span))])
    return shapes


def dimension(ax, p0, p1, text, offset=(0, 0), colour=MUTED, size=8.5, vertical=False):
    arrow = FancyArrowPatch(p0, p1, arrowstyle='<|-|>', mutation_scale=7, lw=0.8, color=colour,
                            shrinkA=0, shrinkB=0)
    ax.add_patch(arrow)
    mx, my = (p0[0] + p1[0]) / 2 + offset[0], (p0[1] + p1[1]) / 2 + offset[1]
    ax.text(mx, my, text, ha='center', va='center', fontsize=size, color=colour, rotation=90 if vertical else 0,
            bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'))


def plan(option: dict, crown: dict, extent: tuple, path: Path):
    """One class's berths in plan at the crown's top level, to scale: the crown, the arms, each berthed ship with its
    fins, dimensions, north and the prevailing wind."""
    lay = option['crown']
    width = crown['width_at_band_bottom_m']
    r0 = width / 2 - 14.0
    fig, ax = plt.subplots(figsize=(7.2, 7.2 * (extent[3] - extent[2]) / (extent[1] - extent[0])), dpi=200)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect('equal')
    ax.axis('off')
    t = crown['terminal']
    r_t = np.sqrt(t['floor_m2'] / t['storeys'] / np.pi)
    ax.add_patch(Circle((0, 0), r_t, fc='#eef0f1', ec=TOWER, lw=1.2, zorder=3))
    ax.add_patch(Circle((0, 0), width / 2, fc='none', ec=TOWER, lw=0.6, ls=(0, (3, 2)), zorder=3))
    ax.text(0, 0, f"terminal\n{t['storeys']} storeys\n(frame dashed)", ha='center', va='center', fontsize=7.5,
            color=MUTED, zorder=4)
    per_level = []
    placed = 0
    for level in range(lay['levels']):
        row = []
        for side in (-1, 1):
            n = min(lay['ships_per_arm'], lay['berths'] - placed)
            if n <= 0:
                continue
            row.append((side, n))
            placed += n
        per_level.append(row)
    across = option['diameter_m'] + 2 * option['fin_span_m']
    full = max(per_level, key=lambda row: sum(n for _, n in row)) if per_level else []
    for side, n in full:
        arm_len = lay['pitch_m'] * (n - 0.5) + crown['arm_beyond'] * across
        y_end = side * (r0 + arm_len)
        ax.add_patch(Rectangle((-13, min(side * (r0 - 20), y_end)), 26, abs(y_end - side * (r0 - 20)),
                               fc=ARM, ec='#8d877e', lw=0.6, zorder=2))
        for k in range(n):
            y = side * (r0 + lay['pitch_m'] * (k + 0.5))
            ax.add_patch(Polygon(hull_outline(option['length_m'], option['diameter_m'], 36, y), closed=True,
                                 fc=HULL, ec='#7b868f', lw=0.7, zorder=5))
            for f in fins_outline(option['length_m'], option['diameter_m'], 36, y):
                ax.add_patch(Polygon(f, closed=True, fc=FIN, ec='#7b868f', lw=0.5, zorder=4))
            ax.add_patch(Rectangle((13, y - 5), 23, 10, fc='#6f777d', ec='none', zorder=5))
    # dimensions: a ship's length, the pitch, the arm's reach
    top_side, n_top = max(full, key=lambda sn: (sn[1], sn[0])) if full else (1, 1)
    y1 = top_side * (r0 + lay['pitch_m'] * (n_top - 0.5))
    ytxt = y1 + top_side * (option['diameter_m'] / 2 + 70)
    dimension(ax, (36, ytxt), (36 + option['length_m'], ytxt), f"{option['length_m']:.0f} m")
    if n_top > 1:
        ya, yb = top_side * (r0 + lay['pitch_m'] * 0.5), top_side * (r0 + lay['pitch_m'] * 1.5)
        xd = 36 + option['length_m'] + 70
        dimension(ax, (xd, ya), (xd, yb), f"{lay['pitch_m']:.0f} m", offset=(55, 0), vertical=True)
    arm_len = lay['pitch_m'] * (n_top - 0.5) + crown['arm_beyond'] * across
    xa = -150
    dimension(ax, (xa, top_side * r0), (xa, top_side * (r0 + arm_len)), f'arm {arm_len:.0f} m', offset=(-60, 0),
              vertical=True)
    # north and wind
    nx, ny = extent[0] + 0.08 * (extent[1] - extent[0]), extent[3] - 0.10 * (extent[3] - extent[2])
    ax.annotate('', xy=(nx, ny + 160), xytext=(nx, ny - 60), arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.2))
    ax.text(nx, ny + 205, 'N', ha='center', va='bottom', fontsize=10, color=INK, fontweight='bold')
    wx, wy = extent[0] + 0.08 * (extent[1] - extent[0]), extent[2] + 0.12 * (extent[3] - extent[2])
    ax.annotate('', xy=(wx + 330, wy), xytext=(wx, wy), arrowprops=dict(arrowstyle='-|>', color=WIND, lw=2.0))
    ax.text(wx, wy + 60, 'prevailing wind', ha='left', va='bottom', fontsize=8.5, color=WIND)
    # scale bar
    sx, sy = extent[1] - 0.08 * (extent[1] - extent[0]) - 500, extent[2] + 0.07 * (extent[3] - extent[2])
    for i in range(5):
        ax.add_patch(Rectangle((sx + 100 * i, sy), 100, 22, fc=INK if i % 2 == 0 else 'white', ec=INK, lw=0.6))
    for d in (0, 250, 500):
        ax.text(sx + d, sy - 30, f'{d}', ha='center', va='top', fontsize=7.5, color=MUTED)
    ax.text(sx + 250, sy + 50, 'metres', ha='center', va='bottom', fontsize=7.5, color=MUTED)
    levels = lay['levels']
    note = (f"{lay['berths']} berth{'s' if lay['berths'] > 1 else ''} on {levels} level{'s' if levels > 1 else ''}"
            + (f", {lay['level_spacing_m']:.0f} m apart; the fullest level shown" if levels > 1 else ''))
    ax.text(extent[0] + 0.02 * (extent[1] - extent[0]), extent[2] + 0.015 * (extent[3] - extent[2]), note,
            ha='left', va='bottom', fontsize=8, color=MUTED)
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(path, facecolor='white')
    plt.close(fig)


def plan_extent(options, crown):
    r0 = crown['width_at_band_bottom_m'] / 2
    ymax = max(r0 + o['crown']['pitch_m'] * min(o['crown']['ships_per_arm'], o['crown']['berths']) + 150
               for o in options)
    xmax = max(36 + o['length_m'] for o in options) + 250
    return (-r0 - 450, xmax, -ymax - 120, ymax + 120)


def facts(o: dict, lh: dict) -> list:
    n = o['network']
    lay = o['crown']
    return [
        ('Passengers per ship', f"{o['passengers']:,.0f}"),
        ('Departures an hour', f"{n['calls_per_hour'][0]:.1f}–{n['calls_per_hour'][1]:.1f}"),
        ('Destinations with a daily ship', f"{n['destinations_daily'][0]:.0f}–{n['destinations_daily'][1]:.0f}"),
        ('Berths at the crown', f"{lay['berths']} on {lay['levels']} level{'s' if lay['levels'] > 1 else ''}"),
        ('Longest arm', f"{lay['arm_length_m']:.0f} m"),
        ('Hydrogen in one ship', f"{o['hydrogen_t']:,.0f} t"),
        ('Energy per passenger-km', f"{o['wh_per_passenger_km']:.0f} Wh"),
        ('Ships cycling through the port', f"{n['fleet'][0]:,.0f}–{n['fleet'][1]:,.0f}"),
    ]


def sheet(options, lh, family):
    """A PNG sheet: the four perspectives, two by two, each with its class and numbers; the elevation below."""
    W = 3000
    pad, gap = 60, 40
    col_w = (W - 2 * pad - gap) // 2
    img_h = int(col_w * 1500 / 2400)
    text_h = 230
    lineup_h = int((W - 2 * pad) * 1270 / 2400)
    H = pad + 110 + 2 * (img_h + text_h + gap) + lineup_h + 150
    canvas = Image.new('RGB', (W, H), 'white')
    draw = ImageDraw.Draw(canvas)
    path = font_manager.findfont(font_manager.FontProperties(family=family))
    bold = font_manager.findfont(font_manager.FontProperties(family=family, weight='bold'))
    f_title, f_head, f_body, f_small = (ImageFont.truetype(bold, 46), ImageFont.truetype(bold, 32),
                                        ImageFont.truetype(path, 25), ImageFont.truetype(path, 22))
    draw.text((pad, pad), 'Long-haul ships at the crown, to scale: four classes', font=f_title, fill=INK)
    draw.text((pad, pad + 60), (f"Rigid hydrogen ships berthed above the eight-storey terminal in the flight band, for the "
                                f"crown's {lh['passengers_per_hour'][0]:,.0f}–{lh['passengers_per_hour'][1]:,.0f} long-haul "
                                f"passengers an hour; the tower as the summit tower study sized its form. Same camera for "
                                f"every class. A massing model: nothing here is structural design."), font=f_small, fill=MUTED)
    y0 = pad + 110
    for i, o in enumerate(options):
        cx = pad + (i % 2) * (col_w + gap)
        cy = y0 + (i // 2) * (img_h + text_h + gap)
        im = Image.open(RENDERS / f"view_{o['length_m']:.0f}.png").convert('RGB').resize((col_w, img_h), Image.LANCZOS)
        canvas.paste(im, (cx, cy))
        draw.text((cx, cy + img_h + 14), f"{o['length_m']:.0f} m ships, {o['diameter_m']:.0f} m across",
                  font=f_head, fill=INK)
        items = facts(o, lh)
        for j, (k, v) in enumerate(items):
            tx = cx + (j % 2) * (col_w // 2)
            ty = cy + img_h + 62 + (j // 2) * 40
            draw.text((tx, ty), k, font=f_small, fill=MUTED)
            draw.text((tx + col_w // 2 - 30, ty), v, font=f_body, fill=INK, anchor='ra')
    ly = y0 + 2 * (img_h + text_h + gap) + 10
    draw.text((pad, ly), 'Side by side, noses aligned: the Hindenburg (245 m) and the four classes',
              font=f_head, fill=INK)
    line = Image.open(RENDERS / 'lineup.png').convert('RGB')
    meta = json.loads((RENDERS / 'lineup_labels.json').read_text())
    lw = W - 2 * pad
    k = lw / line.width
    line = line.resize((lw, int(line.height * k)), Image.LANCZOS)
    ld = ImageDraw.Draw(line)
    names = {245.0: 'LZ 129 Hindenburg, 245 m (1936)'}
    for row in meta['rows']:
        label = names.get(row['length_m'], f"{row['length_m']:.0f} m")
        ld.text((row['x'] * k - 24, row['y'] * k), label, font=f_body, fill=INK, anchor='rm')
    canvas.paste(line, (pad, ly + 50))
    out = BUILD / 'long-haul-classes.png'
    canvas.save(out, optimize=True)
    return out


def page(options, lh, product_hash):
    """An HTML page: each class's perspective and plan side by side, with its numbers, then the elevation."""
    def data(p):
        return 'data:image/png;base64,' + base64.b64encode(Path(p).read_bytes()).decode()
    rows = []
    for o in options:
        items = ''.join(f'<tr><th>{html.escape(k)}</th><td>{html.escape(v)}</td></tr>' for k, v in facts(o, lh))
        rows.append(f"""
<section>
  <h2>{o['length_m']:.0f} m ships, {o['diameter_m']:.0f} m across</h2>
  <div class="pair">
    <figure><img src="{data(RENDERS / f"view_{o['length_m']:.0f}.png")}" alt="{o['length_m']:.0f} m ships at the crown, from the south-east"></figure>
    <figure><img src="{data(BUILD / 'plans' / f"plan_{o['length_m']:.0f}.png")}" alt="berth plan for {o['length_m']:.0f} m ships"></figure>
  </div>
  <table>{items}</table>
</section>""")
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Long-haul ships at the crown</title>
<style>
:root {{ --ink: {INK}; --muted: {MUTED}; --rule: {RULE}; --bg: #ffffff; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --ink: #e6ebef; --muted: #a8b3bc; --rule: #3b4650; --bg: #141a1f; }} }}
:root[data-theme="dark"] {{ --ink: #e6ebef; --muted: #a8b3bc; --rule: #3b4650; --bg: #141a1f; }}
body {{ margin: 0; background: var(--bg); color: var(--ink); font: 15px/1.5 Inter, system-ui, sans-serif; }}
main {{ max-width: 1400px; margin: 0 auto; padding: 24px 16px 64px; }}
h1 {{ font-size: 26px; margin: 8px 0 4px; }} h2 {{ font-size: 19px; margin: 28px 0 10px; }}
p {{ color: var(--muted); max-width: 900px; }}
.pair {{ display: grid; grid-template-columns: 1.6fr 1fr; gap: 12px; }}
@media (max-width: 800px) {{ .pair {{ grid-template-columns: 1fr; }} }}
figure {{ margin: 0; }} img {{ width: 100%; height: auto; display: block; border: 1px solid var(--rule); background: #fff; }}
table {{ border-collapse: collapse; margin-top: 10px; width: 100%; max-width: 900px; }}
th, td {{ text-align: left; padding: 4px 10px 4px 0; border-bottom: 1px solid var(--rule); font-weight: normal; }}
th {{ color: var(--muted); width: 55%; }} td {{ font-variant-numeric: tabular-nums; }}
footer {{ margin-top: 36px; color: var(--muted); font-size: 13px; }}
</style></head>
<body><main>
<h1>Long-haul ships at the crown, to scale</h1>
<p>Rigid hydrogen ships berthed in the flight band at the summit port's crown, one class at a time, for the crown's
long-haul traffic of {lh['passengers_per_hour'][0]:,.0f}–{lh['passengers_per_hour'][1]:,.0f} passengers an hour. Ships
lie nose-in on the lee side of arms running north and south, noses into the steady west wind, spaced so they can swing
{lh['crown']['swing_deg']:.0f}° together before their fins meet. The perspective uses the same camera for every class.
A massing model from the sky-fleet study's numbers: the tower's form follows the author's round diagrid, the arms are
deep trusses, and nothing here is structural design.</p>
{''.join(rows)}
<section><h2>Side by side, noses aligned: the Hindenburg (245 m) and the four classes</h2>
<figure><img src="{data(RENDERS / 'lineup.png')}" alt="side elevation of the Hindenburg and the four classes"></figure></section>
<footer>From research/studies/sky_fleet/results/sky_fleet.json (sha256 {product_hash}); renders by
visualization/sky-fleet/crown.py, plans and page by sheet.py.</footer>
</main></body></html>
"""
    out = BUILD / 'long-haul-classes.html'
    out.write_text(doc)
    return out


def main():
    family = fonts()
    product = json.loads(PRODUCT.read_text())
    lh = product['long_haul']
    options = lh['options']
    (BUILD / 'plans').mkdir(parents=True, exist_ok=True)
    ext = plan_extent(options, lh['crown'])
    for o in options:
        plan(o, lh['crown'], ext, BUILD / 'plans' / f"plan_{o['length_m']:.0f}.png")
    digest = hashlib.sha256(PRODUCT.read_bytes()).hexdigest()[:16]
    print(sheet(options, lh, family))
    print(page(options, lh, digest))


if __name__ == '__main__':
    main()
