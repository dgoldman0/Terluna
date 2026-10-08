"""Drawings of ring 0's design (design.py): a key plan of the ring, Quarter 1's hall-level and roof plans unrolled
along the ring, and a cross-section through the Grand Market. Writes SVG (and PNG previews) into drawings/.

    python3 immersion/world/summit/ring0/draw.py
"""
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Polygon, Rectangle, Wedge  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as D  # noqa: E402

OUT = HERE / 'drawings'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'Inter', 'svg.fonttype': 'none', 'font.size': 8, 'axes.linewidth': 0.6})
INK, MID, LIGHT = '#222222', '#6b6b6b', '#b8b8b8'
FILL = {'gate': '#eadba6', 'square': '#f6f3ec', 'court': '#ffffff', 'node': '#d6d6d6', 'market': '#f0dcbc',
        'crafts': '#dde6cc', 'hotel': '#d3dfec', 'food': '#f1d0c2', 'services': '#e2dcea', 'bazaar': '#f4e3c3',
        'promenade': '#ece8e0', 'gallery': '#f7f5f0'}
ROOF = {'plaza': '#ebe6dd', 'orchard': '#dce9c6', 'beds': '#e8dfc8', 'void': '#ffffff', 'lawn': '#d6eac0',
        'terrace': '#e4dfd5', 'meadow': '#e4eecb', 'water': '#bcd8e6', 'pergola': '#e9e1cf', 'play': '#f0e2bd',
        'pavilion': '#ffffff'}
DATE = '28 September 2026'


def title_block(fig, title, subtitle):
    fig.text(0.015, 0.985, title, fontsize=15, weight='semibold', va='top', color=INK)
    fig.text(0.015, 0.958, subtitle, fontsize=8.5, va='top', color=MID)
    fig.text(0.985, 0.012, f'Summit port · ring 0 · design proposal (artistic), structure informed by the tower '
             f'study\'s hand sizing · {DATE}', fontsize=7, ha='right', color=MID)


def scale_bar(ax, x, y, lengths, unit='m', h=None, label=None):
    h = h or lengths[-1] * 0.025
    for i in range(len(lengths) - 1):
        ax.add_patch(Rectangle((x + lengths[i], y), lengths[i + 1] - lengths[i], h,
                               facecolor=INK if i % 2 == 0 else 'white', edgecolor=INK, lw=0.6))
    for v in lengths:
        ax.text(x + v, y - h * 0.8, f'{v:g}', ha='center', va='top', fontsize=6.5, color=INK)
    ax.text(x + lengths[-1] + lengths[-1] * 0.03, y + h / 2, label or unit, va='center', fontsize=6.5, color=INK)


# ------------------------------------------------------------------------------------------------------ key plan
def key_plan():
    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_axes([0.02, 0.04, 0.62, 0.86])
    ax.set_aspect('equal')
    ax.axis('off')
    ro, ri = D.R_OUTER, D.R_OUTER - D.WIDTH
    shades = ['#f0dcbc', '#e6d8ea', '#cfe4ea', '#e2e2d4', '#f1d0c2', '#d6eac0']
    for q, sh in zip(D.QUARTERS, shades):
        a0 = D.LEG0_DEG + 60 * (q['n'] - 1)
        ax.add_patch(Wedge((0, 0), ro, a0, a0 + 60, width=D.WIDTH, facecolor=sh, edgecolor=INK, lw=0.5))
        am = math.radians(a0 + 30)
        rl = ri - 700
        ax.text(rl * math.cos(am), rl * math.sin(am), f"Quarter {q['n']}\n{q['name']}\n", ha='center', va='center',
                fontsize=9, weight='semibold', color=INK, linespacing=1.3)
        ax.text(rl * math.cos(am), rl * math.sin(am) - 230, 'roof: ' + q['roof'], ha='center', va='center',
                fontsize=6.8, color=MID, style='italic')
    frame_r = 3367.0
    ax.add_patch(Circle((0, 0), frame_r, fill=False, edgecolor=MID, lw=0.5, ls=(0, (4, 3))))
    for k in range(24):
        a = math.radians(D.LEG0_DEG + 15 * k)
        c = (frame_r * math.cos(a), frame_r * math.sin(a))
        s = 45.0
        pts = [(c[0] + s * math.cos(a), c[1] + s * math.sin(a)), (c[0] - s * math.sin(a), c[1] + s * math.cos(a)),
               (c[0] - s * math.cos(a), c[1] - s * math.sin(a)), (c[0] + s * math.sin(a), c[1] - s * math.cos(a))]
        ax.add_patch(Polygon(pts, closed=True, facecolor=MID, edgecolor='none'))
    for g in D.GATES:
        a = math.radians(g['deg'])
        span = math.degrees(D.GATE_HALF / D.R_MID)
        ax.add_patch(Wedge((0, 0), ro + 30, g['deg'] - span, g['deg'] + span, width=D.WIDTH + 60,
                           facecolor='#d9b95c', edgecolor=INK, lw=0.6))
        ax.plot([ro * math.cos(a), (ro + 900) * math.cos(a)], [ro * math.sin(a), (ro + 900) * math.sin(a)],
                color=INK, lw=2.2, solid_capstyle='butt')
        ax.text((ro + 1180) * math.cos(a), (ro + 1180) * math.sin(a), f"Gate {g['n']}\nleg {g['n']} · faces {g['faces']}",
                ha='center', va='center', fontsize=8, color=INK)
    ax.text(0, 250, 'The hollow', ha='center', fontsize=11, color=MID, style='italic')
    ax.text(0, -120, 'inner harbour for sky boats and shuttles,\n6.4 km across, sheltered from the wind',
            ha='center', va='top', fontsize=7.5, color=MID)
    # north arrow
    ax.annotate('', xy=(-4300, 4500), xytext=(-4300, 3900), arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.2))
    ax.text(-4300, 4600, 'N', ha='center', fontsize=10, weight='semibold')
    scale_bar(ax, -4600, -4700, [0, 500, 1000, 2000], label='m')
    ax.set_xlim(-4900, 4900)
    ax.set_ylim(-4900, 4900)
    # the tower's elevation, with ring 0
    ex = fig.add_axes([0.70, 0.40, 0.27, 0.50])
    import numpy as np
    z = np.linspace(0, 24000, 200)
    w = 320 + (8193 - 320) * (1 - z / 24000) ** 1.5
    ex.fill_betweenx(z, -w / 2, w / 2, color='#eeeeee', lw=0)
    ex.plot(w / 2, z, color=MID, lw=0.7)
    ex.plot(-w / 2, z, color=MID, lw=0.7)
    rings = [3063, 4595, 6127, 7658, 9190, 10721]
    disks = [12253, 13785, 15316, 16848, 18380, 19911, 21443]
    for zr in rings:
        ww = 320 + (8193 - 320) * (1 - zr / 24000) ** 1.5
        ex.plot([-ww / 2, -ww / 2 + 300], [zr, zr], color=INK if zr == 3063 else MID, lw=2.4 if zr == 3063 else 1.2)
        ex.plot([ww / 2 - 300, ww / 2], [zr, zr], color=INK if zr == 3063 else MID, lw=2.4 if zr == 3063 else 1.2)
    for zd in disks:
        ww = 320 + (8193 - 320) * (1 - zd / 24000) ** 1.5
        ex.plot([-ww / 2, ww / 2], [zd, zd], color=MID, lw=1.2)
    ex.annotate('ring 0, 3.06 km', xy=(3400, 3063), xytext=(4300, 6500), fontsize=8,
                arrowprops=dict(arrowstyle='-', color=INK, lw=0.6))
    ex.text(0, 24800, 'The tower: six rings, seven disks, the terminal', ha='center', fontsize=8, color=MID)
    ex.set_xlim(-5200, 7200)
    ex.set_ylim(-300, 25500)
    ex.axis('off')
    notes = ('Ring 0 stands on the transfer ring, 3.06 km above the summit:\n'
             '150 m wide and 20.5 km round at mid-width.\n\n'
             'Air 8.9 °C, steady day and night; oxygen like 2,450 m on\n'
             'Earth, so it is open to everyone.\n\n'
             'Two storeys (9 m and 6 m) under a park. It holds the base\n'
             'programme (6 km²), so the ground between the legs stays open.\n\n'
             'The six gates, where the legs arrive, divide it into six\n'
             'quarters of 3.4 km. The frame\'s 24 nodes per level stand at\n'
             'the outer edge (grey diamonds), 853 m apart at mid-width.')
    fig.text(0.66, 0.30, notes, fontsize=8, va='top', color=INK, linespacing=1.6)
    title_block(fig, 'Ring 0 · key plan', 'The six gates and quarters; Quarter 1 is drawn in detail on the next sheets')
    save(fig, 'ring0_key_plan')


# ----------------------------------------------------------------------------------------------- unrolled plans
ROW_H = D.WIDTH
GAP = 95.0


def row_frame(ax, y0, s0, s1, label):
    """One row of an unrolled plan: s0..s1 along the ring (drawn from x = 0), the outer edge at the top."""
    ax.text(0, y0 + ROW_H + 58, label, fontsize=8.5, weight='semibold', color=INK, va='bottom')
    names = {0: 'Gate 1: leg 1 and a frame node', 1: 'frame node A', 2: 'frame node B', 3: 'frame node C',
             4: 'Gate 2: leg 2 and a frame node'}
    for k in range(0, 5):
        s = k * D.NODE_SPACING
        if s0 - 1 <= s <= s1 + 1:
            x = s - s0
            c = (x, y0 + ROW_H + 30)
            ax.add_patch(Polygon([(c[0], c[1] + 22), (c[0] + 15, c[1]), (c[0], c[1] - 22), (c[0] - 15, c[1])],
                                 closed=True, facecolor=MID, edgecolor='none'))
            left = x > (s1 - s0) / 2
            ax.text(x + (-22 if left else 22), y0 + ROW_H + 30, names[k], fontsize=6, color=MID, va='center',
                    ha='right' if left else 'left')


def to_y(y0, d):
    return y0 + ROW_H - d


def place_rect(ax, y0, s0, s1, d0, d1, face, lw=0.6, hatch=None, ec=INK):
    ax.add_patch(Rectangle((s0, to_y(y0, d1)), s1 - s0, d1 - d0, facecolor=face, edgecolor=ec, lw=lw, hatch=hatch))


def label_in(ax, y0, s0, s1, d0, d1, text, note='', size=7):
    cx, cy = (s0 + s1) / 2, to_y(y0, (d0 + d1) / 2)
    if s1 - s0 < 75:
        ax.text(cx, cy, text, rotation=90, ha='center', va='center', fontsize=6.2, color=INK)
        return
    ax.text(cx, cy + (5 if note else 0), text, ha='center', va='center', fontsize=size, color=INK, weight='medium')
    if note and s1 - s0 > 180:
        ax.text(cx, cy - 9, note, ha='center', va='center', fontsize=5.6, color=MID, wrap=True)


def hall_level():
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.03, 0.05, 0.94, 0.86])
    ax.set_aspect('equal')
    ax.axis('off')
    rows = [(k * D.NODE_SPACING, (k + 1) * D.NODE_SPACING) for k in range(4)]
    labels = ['Gate 1 to node A', 'Node A to node B', 'Node B to node C', 'Node C to Gate 2']
    for r, ((s0, s1), lab) in enumerate(zip(rows, labels)):
        y0 = (3 - r) * (ROW_H + GAP)
        off = -s0
        row_frame(ax, y0, s0, s1, f'{r + 1}. {lab}  (s = {s0:.0f}–{s1:.0f} m)')
        # galleries and the promenade, continuous
        place_rect(ax, y0, 0, s1 - s0, 0, 8, FILL['gallery'], lw=0.4)
        place_rect(ax, y0, 0, s1 - s0, 142, 150, FILL['gallery'], lw=0.4)
        place_rect(ax, y0, 0, s1 - s0, 62, 88, FILL['promenade'], lw=0.4)
        for name, kind, a, b, d0, d1, note in D.Q1_L1:
            a2, b2 = max(a, s0), min(b, s1)
            if b2 <= a2:
                continue
            hatch = '////' if kind == 'court' else None
            place_rect(ax, y0, a2 + off, b2 + off, d0, d1, FILL[kind], hatch=hatch,
                       ec=LIGHT if kind == 'court' else INK)
            label_in(ax, y0, a2 + off, b2 + off, d0, d1, name, note)
        # promenade label and flow arrows
        open_ = [(a, b) for (_, k, a, b, d0, d1, _) in D.Q1_L1 if k in ('gate', 'square', 'court') or (d0 <= 62 and d1 >= 88)]
        for sx in range(int(s0) + 140, int(s1) - 60, 300):
            if not any(a - 20 <= sx <= b + 20 for a, b in open_):
                ax.text(sx - s0, to_y(y0, 75), 'promenade', fontsize=5.8, color=MID, ha='center', va='center',
                        style='italic')
        # facades: outer and inner glass lines
        for d in (0.0, 150.0):
            ax.plot([0, s1 - s0], [to_y(y0, d)] * 2, color=INK, lw=1.6)
            ax.plot([0, s1 - s0], [to_y(y0, d + (1.2 if d == 0 else -1.2))] * 2, color=INK, lw=0.4)
        # column grid: dots at the column lines every 15 m of the outer edge (14.7 m at mid-width)
        import numpy as np
        ss = np.arange(0, s1 - s0, 14.66)
        for d in D.COLUMN_LINES:
            ax.scatter(ss, [to_y(y0, d)] * len(ss), s=0.6, color=INK, marker='s', linewidths=0)
        for st in D.Q1_STATIONS:
            if s0 - 1 <= st <= s1 + 1:
                c = (st + off, to_y(y0, 75))
                ax.add_patch(Circle(c, 7.5, facecolor='white', edgecolor=INK, lw=0.8))
                ax.text(c[0], c[1], 'U', ha='center', va='center', fontsize=6, weight='bold')
        ax.text(-12, to_y(y0, 4), 'outer gallery', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(-12, to_y(y0, 35), 'outer band', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(-12, to_y(y0, 75), 'promenade', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(-12, to_y(y0, 115), 'inner band', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(-12, to_y(y0, 146), 'inner gallery', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(s1 - s0 + 10, to_y(y0, -10), 'frame side', fontsize=5.5, color=MID)
        ax.text(s1 - s0 + 10, to_y(y0, 162), 'the hollow', fontsize=5.5, color=MID)
    scale_bar(ax, 0, -70, [0, 25, 50, 100, 200], label='m')
    ax.text(260, -70, 'U  Undercroft Line station (the people mover runs under the promenade)     '
            'hatched: light courts open to the park above     dots: columns on the trusses\' panel points',
            fontsize=6.5, color=MID, va='bottom')
    legend = [('gate', 'gate hall'), ('market', 'market'), ('crafts', 'crafts, flowers'), ('hotel', 'hotel'),
              ('food', 'food, restaurants'), ('bazaar', 'bazaar'), ('services', 'services'), ('node', 'node room')]
    for i, (k, lab) in enumerate(legend):
        x = 260 + i * 70
        ax.add_patch(Rectangle((x, -118), 14, 12, facecolor=FILL[k], edgecolor=INK, lw=0.5))
        ax.text(x + 18, -112, lab, fontsize=6.3, va='center')
    ax.set_xlim(-75, D.NODE_SPACING + 55)
    ax.set_ylim(-135, 4 * (ROW_H + GAP) + 20)
    title_block(fig, 'Ring 0 · Quarter 1 (Market Quarter) · hall level (L1)',
                'Unrolled along the mid-width circle (radius 3.26 km): widths true, lengths along the ring. '
                'Outer edge (frame) at the top of each row, the hollow at the bottom. L2 above each place is noted.')
    save(fig, 'ring0_q1_hall_level')


def roof_level():
    import numpy as np
    rng = np.random.default_rng(7)
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.03, 0.05, 0.94, 0.86])
    ax.set_aspect('equal')
    ax.axis('off')
    rows = [(k * D.NODE_SPACING, (k + 1) * D.NODE_SPACING) for k in range(4)]
    labels = ['Gate 1 to node A', 'Node A to node B', 'Node B to node C', 'Node C to Gate 2']
    for r, ((s0, s1), lab) in enumerate(zip(rows, labels)):
        y0 = (3 - r) * (ROW_H + GAP)
        off = -s0
        L = s1 - s0
        row_frame(ax, y0, s0, s1, f'{r + 1}. {lab}')
        place_rect(ax, y0, 0, L, 0, 10, ROOF['plaza'], lw=0.4)           # outer walk
        place_rect(ax, y0, 0, L, 130, 150, ROOF['plaza'], lw=0.4)        # inner walk
        place_rect(ax, y0, 0, L, 62, 88, ROOF['plaza'], lw=0.4)          # skylight walks
        place_rect(ax, y0, 0, L, 69, 81, '#d2e5ef', lw=0.5, hatch='----', ec=MID)   # the skylight
        order = sorted(D.Q1_ROOF, key=lambda e: e[1] in ('water', 'pavilion'))   # small features drawn last
        for name, kind, a, b, d0, d1 in order:
            if kind == 'plaza' and 'Gate' in name:
                a = max(a, D.GATE_HALF) if a < D.GATE_HALF else a
                b = min(b, D.QUARTER_LENGTH - D.GATE_HALF) if b > D.QUARTER_LENGTH - D.GATE_HALF else b
            a2, b2 = max(a, s0), min(b, s1)
            if b2 <= a2:
                continue
            x0, x1 = a2 + off, b2 + off
            hatch = {'void': 'xxxx', 'pavilion': '....', 'pergola': '||||'}.get(kind)
            place_rect(ax, y0, x0, x1, d0, d1, ROOF[kind], hatch=hatch, ec=LIGHT if kind == 'void' else INK)
            if kind == 'orchard':
                for s in np.arange(x0 + 6, x1 - 3, 7.5):
                    for d in np.arange(d0 + 6, d1 - 3, 7.5):
                        ax.add_patch(Circle((s, to_y(y0, d)), 2.9, facecolor='#a9c98a', edgecolor='#4f7a3a', lw=0.3))
            if kind == 'beds':
                for d in np.arange(d0 + 4, d1 - 2, 4.5):
                    ax.plot([x0 + 3, x1 - 3], [to_y(y0, d)] * 2, color='#9d8f6c', lw=0.9)
            if kind == 'meadow':
                n = int((x1 - x0) * (d1 - d0) / 90)
                ax.scatter(rng.uniform(x0 + 2, x1 - 2, n), [to_y(y0, v) for v in rng.uniform(d0 + 2, d1 - 2, n)],
                           s=0.8, color='#7ea05a', linewidths=0)
            if kind in ('lawn', 'play', 'plaza', 'terrace') and x1 - x0 > 60:
                for s in np.arange(x0 + 20, x1 - 10, 38):
                    ax.add_patch(Circle((s, to_y(y0, d0 + 8)), 4.5, facecolor='#a9c98a', edgecolor='#4f7a3a', lw=0.3))
            if kind in ('water', 'pavilion'):
                ax.text((x0 + x1) / 2, to_y(y0, d1 + 3.5), name, fontsize=5.8, ha='center', va='top', color=INK)
            elif kind == 'void':
                for dd in (30.0, 105.0):
                    place_rect(ax, y0, x0, x1, dd - 3, dd + 3, ROOF['plaza'], lw=0.5)
                ax.text((x0 + x1) / 2, to_y(y0, 70), name.replace(' void', ''), rotation=90, fontsize=6, ha='center',
                        va='center', color=MID)
            else:
                small = [(a2_, b2_, e0, e1) for (_, k2, a2_, b2_, e0, e1) in D.Q1_ROOF if k2 in ('water', 'pavilion')]
                inside = any(a <= sa and sb <= b and d0 <= e0 and e1 <= d1 for sa, sb, e0, e1 in small)
                dshift = 14 if inside else 0
                cx = (x0 + x1) / 2
                if x1 - x0 < 75:
                    ax.text(cx, to_y(y0, (d0 + d1) / 2), name, rotation=90, ha='center', va='center', fontsize=6.2)
                else:
                    ax.text(cx, to_y(y0, d0 + 12 + dshift), name, ha='center', va='center', fontsize=6.8,
                            weight='medium', color=INK)
        # avenue along the outer walk
        for s in np.arange(9, L - 5, 15):
            ax.add_patch(Circle((s, to_y(y0, 12.5)), 3.8, facecolor='#9fc27f', edgecolor='#4f7a3a', lw=0.3))
        # kiosks and pavilions
        for name, s, d, count in D.Q1_KIOSKS:
            if not (s0 - 1 <= s <= s1 + 1):
                continue
            if name == 'kiosks':
                cols = int(math.ceil(count / 2))
                for i in range(count):
                    x = s + off - cols * 4.5 + (i % cols) * 9
                    yy = to_y(y0, d + (i // cols) * 9)
                    ax.add_patch(Rectangle((x, yy), 5, 5, facecolor='#c0543f', edgecolor=INK, lw=0.4))
                ax.text(s + off, to_y(y0, d + (count // cols) * 9 + 4), 'open-air shops', fontsize=5.4, ha='center',
                        va='top', color=INK)
            else:
                ax.add_patch(Rectangle((s + off - 12, to_y(y0, d + 8)), 24, 16, facecolor='white', edgecolor=INK,
                                       lw=0.8, hatch='....'))
                ax.text(s + off, to_y(y0, d - 12), name, fontsize=5.6, ha='center', color=INK)
        for d in (0.6, 149.4):
            ax.plot([0, L], [to_y(y0, d)] * 2, color=INK, lw=0.8)
        for g in (0.0, D.QUARTER_LENGTH):                  # the gate halls' glass vaults rise through the plazas
            a, b = max(g - D.GATE_HALF, s0), min(g + D.GATE_HALF, s1)
            if b > a:
                place_rect(ax, y0, a - s0, b - s0, 18, 132, '#d2e5ef', lw=0.9, hatch='++', ec='#3d6f8e')
                ax.text((a + b) / 2 - s0, to_y(y0, 75), 'gate vault\n(glass)', ha='center', va='center', fontsize=6.2,
                        color='#1f4e6b')
        ax.text(-12, to_y(y0, 5), 'outer walk', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(-12, to_y(y0, 75), 'skylight', fontsize=5.5, ha='right', va='center', color=MID)
        ax.text(-12, to_y(y0, 140), 'inner walk', fontsize=5.5, ha='right', va='center', color=MID)
    scale_bar(ax, 0, -70, [0, 25, 50, 100, 200], label='m')
    ax.text(260, -70, 'Railings 2.2 m high at both edges (lunar gravity). Lanes\' voids are railed and bridged. '
            'Trees and plants are placeholders: the Open Moon\'s plants are not yet designed.',
            fontsize=6.5, color=MID, va='bottom')
    ax.set_xlim(-75, D.NODE_SPACING + 55)
    ax.set_ylim(-135, 4 * (ROW_H + GAP) + 20)
    title_block(fig, 'Ring 0 · Quarter 1 (Market Quarter) · roof park',
                'Orchards and kitchen gardens; open to everyone. Air 8.9 °C, steady. Unrolled as the hall-level sheet.')
    save(fig, 'ring0_q1_roof')


# --------------------------------------------------------------------------------------------------------- section
def person(ax, d, z, h=1.75, color=INK):
    ax.add_patch(Circle((d, z + h - 0.13), 0.13, facecolor=color, edgecolor='none'))
    ax.add_patch(Polygon([(d - 0.2, z), (d + 0.2, z), (d + 0.18, z + h - 0.3), (d - 0.18, z + h - 0.3)],
                         closed=True, facecolor=color, edgecolor='none'))


def tree(ax, d, z, h=8.0, w=5.0):
    ax.plot([d, d], [z, z + h * 0.45], color='#6b5238', lw=1.2)
    ax.add_patch(Circle((d, z + h * 0.7), w / 2, facecolor='#b7d39a', edgecolor='#4f7a3a', lw=0.5))


def section():
    fig = plt.figure(figsize=(17, 9))
    ax = fig.add_axes([0.03, 0.08, 0.94, 0.80])
    ax.set_aspect('equal')
    ax.axis('off')
    poche = '#3a3a3a'
    # the transfer ring's top: its inner chords at d = 27 (the outer ones are 120 m further out, off the sheet),
    # the bracing between them at its top, on which the trusses' roots bear; cut short below
    ax.add_patch(Circle((27.0, -27.0 - 2.35), 2.35, facecolor=poche, edgecolor=INK))
    ax.plot([-60, 27], [-27.0, -27.0], color=INK, lw=1.2)
    ax.plot([27, 27], [-27 - 4.7, -46], color=INK, lw=1.4)
    for k in range(5):
        x0 = 27 - 18.0 * k
        ax.plot([x0, x0 - 18.0], [-27.0, -46], color=MID, lw=0.7)
        ax.plot([x0, x0], [-27.0, -46], color=LIGHT, lw=0.5)
    ax.plot([-62, 32], [-46.5, -46.5], color=INK, lw=0.6, ls=(0, (6, 2, 1, 2)))
    ax.text(-20, -40, 'transfer ring: 280 m deep, 120 m across\n(its outer chords are off the sheet; cut short)',
            ha='center', fontsize=6.5, color=INK)
    ax.add_patch(Rectangle((0, -27.0), 27, 1.2, facecolor=MID, edgecolor='none'))
    ax.text(13, -24.5, 'truss roots bear here', ha='center', fontsize=5.8, color=MID)
    # the frame beyond the outer edge: two diagrid members leaving the node below and the node above (the cut is
    # between nodes, so only their far parts show)
    for x0, z0, x1, z1 in ((-48, -30, -20, 35), (-20, -30, -48, 35)):
        ax.plot([x0, x1], [z0, z1], color=LIGHT, lw=5, solid_capstyle='butt')
    ax.text(-34, 37, 'diagrid members beyond\n(48 m wide, 873 m between nodes)', ha='center', fontsize=6.2, color=MID)
    # radial truss: top chord under L1, bottom chord 25 m deep at the frame to 2 m at the edge
    pts = [0.0] + D.COLUMN_LINES + [150.0]
    depth = lambda d: 25.0 + (2.0 - 25.0) * d / 150.0
    ax.plot([0, 150], [-0.6, -0.6], color=INK, lw=1.6)
    ax.plot([0, 150], [-0.6 - depth(0), -0.6 - depth(150)], color=INK, lw=1.6)
    for d in pts:
        ax.plot([d, d], [-0.6, -0.6 - depth(d)], color=INK, lw=0.9)
    for a, b in zip(pts[:-1], pts[1:]):
        if (a, b) == (62.0, 88.0):
            continue                                   # a Vierendeel panel round the people mover
        ax.plot([a, b], [-0.6, -0.6 - depth(b)], color=INK, lw=0.8)
    # the Undercroft Line
    lo, hi = D.UNDERCROFT_LINE
    ax.add_patch(Rectangle((lo, -0.6 - depth(75) + 1.0), hi - lo, 4.2, facecolor='#e7e2d6', edgecolor=INK, lw=0.7))
    ax.text(75, -0.6 - depth(75) - 1.6, 'Undercroft Line (people mover)', ha='center', va='top', fontsize=6.5)
    ax.text(110, -12, 'radial trusses every 15 m,\n25 m deep at the frame to 2 m at the edge', fontsize=6.5, color=INK)
    # slabs
    for z0, z1 in ((-0.6, 0.0), (8.6, 9.0), (14.4, 15.0)):
        ax.add_patch(Rectangle((0, z0), 150, z1 - z0, facecolor=poche, edgecolor='none'))
    # void over the promenade in L2 and the roof: cut the slabs there
    ax.add_patch(Rectangle((62, 8.55), 26, 0.5, facecolor='white', edgecolor='none'))
    ax.add_patch(Rectangle((66, 14.35), 18, 0.7, facecolor='white', edgecolor='none'))
    # columns
    for d in D.COLUMN_LINES:
        if d in (62.0, 88.0):
            ax.add_patch(Rectangle((d - 0.4, 0), 0.8, 15, facecolor=MID, edgecolor='none'))
        else:
            ax.add_patch(Rectangle((d - 0.35, 0), 0.7, 15, facecolor=MID, edgecolor='none'))
    # facades: outer and inner glass, full height; mullion ticks
    for d in (0.5, 149.5):
        ax.plot([d, d], [0, 15], color='#3d6f8e', lw=1.6)
        for z in (3.0, 6.0, 12.0):
            ax.plot([d - 0.6, d + 0.6], [z, z], color='#3d6f8e', lw=0.6)
    # outer gallery mezzanine at L2, inner gallery balcony over the hollow
    ax.add_patch(Rectangle((0.5, 8.6), 3.5, 0.4, facecolor=poche))
    ax.add_patch(Rectangle((149.5, 8.6), 3.5, 0.4, facecolor=poche))
    ax.plot([153, 153], [9, 10.2], color=INK, lw=0.6)
    # L2 galleries' balustrades along the promenade
    for d in (62.0, 88.0):
        ax.plot([d, d], [9, 10.1], color=INK, lw=0.6)
    # market stalls (L1) and restaurant tables (L2) in both bands
    for band in ((10, 60), (90, 140)):
        for d in range(band[0], band[1], 12):
            ax.add_patch(Rectangle((d, 0), 3.0, 2.4, facecolor='#e8c89a', edgecolor=INK, lw=0.4))
            ax.plot([d - 0.4, d + 3.4], [2.9, 2.6], color='#c0543f', lw=1.4)
        for d in range(band[0] + 3, band[1], 9):
            ax.plot([d, d + 1.2], [9.75, 9.75], color=INK, lw=0.8)
            ax.plot([d + 0.6, d + 0.6], [9.0, 9.75], color=INK, lw=0.5)
    # the skylight: a shallow glazed vault over the promenade, pergolas along its walks
    import numpy as np
    a0, a1 = D.SKYLIGHT
    xs = np.linspace(a0 - 3, a1 + 3, 40)
    ax.plot(xs, 15.0 + 3.5 * np.sin(np.pi * (xs - (a0 - 3)) / (a1 - a0 + 6)), color='#3d6f8e', lw=1.2)
    for d in (64.0, 86.0):
        ax.plot([d, d], [15.6, 18.6], color='#6b5238', lw=0.8)
        ax.plot([d - 2, d + 2], [18.6, 18.6], color='#6b5238', lw=1.0)
    # roof park: soil, orchard trees, kitchen beds, walks, railings
    ax.add_patch(Rectangle((10, 15.0), 52, 1.2, facecolor='#8a6d4b', edgecolor='none'))
    ax.add_patch(Rectangle((88, 15.0), 42, 0.9, facecolor='#8a6d4b', edgecolor='none'))
    for d in np.arange(16, 60, 7.5):
        tree(ax, d, 16.2, h=8.5, w=5.2)
    tree(ax, 12.5, 16.2, h=11, w=6.5)
    for d in np.arange(91, 128, 4.5):
        ax.add_patch(Rectangle((d, 15.9), 3.0, 0.6, facecolor='#9d8f6c', edgecolor=INK, lw=0.3))
    for d in (0.8, 149.2):
        ax.plot([d, d], [15.0, 17.2], color=INK, lw=0.8)
    # people for scale
    for d, z in ((5, 0), (30, 0), (70, 0), (76, 0), (82, 0), (118, 0), (145, 0), (66, 9), (95, 9), (151, 9),
                 (4, 15), (74, 15.2), (100, 15.9), (140, 15), (2, 9)):
        person(ax, d, z)
    # dimensions
    for (name, a, b) in D.ZONES:
        ax.annotate('', xy=(a, 24.5), xytext=(b, 24.5), arrowprops=dict(arrowstyle='<->', lw=0.5, color=INK))
        ax.text((a + b) / 2, 25.3, f'{name}\n{b - a:g} m', ha='center', va='bottom', fontsize=6.3)
    for z0, z1, lab in ((0, 9, 'L1  9 m'), (9, 15, 'L2  6 m')):
        ax.annotate('', xy=(158, z0), xytext=(158, z1), arrowprops=dict(arrowstyle='<->', lw=0.5, color=INK))
        ax.text(159, (z0 + z1) / 2, lab, va='center', fontsize=6.5)
    ax.text(0, -52, '← frame', fontsize=7, color=MID)
    ax.text(150, -52, 'the hollow →', fontsize=7, color=MID, ha='right')
    ax.text(36, 27.5 + 6, 'Orchard', fontsize=7, ha='center', style='italic')
    ax.text(109, 21.5, 'Kitchen gardens', fontsize=7, ha='center', style='italic')
    ax.text(75, 20.5, 'skylight', fontsize=6.5, ha='center', color='#3d6f8e')
    ax.text(35, 4.5, 'market stalls', fontsize=6.5, ha='center')
    ax.text(115, 4.5, 'market stalls', fontsize=6.5, ha='center')
    ax.text(35, 11.7, 'market restaurants', fontsize=6.5, ha='center')
    ax.text(75, 6.5, 'Grand Market nave\n(the promenade)', fontsize=6.5, ha='center')
    scale_bar(ax, 110, -52, [0, 5, 10, 20], label='m')
    ax.set_xlim(-62, 170)
    ax.set_ylim(-58, 44)
    title_block(fig, 'Ring 0 · cross-section through the Grand Market (Quarter 1, s = 400 m)',
                'Cut between two frame nodes, looking along the ring. Levels above ring 0\'s floor, 3,063 m above '
                'the summit.')
    save(fig, 'ring0_section')


def save(fig, name):
    fig.savefig(OUT / f'{name}.svg')
    fig.savefig(OUT / f'{name}.png', dpi=110)
    plt.close(fig)
    print('wrote', name)


if __name__ == '__main__':
    key_plan()
    hall_level()
    roof_level()
    section()
