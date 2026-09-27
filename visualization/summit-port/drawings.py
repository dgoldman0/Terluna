"""First-pass drawings of the summit port as SVG sheets, from the summit tower study's results: the stacking plan
(used by page.py), elevations, plans, band details and the crown.

    python3 visualization/summit-port/drawings.py   # writes build/sheets/*.svg
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
TOWER = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'
R = json.loads(TOWER.read_text())
PORT = R['port']
D = PORT['design']
H = D['height_km'] * 1e3
B0, BT, P = D['base_width_m'], 320.0, D['flare_exponent']
GROUND_ASL = R['site']['ground_above_sea_m']
FLIGHT_LO = 35000.0 - GROUND_ASL
GOLDEN = 137.50776
BAND, DEPTH, CAP, SHARE = 200.0, 40.0, 100_000.0, 0.02

# Blueprint colours.
BG, GRID, INK, FAINT, DIM = '#0f3a5c', '#1a4a70', '#e2eef7', '#7fa6c6', '#4f7ea6'
AMBER, AMBER_DIM, CYAN, GREEN = '#ffd37a', '#b9955a', '#9fe3ff', '#a8e6a1'
MONO = "'IBM Plex Mono', ui-monospace, SFMono-Regular, Menlo, monospace"
COND = "'Barlow Condensed', 'Arial Narrow', sans-serif-condensed, sans-serif"


def width(z):
    return BT + (B0 - BT) * max(0.0, 1.0 - z / H) ** P


def storey(z):
    return min(CAP, SHARE * width(z) ** 2)


def f(x):
    return f'{x:.1f}'


class Svg:
    def __init__(self, w, h, label):
        self.w, self.h, self.label, self.parts = w, h, label, []

    def add(self, s):
        self.parts.append(s)

    def line(self, x1, y1, x2, y2, stroke=INK, width=1.0, dash=None, opacity=1.0, cap='butt'):
        d = f' stroke-dasharray="{dash}"' if dash else ''
        self.add(f'<line x1="{f(x1)}" y1="{f(y1)}" x2="{f(x2)}" y2="{f(y2)}" stroke="{stroke}" '
                 f'stroke-width="{width}" stroke-opacity="{opacity}" stroke-linecap="{cap}"{d}/>')

    def poly(self, pts, stroke=INK, width=1.0, fill='none', opacity=1.0, dash=None, close=False, fill_opacity=1.0,
             join='round'):
        if len(pts) < 2:
            return
        tag = 'polygon' if close else 'polyline'
        d = f' stroke-dasharray="{dash}"' if dash else ''
        pts_s = ' '.join(f'{f(x)},{f(y)}' for x, y in pts)
        self.add(f'<{tag} points="{pts_s}" fill="{fill}" fill-opacity="{fill_opacity}" stroke="{stroke}" '
                 f'stroke-width="{width}" stroke-opacity="{opacity}" stroke-linejoin="{join}"{d}/>')

    def rect(self, x, y, w, h, fill='none', stroke='none', width=1.0, opacity=1.0, dash=None, rx=0):
        d = f' stroke-dasharray="{dash}"' if dash else ''
        self.add(f'<rect x="{f(x)}" y="{f(y)}" width="{f(w)}" height="{f(h)}" rx="{rx}" fill="{fill}" '
                 f'fill-opacity="{opacity}" stroke="{stroke}" stroke-width="{width}"{d}/>')

    def circle(self, x, y, r, fill='none', stroke=INK, width=1.0, opacity=1.0):
        self.add(f'<circle cx="{f(x)}" cy="{f(y)}" r="{f(r)}" fill="{fill}" stroke="{stroke}" '
                 f'stroke-width="{width}" opacity="{opacity}"/>')

    def text(self, x, y, s, size=12, fill=INK, anchor='start', family=MONO, weight=400, spacing=0.0, opacity=1.0,
             rotate=None):
        s = s.replace('&', '&amp;').replace('<', '&lt;')
        t = f' transform="rotate({rotate} {f(x)} {f(y)})"' if rotate is not None else ''
        self.add(f'<text x="{f(x)}" y="{f(y)}" font-family="{family}" font-size="{size}" font-weight="{weight}" '
                 f'letter-spacing="{spacing}" fill="{fill}" fill-opacity="{opacity}" text-anchor="{anchor}"{t}>{s}</text>')

    def svg(self):
        grid = []
        for x in range(0, int(self.w) + 1, 40):
            grid.append(f'<line x1="{x}" y1="0" x2="{x}" y2="{self.h}" stroke="{GRID}" stroke-width="0.6"/>')
        for y in range(0, int(self.h) + 1, 40):
            grid.append(f'<line x1="0" y1="{y}" x2="{self.w}" y2="{y}" stroke="{GRID}" stroke-width="0.6"/>')
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" role="img" '
                f'aria-label="{self.label}"><title>{self.label}</title>'
                f'<rect width="{self.w}" height="{self.h}" fill="{BG}"/>' + ''.join(grid) +
                f'<rect x="8" y="8" width="{self.w - 16}" height="{self.h - 16}" fill="none" stroke="{FAINT}" '
                f'stroke-width="1.2"/>' + ''.join(self.parts) + '</svg>')


def title_block(s, sheet, name, note):
    x, y, w, h = s.w - 430, s.h - 70, 414, 54
    s.rect(x, y, w, h, fill=BG, stroke=FAINT, width=1.2)
    s.line(x + 88, y, x + 88, y + h, stroke=FAINT, width=1)
    s.text(x + 44, y + 34, sheet, 20, anchor='middle', family=COND, weight=600, spacing=1)
    s.text(x + 100, y + 22, name.upper(), 17, family=COND, weight=600, spacing=1.5)
    s.text(x + 100, y + 42, note, 10.5, fill=FAINT)


# Geometry of the two frames. Both use the frame model's width: across the corners for the hexagon, the diameter for
# the round diagrid.
def hexagon_point(radius, s):
    """Point at arc length s (m) along a hexagon of circumradius `radius`, starting at the vertex at angle 0 and
    running anticlockwise."""
    side = radius
    s %= 6 * side
    k = int(s // side)
    t = (s - k * side) / side
    a0, a1 = math.radians(60 * k), math.radians(60 * (k + 1))
    x0, y0 = radius * math.cos(a0), radius * math.sin(a0)
    x1, y1 = radius * math.cos(a1), radius * math.sin(a1)
    return x0 + (x1 - x0) * t, y0 + (y1 - y0) * t


def hexagon_arc_at(radius, angle_deg):
    """Arc length along the hexagon to the point at this polar angle."""
    a = angle_deg % 360.0
    k = int(a // 60)
    a0 = math.radians(60 * k)
    phi = math.radians(a) - a0
    # distance along the side from vertex k, by the sine rule in the triangle (centre, vertex k, point)
    t = math.sin(phi) / math.sin(math.radians(120) - phi)
    return (k + t) * radius


def strip_points(shape, radius, centre_deg, length, n=60):
    if shape == 'hex':
        s0 = hexagon_arc_at(radius, centre_deg)
        return [hexagon_point(radius, s0 - length / 2 + length * i / n) for i in range(n + 1)]
    half = length / (2 * radius)
    c = math.radians(centre_deg)
    return [(radius * math.cos(c - half + 2 * half * i / n), radius * math.sin(c - half + 2 * half * i / n))
            for i in range(n + 1)]


def split_near_far(points):
    """Runs of points on the near (y < 0) and far sides, for an elevation seen from -y."""
    runs, cur, side = [], [], None
    for x, y in points:
        s = y < 0
        if side is None or s == side:
            cur.append((x, y))
        else:
            runs.append((side, cur))
            cur = [cur[-1], (x, y)]
        side = s
    if cur:
        runs.append((side, cur))
    return runs


# ---------------------------------------------------------------------------------------------------------- sheet 1
def elevations():
    s = Svg(1600, 1340, 'Elevations of the summit port: a hexagonal lattice and a round diagrid, to scale, with Mount Everest')
    S = 44.0 / 1000.0            # px per m
    left, ground_y = 78.0, 1238.0

    def X(xkm):
        return left + xkm * 1000 * S

    def Y(z):
        return ground_y - z * S

    top_z = 26600.0
    # flight band
    s.rect(left, Y(top_z), 1600 - left - 20, Y(FLIGHT_LO) - Y(top_z), fill=CYAN, opacity=0.07)
    s.line(left, Y(FLIGHT_LO), 1580, Y(FLIGHT_LO), stroke=CYAN, width=1, dash='6 4', opacity=0.8)
    s.text(X(26.9) + 14, Y(FLIGHT_LO) - 50, 'FLIGHT BAND', 15, fill=CYAN, family=COND, weight=600, spacing=2)
    s.text(X(26.9) + 14, Y(FLIGHT_LO) - 34, 'from 35 km above sea level,', 10.5, fill=CYAN)
    s.text(X(26.9) + 14, Y(FLIGHT_LO) - 20, '23.4 km above the summit', 10.5, fill=CYAN)
    # height axis
    for km in range(0, 27, 2):
        y = Y(km * 1000)
        s.line(left - 8, y, left, y, stroke=FAINT)
        s.text(left - 12, y + 4, f'{km}', 11, fill=FAINT, anchor='end')
    s.text(30, Y(13000), 'KM ABOVE THE SUMMIT GROUND (11.6 KM ABOVE SEA LEVEL)', 12, fill=FAINT, anchor='middle',
           family=COND, weight=600, spacing=1.5, rotate=-90)
    # ground
    s.line(left, ground_y, 1580, ground_y, stroke=INK, width=1.6)
    for i in range(0, 76):
        x = left + i * 20
        s.line(x, ground_y, x - 8, ground_y + 8, stroke=FAINT, width=0.8)

    zones = [(0, 3000, 'OPEN FLOORS'), (3000, 10000, 'ENCLOSED'), (10000, 15000, 'ENCLOSED'), (15000, 24000, 'SEALED')]
    zs = PORT['zones']

    def tower(cx_km, shape, label, sub):
        cx = X(cx_km)
        zz = [H * i / 240 for i in range(241)]
        # zone tints
        for (lo, hi, _), op in zip(zones, (0.10, 0.05, 0.08, 0.12)):
            pts = [(cx - width(z) / 2 * S, Y(z)) for z in zz if lo <= z <= hi] + \
                  [(cx + width(z) / 2 * S, Y(z)) for z in reversed(zz) if lo <= z <= hi]
            s.poly(pts, stroke='none', fill=INK, fill_opacity=op, close=True)
        if shape == 'hex':
            for frac, w_, op in ((1.0, 1.6, 1.0), (0.5, 1.2, 0.9), (-0.5, 1.2, 0.9), (-1.0, 1.6, 1.0)):
                s.poly([(cx + frac * width(z) / 2 * S, Y(z)) for z in zz], width=w_, opacity=op)
            # braced panels: front face between the legs at +-R/2, side faces foreshortened
            z = 0.0
            while z < H - 1:
                r = width(z) / 2
                z2 = min(H, z + r * 0.9)
                r2 = width(z2) / 2
                s.line(cx - r / 2 * S, Y(z), cx + r2 / 2 * S, Y(z2), width=0.9, opacity=0.85)
                s.line(cx + r / 2 * S, Y(z), cx - r2 / 2 * S, Y(z2), width=0.9, opacity=0.85)
                for sgn in (-1, 1):
                    s.line(cx + sgn * r / 2 * S, Y(z), cx + sgn * r2 * S, Y(z2), stroke=FAINT, width=0.6, opacity=0.7)
                    s.line(cx + sgn * r * S, Y(z), cx + sgn * r2 / 2 * S, Y(z2), stroke=FAINT, width=0.6, opacity=0.7)
                s.line(cx - r2 * S, Y(z2), cx + r2 * S, Y(z2), stroke=FAINT, width=0.6, opacity=0.8)
                z = z2
        else:
            # round diagrid: 24 members each way; the spiral relaxes towards the crown (a constant twist rate, so the
            # members lean 35 degrees at the base and stand nearly upright at the top)
            n, c = 24, math.tan(math.radians(35)) / (B0 / 2)
            for sgn in (1, -1):
                for i in range(n):
                    pts3 = []
                    for z in zz:
                        a = 2 * math.pi * i / n + sgn * c * z
                        r = width(z) / 2
                        pts3.append((r * math.cos(a), r * math.sin(a), z))
                    run, side = [], None
                    for x, y, z in pts3:
                        near = y < 0
                        if side is not None and near != side:
                            s.poly(run, width=0.9 if side else 0.5, stroke=INK if side else DIM, opacity=0.9 if side else 0.6)
                            run = run[-1:]
                        run.append((cx + x * S, Y(z)))
                        side = near
                    s.poly(run, width=0.9 if side else 0.5, stroke=INK if side else DIM, opacity=0.9 if side else 0.6)
            # six lift spines along the meridians
            for k in range(6):
                a = math.radians(60 * k + 30)
                near = math.sin(a) < 0
                s.poly([(cx + width(z) / 2 * math.cos(a) * S, Y(z)) for z in zz], width=1.6 if near else 1.0,
                       stroke=CYAN if near else DIM, opacity=0.9 if near else 0.7)
            s.poly([(cx - width(z) / 2 * S, Y(z)) for z in zz], width=1.4)
            s.poly([(cx + width(z) / 2 * S, Y(z)) for z in zz], width=1.4)
            for zr in (3000, 6000, 9000, 10000, 12000, 15000, 18000, 21000):
                s.line(cx - width(zr) / 2 * S, Y(zr), cx + width(zr) / 2 * S, Y(zr), stroke=FAINT, width=0.8)
        # floor bands: each a 40 m strip, set 137.5 degrees round from the band below
        for k in range(1, int(H // BAND)):
            z = k * BAND
            r = width(z) / 2
            pts = strip_points(shape, r, k * GOLDEN, storey(z) / DEPTH)
            hub = k % 3 == 0
            for near, run in split_near_far(pts):
                s.poly([(cx + x * S, Y(z)) for x, y in run], stroke=AMBER if near else AMBER_DIM,
                       width=(2.4 if hub else 1.8) if near else 1.0, opacity=1.0 if near else 0.7, join='miter')
        # crown: docking arms in the flight band
        for j in range(8):
            z = FLIGHT_LO + 60 + j * 65
            a = math.radians(j * GOLDEN + 20)
            length = 900 + 180 * (j % 3)
            r0 = width(z) / 2
            x0, x1 = r0 * math.cos(a), (r0 + length) * math.cos(a)
            near = math.sin(a) < 0
            s.line(cx + x0 * S, Y(z), cx + x1 * S, Y(z), stroke=CYAN if near else DIM, width=1.6 if near else 1.0)
            if near:
                s.add(f'<ellipse cx="{f(cx + x1 * S)}" cy="{f(Y(z) - 3)}" rx="6" ry="2.4" fill="{CYAN}" fill-opacity="0.9"/>')
        s.text(cx, Y(H) - 64, label, 18, anchor='middle', family=COND, weight=600, spacing=2)
        s.text(cx, Y(H) - 48, sub, 10.5, fill=FAINT, anchor='middle')

    tower(4.95, 'hex', 'A  HEXAGONAL LATTICE', 'six straight legs, braced faces')
    tower(14.6, 'round', 'B  ROUND DIAGRID', 'two opposite spirals, six lift spines')
    # dimension: base width
    for cx_km in (4.95, 14.6):
        y = ground_y + 30
        s.line(X(cx_km) - B0 / 2 * S, y, X(cx_km) + B0 / 2 * S, y, stroke=FAINT, width=1)
        for sgn in (-1, 1):
            s.line(X(cx_km) + sgn * B0 / 2 * S, y - 5, X(cx_km) + sgn * B0 / 2 * S, y + 5, stroke=FAINT)
        s.text(X(cx_km), y + 18, f'{B0 / 1000:.1f} km at the base, narrowing to 320 m', 11, fill=FAINT, anchor='middle')
    # Everest and the Burj Khalifa, to the same scale
    ex = 22.6
    peak = [(-3.8, 0), (-3.0, 1.4), (-2.2, 2.6), (-1.6, 4.4), (-0.9, 6.3), (-0.4, 7.9), (0, 8.85), (0.5, 8.2), (1.1, 7.4),
            (1.7, 5.9), (2.4, 4.1), (3.2, 2.2), (4.0, 0)]
    s.poly([(X(ex + a), Y(b * 1000)) for a, b in peak], stroke=FAINT, width=1.2, fill=INK, fill_opacity=0.06, close=True)
    s.text(X(ex), Y(8850) - 12, 'MOUNT EVEREST, 8.85 KM', 13, fill=FAINT, anchor='middle', family=COND, weight=600, spacing=1.5)
    bx = 9.8
    s.line(X(bx), ground_y, X(bx), Y(828), stroke=INK, width=2.2)
    s.text(X(bx), Y(828) - 20, 'Burj Khalifa', 10, fill=FAINT, anchor='middle')
    s.text(X(bx), Y(828) - 8, '0.83 km', 10, fill=FAINT, anchor='middle')
    # zone labels
    lx = X(26.9)
    for (lo, hi, name), z in zip(zones, zs):
        yl, yh = Y(lo), Y(hi)
        s.line(lx, yl - 2, lx, yh + 2, stroke=FAINT, width=1)
        s.line(lx, yl - 2, lx + 6, yl - 2, stroke=FAINT)
        s.line(lx, yh + 2, lx + 6, yh + 2, stroke=FAINT)
        ym = (yl + yh) / 2
        busy = z['present_at_busy_hour']
        uses = z['uses_km2']
        top3 = sorted(uses.items(), key=lambda kv: -kv[1])[:3]
        short = {'travel: terminals, docks, customs, lounges, ship servicing': 'travel',
                 'hotels': 'hotels', 'short stay: compact and transit rooms': 'short stay',
                 'commerce and trade': 'commerce', 'shops, markets, food and drink': 'shops and food',
                 'events, conferences, leisure, observation': 'events', 'services and back of house': 'services'}
        s.text(lx + 14, ym - 22, f'{name}  {lo / 1000:g}–{hi / 1000:g} KM', 14, family=COND, weight=600, spacing=1.5)
        air = z['air_at_top']
        s.text(lx + 14, ym - 6, f"{z['floor_km2']:.1f} km² · {air['pressure_atm']:.2f} atm at the top", 10.5, fill=FAINT)
        s.text(lx + 14, ym + 8, f"air like {air['oxygen_like_earth_at_m']:,.0f} m on Earth", 10.5, fill=FAINT)
        s.text(lx + 14, ym + 22, ', '.join(f'{short[k]} {v:.1f}' for k, v in top3), 10.5, fill=AMBER)
        s.text(lx + 14, ym + 36, f'{busy[0] / 1e3:,.0f}–{busy[1] / 1e3:,.0f} thousand at a busy hour', 10.5, fill=FAINT)
    # legend
    lg_y = 70
    s.line(X(26.9), lg_y, X(26.9) + 26, lg_y, stroke=AMBER, width=2.2)
    s.text(X(26.9) + 34, lg_y + 4, 'floor band on the near side (hubs bolder)', 10.5)
    s.line(X(26.9), lg_y + 18, X(26.9) + 26, lg_y + 18, stroke=AMBER_DIM, width=1)
    s.text(X(26.9) + 34, lg_y + 22, 'floor band on the far side', 10.5)
    s.line(X(26.9), lg_y + 36, X(26.9) + 26, lg_y + 36, stroke=CYAN, width=1.6)
    s.text(X(26.9) + 34, lg_y + 40, 'docking arms, lift spines', 10.5)
    # scale bar
    sb_y = ground_y + 62
    for i in range(6):
        s.rect(X(0.4 + i), sb_y, 1000 * S, 6, fill=INK if i % 2 == 0 else BG, stroke=INK, width=0.8)
    s.text(X(0.4), sb_y + 20, '0', 10, fill=FAINT, anchor='middle')
    s.text(X(6.4), sb_y + 20, '6 km', 10, fill=FAINT, anchor='middle')
    title_block(s, 'A-201', 'Elevations', 'to scale; 120 bands of four storeys every 200 m')
    return s.svg()


# ---------------------------------------------------------------------------------------------------------- sheet 2
def plans():
    heights = [1500.0, 6500.0, 12500.0, 19500.0]
    s = Svg(1600, 880, 'Plans of both frames at four heights, with the floor bands turning 137.5 degrees each')
    pw, ph, x0, y0 = 330.0, 330.0, 26.0, 56.0
    rows = [('hex', 'A  HEXAGONAL LATTICE'), ('round', 'B  ROUND DIAGRID')]
    for ri, (shape, name) in enumerate(rows):
        s.text(x0, y0 + ri * (ph + 22) - 14, name, 16, family=COND, weight=600, spacing=2)
        for ci, z in enumerate(heights):
            px, py = x0 + ci * (pw + 10), y0 + ri * (ph + 22)
            s.rect(px, py, pw, ph, stroke=DIM, width=0.8)
            cx, cy = px + pw / 2, py + 64 + (ph - 100) / 2
            wz = width(z)
            r = wz / 2
            sc = (ph - 128) / wz
            k = int(round(z / BAND))

            def P_(x, y):
                return cx + x * sc, cy - y * sc

            if shape == 'hex':
                verts = [P_(r * math.cos(math.radians(60 * i)), r * math.sin(math.radians(60 * i))) for i in range(6)]
                s.poly(verts, close=True, width=1.4)
                for vx, vy in verts:
                    s.rect(vx - 4, vy - 4, 8, 8, fill=CYAN, stroke=BG, width=0.8)
                # face bracing, one X per face
                for i in range(6):
                    a0, a1 = verts[i], verts[(i + 1) % 6]
                    mx, my = (a0[0] + a1[0]) / 2, (a0[1] + a1[1]) / 2
                    s.line(mx, my, cx + (mx - cx) * 0.93, cy + (my - cy) * 0.93, stroke=FAINT, width=0.6)
            else:
                s.circle(cx, cy, r * sc, width=1.4)
                n = 24
                for i in range(n):
                    a = 2 * math.pi * i / n
                    s.line(cx + r * sc * math.cos(a) * 0.96, cy - r * sc * math.sin(a) * 0.96,
                           cx + r * sc * math.cos(a) * 1.04, cy - r * sc * math.sin(a) * 1.04, stroke=FAINT, width=0.8)
                for i in range(6):
                    a = math.radians(60 * i + 30)
                    x, y = P_(r * math.cos(a), r * math.sin(a))
                    s.rect(x - 4, y - 4, 8, 8, fill=CYAN, stroke=BG, width=0.8)
            # the next bands above, faint, numbered by how many bands up
            for j in range(7, 0, -1):
                pts = strip_points(shape, r, (k + j) * GOLDEN, storey(z) / DEPTH)
                s.poly([P_(x, y) for x, y in pts], stroke=AMBER_DIM, width=3.2, opacity=0.55, join='miter')
                mx, my = pts[len(pts) // 2]
                lx, ly = P_(mx * 1.13, my * 1.13)
                s.text(lx, ly + 4, f'+{j}', 10, fill=AMBER_DIM, anchor='middle')
            pts = strip_points(shape, r, k * GOLDEN, storey(z) / DEPTH)
            s.poly([P_(x, y) for x, y in pts], stroke=AMBER, width=5, join='miter')
            # berths on the outer edge of this band, every 60 m
            length = storey(z) / DEPTH
            nb = int(length // 60)
            for b in range(nb):
                frac = (b + 0.5) / nb
                i = frac * (len(pts) - 1)
                i0 = int(i)
                x, y = pts[i0]
                rr = math.hypot(x, y)
                ox, oy = x / rr, y / rr
                a1, b1 = P_(x + ox * 0.012 * wz, y + oy * 0.012 * wz)
                a2, b2 = P_(x + ox * 0.03 * wz, y + oy * 0.03 * wz)
                s.line(a1, b1, a2, b2, stroke=CYAN, width=1)
            # scale bar
            bar = 10 ** math.floor(math.log10(wz / 3))
            if wz / bar > 6:
                bar *= 2
            s.rect(px + 14, py + ph - 26, bar * sc, 5, fill=INK)
            s.text(px + 14 + bar * sc + 6, py + ph - 20, f'{bar / 1000:g} km', 10, fill=FAINT)
            s.text(px + 12, py + 20, f'{z / 1000:g} KM UP', 14, family=COND, weight=600, spacing=1.5)
            s.text(px + 12, py + 36, f'frame {wz / 1000:.1f} km across', 10.5, fill=FAINT)
            st = storey(z)
            s.text(px + pw - 12, py + 20, f'storey {st:,.0f} m²', 10.5, fill=AMBER, anchor='end')
            s.text(px + pw - 12, py + 34, f'strip {length / 1000:.2f} km x 40 m', 10.5, fill=FAINT, anchor='end')
            s.text(px + pw - 12, py + 48, f'{st / 2000:.0f} compartments, {nb} berths', 10.5, fill=FAINT, anchor='end')
    note_x = x0 + 4 * (pw + 10) + 14
    s.text(note_x, 70, 'READING THE PLANS', 14, family=COND, weight=600, spacing=1.5)
    lines = ['Bold amber: the band at this', 'height, a 40 m strip drawn', 'thicker than scale.', '',
             'Faint amber, +1 to +7: the', 'bands above, each turned', '137.5° (the golden angle,', 'as leaves round a stem).', '',
             'Cyan squares: legs (A) or', 'lift spines (B), carrying', 'lifts, stairs and risers.', '',
             'Cyan ticks: berths, one', 'every 60 m of outer rim.', '',
             'Inside the frame: open air.', 'Every strip has open air', 'on both faces.']
    for i, t in enumerate(lines):
        s.text(note_x, 92 + i * 15, t, 10.5, fill=FAINT if t else FAINT)
    title_block(s, 'A-101', 'Plans', 'each plan to its own scale bar')
    return s.svg()


# ---------------------------------------------------------------------------------------------------------- sheet 3
def band_detail():
    s = Svg(1600, 1000, 'A stretch of one floor band in plan and in section')
    k = 4.2                      # px per m, plan
    x0, y0 = 60.0, 214.0
    L, Dp = 300.0, DEPTH
    s.text(x0, 64, 'PLAN OF 300 M OF A FLOOR BAND', 16, family=COND, weight=600, spacing=2)
    s.text(x0, 82, 'six compartments of 2,000 m² (50 x 40 m) along a strip with open air on both faces', 10.5, fill=FAINT)
    # outer edge at the top (y0), inner at the bottom
    s.text(x0 + L * k + 14, y0 - 40, 'OUTSIDE THE FRAME', 12, fill=CYAN, family=COND, weight=600, spacing=1.5)
    s.text(x0 + L * k + 14, y0 - 24, 'berths every 60 m', 10.5, fill=CYAN)
    s.text(x0 + L * k + 14, y0 + Dp * k + 44, 'OPEN INTERIOR OF THE FRAME', 12, fill=CYAN, family=COND, weight=600,
           spacing=1.5)
    # gardens along both rims (5 m)
    for yy in (y0, y0 + (Dp - 5) * k):
        s.rect(x0, yy, L * k, 5 * k, fill=GREEN, opacity=0.14, stroke=GREEN, width=0.8)
        for i in range(0, int(L), 6):
            s.circle(x0 + (i + 3) * k, yy + 2.5 * k, 1.6 * k / 2, fill=GREEN, stroke='none', opacity=0.55)
    s.rect(x0, y0, L * k, Dp * k, stroke=INK, width=1.6)
    # spine corridor along the middle
    cy1, cy2 = y0 + 18.5 * k, y0 + 21.5 * k
    s.line(x0, cy1, x0 + L * k, cy1, stroke=FAINT, width=0.8)
    s.line(x0, cy2, x0 + L * k, cy2, stroke=FAINT, width=0.8)
    uses = ['HOTEL ROOMS', 'HOTEL ROOMS', 'LEG AND CORE', 'COMMERCE', 'SHOPS AND FOOD', 'EVENTS HALL']
    for i in range(6):
        cx0 = x0 + i * 50 * k
        if i > 0:
            # fire wall with double doors on the spine
            s.line(cx0, y0 + 5 * k, cx0, cy1, stroke=INK, width=3)
            s.line(cx0, cy2, cx0, y0 + (Dp - 5) * k, stroke=INK, width=3)
            s.add(f'<path d="M {f(cx0)} {f(cy1)} a {f(3 * k)} {f(3 * k)} 0 0 1 {f(3 * k)} {f(3 * k)}" fill="none" '
                  f'stroke="{INK}" stroke-width="0.8"/>')
        name = uses[i]
        if name == 'HOTEL ROOMS':
            for side_y0, side_y1 in ((y0 + 5 * k, cy1), (cy2, y0 + (Dp - 5) * k)):
                for rj in range(10):
                    rx = cx0 + rj * 5 * k
                    s.rect(rx, side_y0, 5 * k, side_y1 - side_y0, stroke=FAINT, width=0.6)
            s.text(cx0 + 25 * k, y0 + 20.4 * k, 'HOTEL ROOMS, 5 x 13.5 M', 10, anchor='middle', fill=INK)
        elif name == 'LEG AND CORE':
            lx0 = cx0 + 10 * k
            s.rect(lx0, y0 + 5 * k, 30 * k, 30 * k, fill=INK, opacity=0.06, stroke=CYAN, width=1.6)
            for a in range(4):
                for b in range(2):
                    s.rect(lx0 + (3 + a * 6) * k, y0 + (8 + b * 7) * k, 5 * k, 5 * k, stroke=CYAN, width=0.9)
            s.rect(lx0 + 3 * k, y0 + 24 * k, 11 * k, 9 * k, stroke=INK, width=0.9)
            s.rect(lx0 + 16 * k, y0 + 24 * k, 11 * k, 9 * k, stroke=INK, width=0.9)
            for t in range(8):
                s.line(lx0 + (3 + t * 1.4) * k, y0 + 24 * k, lx0 + (3 + t * 1.4) * k, y0 + 33 * k, stroke=FAINT, width=0.5)
                s.line(lx0 + (16 + t * 1.4) * k, y0 + 24 * k, lx0 + (16 + t * 1.4) * k, y0 + 33 * k, stroke=FAINT, width=0.5)
            s.text(cx0 + 25 * k, y0 + 44 * k + 16, 'leg: 8 lifts, 2 stairs, risers, pressurised lobby', 10, anchor='middle', fill=CYAN)
        else:
            s.text(cx0 + 25 * k, y0 + 13 * k, name, 12, anchor='middle', family=COND, weight=600, spacing=1.5)
            s.text(cx0 + 25 * k, y0 + 13 * k + 14, 'open plan, sprinklered', 10, anchor='middle', fill=FAINT)
            s.text(cx0 + 25 * k, y0 + 28 * k, '2,000 m²', 10, anchor='middle', fill=FAINT)
    # berths every 60 m on the outer rim, with ships
    for b in range(5):
        bx = x0 + (30 + b * 60) * k
        s.line(bx, y0, bx, y0 - 8 * k, stroke=CYAN, width=1.4)
        s.add(f'<ellipse cx="{f(bx)}" cy="{f(y0 - 18 * k)}" rx="{f(22 * k)}" ry="{f(6 * k)}" fill="none" stroke="{CYAN}" stroke-width="1.2"/>')
        s.text(bx, y0 - 18 * k + 4, 'air taxi' if b % 2 == 0 else 'regional ship', 9.5, fill=CYAN, anchor='middle')
    # dimensions
    yd = y0 + Dp * k + 16
    s.line(x0, yd, x0 + 50 * k, yd, stroke=FAINT)
    s.text(x0 + 25 * k, yd + 14, '50 m', 10, fill=FAINT, anchor='middle')
    xd = x0 - 16
    s.line(xd, y0, xd, y0 + Dp * k, stroke=FAINT)
    s.text(xd - 4, y0 + Dp * k / 2, '40 m', 10, fill=FAINT, anchor='end')

    # section across the band
    ks = 7.0
    sx0, sy0 = 90.0, 610.0
    s.text(60, 520, 'SECTION ACROSS A BAND', 16, family=COND, weight=600, spacing=2)
    s.text(60, 538, 'four storeys of 4 m; the next band is 184 m of open frame above, turned 137.5°', 10.5, fill=FAINT)
    base = sy0 + 16 * ks
    for st in range(5):
        y = base - st * 4 * ks
        s.line(sx0, y, sx0 + 40 * ks, y, stroke=INK, width=2 if st in (0, 4) else 1)
    s.line(sx0, base, sx0, base - 16 * ks, stroke=INK, width=1.2)
    s.line(sx0 + 40 * ks, base, sx0 + 40 * ks, base - 16 * ks, stroke=INK, width=1.2)
    for yy in range(4):
        y = base - yy * 4 * ks
        for i in range(8):
            s.circle(sx0 + (2.5 + i * 5) * ks, y - 3.6 * ks, 1.2, fill=CYAN, stroke='none')
    # gardens at the rims on the top storey terraces
    for gx in (sx0 - 6 * ks, sx0 + 40 * ks):
        s.rect(gx, base - 16 * ks - 1.2 * ks, 6 * ks, 1.2 * ks, fill=GREEN, opacity=0.35, stroke=GREEN, width=0.8)
        for i in range(3):
            s.circle(gx + (1 + i * 2) * ks, base - 16 * ks - 3 * ks, 1.3 * ks, fill=GREEN, stroke='none', opacity=0.5)
    # face lattice member through the band
    s.line(sx0 + 38.5 * ks, base + 38, sx0 + 42.5 * ks, base - 16 * ks - 38, stroke=FAINT, width=4, opacity=0.6)
    s.text(sx0 + 38.5 * ks - 6, base + 34, 'face lattice', 10, fill=FAINT, anchor='end')
    # ship at the outer berth
    shx = sx0 + 40 * ks + 16 * ks
    s.add(f'<ellipse cx="{f(shx + 14 * ks)}" cy="{f(base - 9 * ks)}" rx="{f(16 * ks)}" ry="{f(4.5 * ks)}" fill="none" stroke="{CYAN}" stroke-width="1.4"/>')
    s.line(sx0 + 40 * ks, base - 8 * ks, shx - 2 * ks, base - 8 * ks, stroke=CYAN, width=2)
    s.text(shx + 14 * ks, base - 9 * ks + 4, 'ship at berth', 10, fill=CYAN, anchor='middle')
    # break lines to the next band
    for yy, lab in ((base - 16 * ks - 40, 'next band, 200 m up'), (base + 40, 'band below, 200 m down')):
        s.line(sx0 - 10, yy, sx0 + 50 * ks, yy, stroke=FAINT, dash='8 5')
        s.text(sx0 + 50 * ks + 8, yy + 4, lab, 10, fill=FAINT)
    s.text(sx0 + 20 * ks, base + 18, 'inner rim', 10, fill=FAINT, anchor='middle')

    # a hub band, 80 m across
    hx0, hy0 = 820.0, 610.0
    s.text(790, 520, 'HUB BAND, EVERY THIRD BAND (600 M)', 16, family=COND, weight=600, spacing=2)
    s.text(790, 538, 'shafts break and step; water is stored and pumped; every fifth hub keeps a fire station', 10.5, fill=FAINT)
    hb = hy0 + 16 * ks
    for st in range(5):
        y = hb - st * 4 * ks
        s.line(hx0, y, hx0 + 80 * ks, y, stroke=INK, width=2 if st in (0, 4) else 1)
    s.line(hx0, hb, hx0, hb - 16 * ks, stroke=INK, width=1.2)
    s.line(hx0 + 80 * ks, hb, hx0 + 80 * ks, hb - 16 * ks, stroke=INK, width=1.2)
    s.rect(hx0 + 4 * ks, hb - 8 * ks + 1, 26 * ks, 8 * ks - 2, fill=CYAN, opacity=0.18, stroke=CYAN, width=1)
    s.text(hx0 + 17 * ks, hb - 5.2 * ks, 'water, 186 m³', 10, fill=CYAN, anchor='middle')
    for i in range(3):
        s.circle(hx0 + (36 + i * 4) * ks, hb - 2 * ks, 1.5 * ks, stroke=CYAN, width=1)
    s.text(hx0 + 40 * ks, hb - 5.2 * ks, 'pumps', 10, fill=CYAN, anchor='middle')
    s.rect(hx0 + 56 * ks, hb - 8 * ks + 1, 22 * ks, 8 * ks - 2, stroke=AMBER, width=1)
    s.text(hx0 + 67 * ks, hb - 5.2 * ks, 'fire station', 10, fill=AMBER, anchor='middle')
    s.rect(hx0 + 4 * ks, hb - 16 * ks + 1, 36 * ks, 8 * ks - 2, stroke=FAINT, width=0.8, dash='4 3')
    s.text(hx0 + 22 * ks, hb - 13.2 * ks, 'stores, protected compartments', 10, fill=FAINT, anchor='middle')
    s.rect(hx0 + 50 * ks, hb - 16 * ks + 1, 26 * ks, 8 * ks - 2, stroke=INK, width=1)
    s.text(hx0 + 63 * ks, hb - 13.2 * ks, 'lift transfer lobby', 10, anchor='middle')
    s.add(f'<ellipse cx="{f(hx0 + 92 * ks)}" cy="{f(hb - 5 * ks)}" rx="{f(9 * ks)}" ry="{f(3 * ks)}" fill="none" stroke="{AMBER}" stroke-width="1.4"/>')
    s.line(hx0 + 80 * ks, hb - 4.5 * ks, hx0 + 83 * ks, hb - 4.5 * ks, stroke=AMBER, width=2)
    s.text(hx0 + 92 * ks, hb - 5 * ks + 4, 'fire ship', 10, fill=AMBER, anchor='middle')
    s.line(hx0 + 45 * ks, hb - 16 * ks - 34, hx0 + 45 * ks, hb - 16 * ks, stroke=INK, width=3)
    s.line(hx0 + 47 * ks, hb, hx0 + 47 * ks, hb + 34, stroke=INK, width=3)
    s.text(hx0 + 46 * ks, hb + 50, 'stair and lift shafts break and step here', 10, fill=FAINT, anchor='middle')
    title_block(s, 'A-301', 'Band details', 'plan 1 m = 4.2 px; sections 1 m = 7 px')
    return s.svg()


# ---------------------------------------------------------------------------------------------------------- sheet 4
def crown():
    s = Svg(1600, 900, 'The crown in the flight band: docking arms set by the golden angle, in plan and elevation')
    k = 0.21                     # px per m
    cx, cy = 430.0, 460.0
    s.text(60, 64, 'CROWN PLAN', 16, family=COND, weight=600, spacing=2)
    s.text(60, 82, 'long-haul docking arms in the flight band, each turned 137.5° from the one below', 10.5, fill=FAINT)
    r = BT / 2
    s.poly([(cx + r * k * math.cos(math.radians(60 * i)), cy - r * k * math.sin(math.radians(60 * i))) for i in range(6)],
           close=True, width=1.6, fill=INK, fill_opacity=0.08)
    for j in range(8):
        a = math.radians(j * GOLDEN + 20)
        length = 900 + 180 * (j % 3)
        x0, y0 = cx + r * k * math.cos(a), cy - r * k * math.sin(a)
        x1, y1 = cx + (r + length) * k * math.cos(a), cy - (r + length) * k * math.sin(a)
        # arm drawn 40 m wide
        nx, ny = -math.sin(a) * 20 * k, -math.cos(a) * 20 * k
        s.poly([(x0 + nx, y0 + ny), (x1 + nx, y1 + ny), (x1 - nx, y1 - ny), (x0 - nx, y0 - ny)], close=True,
               stroke=CYAN, width=1.2, fill=CYAN, fill_opacity=0.12)
        s.text(x1 + math.cos(a) * 18, y1 - math.sin(a) * 18 + 4, f'{j + 1}', 11, fill=CYAN, anchor='middle')
        # berths on the sheltered east side of each arm: long-haul ships drawn 220 m long
        for b in range(3):
            t = 0.3 + 0.3 * b
            px_m, py_m = (r + t * length) * math.cos(a), (r + t * length) * math.sin(a)
            nx, ny = -math.sin(a), math.cos(a)
            if nx < 0:
                nx, ny = -nx, -ny
            sx_, sy_ = cx + (px_m + nx * 60) * k, cy - (py_m + ny * 60) * k
            s.add(f'<ellipse cx="{f(sx_)}" cy="{f(sy_)}" rx="{f(110 * k)}" ry="{f(24 * k)}" '
                  f'transform="rotate({f(-math.degrees(a))} {f(sx_)} {f(sy_)})" fill="none" stroke="{INK}" '
                  f'stroke-width="1" stroke-opacity="0.85"/>')
    # wind
    s.line(80, 180, 190, 180, stroke=AMBER, width=2)
    s.poly([(190, 174), (202, 180), (190, 186)], close=True, stroke=AMBER, fill=AMBER)
    s.text(80, 168, 'PREVAILING WIND ALOFT, FROM THE WEST', 12, fill=AMBER, family=COND, weight=600, spacing=1.5)
    s.text(80, 202, 'steadiness 0.95 at the crown in the design run;', 10, fill=FAINT)
    s.text(80, 216, 'ships berth on the sheltered east side of each arm', 10, fill=FAINT)
    s.rect(60, 800, 500 * k, 5, fill=INK)
    s.text(60 + 500 * k + 6, 806, '500 m', 10, fill=FAINT)
    s.text(cx, cy + 4, 'crown', 10, fill=FAINT, anchor='middle')
    s.text(cx, cy + 16, '320 m', 10, fill=FAINT, anchor='middle')

    # elevation of the top of the tower
    ke = 0.21
    ex0, ey0 = 1180.0, 820.0          # centre x, y of z = 21 km
    z0 = 21000.0

    def EY(z):
        return ey0 - (z - z0) * ke

    s.text(880, 64, 'ELEVATION OF THE TOP 3.3 KM', 16, family=COND, weight=600, spacing=2)
    s.text(880, 82, 'sealed terminal bands below, arms above the flight band\'s floor', 10.5, fill=FAINT)
    zz = [z0 + i * 10 for i in range(0, 301)]
    s.poly([(ex0 - width(z) / 2 * ke, EY(z)) for z in zz], width=1.4)
    s.poly([(ex0 + width(z) / 2 * ke, EY(z)) for z in zz], width=1.4)
    s.line(870, EY(FLIGHT_LO), 1580, EY(FLIGHT_LO), stroke=CYAN, dash='6 4')
    s.text(1580 - 6, EY(FLIGHT_LO) - 6, 'flight band from 35 km above sea level', 10, fill=CYAN, anchor='end')
    s.rect(870, EY(24300), 710, EY(FLIGHT_LO) - EY(24300), fill=CYAN, opacity=0.06)
    for kb in range(int(z0 // BAND) + 1, int(H // BAND)):
        z = kb * BAND
        pts = strip_points('hex', width(z) / 2, kb * GOLDEN, storey(z) / DEPTH)
        for near, run in split_near_far(pts):
            s.poly([(ex0 + x * ke, EY(z)) for x, y in run], stroke=AMBER if near else AMBER_DIM,
                   width=3 if near else 1.4, opacity=1 if near else 0.7)
    for j in range(8):
        z = FLIGHT_LO + 60 + j * 65
        a = math.radians(j * GOLDEN + 20)
        length = 900 + 180 * (j % 3)
        r0 = width(z) / 2
        x0, x1 = r0 * math.cos(a), (r0 + length) * math.cos(a)
        near = math.sin(a) < 0
        s.line(ex0 + x0 * ke, EY(z), ex0 + x1 * ke, EY(z), stroke=CYAN if near else DIM, width=2.4 if near else 1.2)
        if near:
            s.add(f'<ellipse cx="{f(ex0 + x1 * ke - 50 * ke * (1 if x1 > 0 else -1))}" cy="{f(EY(z) - 12)}" rx="{f(110 * ke)}" ry="{f(20 * ke)}" '
                  f'fill="none" stroke="{INK}" stroke-width="1"/>')
    s.text(ex0, EY(H) - 40, '24 km, 35.6 km above sea level', 10.5, fill=FAINT, anchor='middle')
    for zt in (21000, 22000, 23000, 24000):
        s.line(880, EY(zt), 890, EY(zt), stroke=FAINT)
        s.text(894, EY(zt) + 4, f'{zt / 1000:g} km', 10, fill=FAINT)
    title_block(s, 'A-401', 'Crown', 'plan and elevation, 1 km = 210 px')
    return s.svg()


# ---------------------------------------------------------------------------------------------------------- sheet 5
def stacking():
    s = Svg(1600, 600, 'Stacking plan: floor by use in each zone, after the travel split')
    names = ['travel: terminals, docks, customs, lounges, ship servicing', 'hotels', 'short stay: compact and transit rooms',
             'commerce and trade', 'shops, markets, food and drink', 'events, conferences, leisure, observation',
             'services and back of house']
    labels = ['Travel', 'Hotels', 'Short stay', 'Commerce', 'Shops and food', 'Events', 'Services']
    colours = ['#9fe3ff', '#ffd37a', '#ffb38a', '#b8f2b0', '#f5a3c7', '#d1b3ff', '#c9d6e2']
    s.text(60, 64, 'STACKING PLAN', 16, family=COND, weight=600, spacing=2)
    s.text(60, 82, 'floor by use in each zone (km²), placed by where people arrive and what air they need', 10.5, fill=FAINT)
    k = 64.0                     # px per km2
    x0 = 300.0
    rows = list(reversed(PORT['zones']))
    for i, z in enumerate(rows):
        y = 118 + i * 88
        s.text(60, y + 20, f"{z['from_km']:g}–{z['to_km']:g} KM", 16, family=COND, weight=600, spacing=1.5)
        s.text(60, y + 36, z['use'], 10, fill=FAINT)
        s.text(60, y + 50, f"air like {z['air_at_top']['oxygen_like_earth_at_m']:,.0f} m on Earth at the top", 10, fill=FAINT)
        x = x0
        for n, lab, col in zip(names, labels, colours):
            v = z['uses_km2'][n]
            if v <= 0:
                continue
            s.rect(x, y, v * k, 46, fill=col, opacity=0.9, stroke=BG, width=1.5)
            if v * k > 34:
                s.text(x + v * k / 2, y + 20, f'{v:.1f}', 12, fill=BG, anchor='middle', weight=600)
                if v * k > 70:
                    s.text(x + v * k / 2, y + 36, lab.lower(), 10, fill=BG, anchor='middle')
            x += v * k
        busy = z['present_at_busy_hour']
        s.text(x + 12, y + 20, f"{z['floor_km2']:.1f} km²", 13, weight=600)
        s.text(x + 12, y + 36, f"{busy[0] / 1e3:,.0f}–{busy[1] / 1e3:,.0f} thousand at a busy hour", 10.5, fill=FAINT)
    # legend
    lx = x0
    for lab, col in zip(labels, colours):
        s.rect(lx, 480, 14, 14, fill=col, stroke='none')
        s.text(lx + 20, 492, lab, 11)
        lx += 20 + 9 * len(lab) + 26
    t = PORT['occupancy']['travel_by_trip']
    s.text(x0, 522, f"Travel by trip: long-haul {t['long_haul_km2']:.1f} km² at the top, regional {t['regional_km2']:.1f} km² at 10–15 km, "
           f"local {t['local_km2']:.1f} km² at the base and hubs", 11, fill=CYAN)
    title_block(s, 'A-001', 'Stacking plan', 'fitted so every use and zone adds up')
    return s.svg()


SHEETS = [
    ('A-001', 'Stacking plan', stacking,
     'Where each use sits, after the fix. Travel is split by trip: long-haul terminals fill the sealed zone under the '
     'crown, regional docks sit at 10–15 km, and local interchange is at the base and in the hub bands. Hotels sit low, '
     'where the air is best for sleep.'),
    ('A-201', 'Elevations', elevations,
     'Both frames to one scale, with Mount Everest and the Burj Khalifa beside them. Each amber line is a band of four '
     'storeys; each band sits 137.5° round from the one below, so the floors spiral up the tower while the frame '
     'stays straight. The round diagrid\'s two sets of spirals cancel each other\'s twist; they lean 35° at the base '
     'and straighten towards the crown.'),
    ('A-101', 'Plans', plans,
     'Each frame at four heights. A band\'s floor is a 40 m strip set on the frame\'s face, with open air on both '
     'sides. The numbered faint strips are the next seven bands up, each turned by the golden angle, so the bands '
     'nearest above a floor sit elsewhere round the tower and leave it sun, wind and clear air for ships.'),
    ('A-301', 'Band details', band_detail,
     'A stretch of floor band: 2,000 m² compartments with fire walls and double doors along a spine, gardens on both '
     'rims, and berths every 60 m. The legs carry the lifts and stairs. Every third band is a hub, where shafts break, '
     'water is stored and pumped, and every fifth hub keeps a fire station.'),
    ('A-401', 'Crown', crown,
     'The top of the tower reaches into the flight band. Long-haul ships berth on docking arms that branch from the '
     'crown, each turned 137.5° from the one below, on the sheltered east side of the steady westerly aloft.'),
]


def main():
    out = HERE / 'build' / 'sheets'
    out.mkdir(parents=True, exist_ok=True)
    for sheet, _, fn, _ in SHEETS:
        (out / f'{sheet}.svg').write_text(fn())
    print(out)


if __name__ == '__main__':
    main()
