"""Figures for the fabric and cut study: the screen's figures for every fabric, the trials against the estimates, and
the trial runs' final cloth seen in three dimensions.

    python -m humanities.dance.garment.figures TRIALS_DIR      # PNGs into TRIALS_DIR/figures

Previsualisation-grade plots of simulated cloth (flat-shaded triangles), for reading the study, not renders.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

from .screen import OUT as SCREEN

MOON_C, EARTH_C = '#3b6fb6', '#c8553d'


def _short(name):
    return name.replace('silk_', '').replace('_', ' ')


def screen_figure(G, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = sorted(G['fabrics'], key=lambda r: r['areal_density_kg_m2'])
    mark = {'measured': '', 'inferred': ' (i)', 'bracket': ' (b)', 'chosen': ' (c)'}
    names = [f"{_short(r['fabric'])}{mark.get(r.get('bending_status'), '')}\n{1000 * r['areal_density_kg_m2']:.0f} g/m²"
             for r in rows]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(1, 5, figsize=(17, 0.42 * len(rows) + 2.0), sharey=True)
    panels = [('bending_length_cm', 'bending length (cm):\nstiffness against weight'),
              ('face_on_fall_m_s', 'falling face-on (m/s)'),
              (('trail_deg', '2_m_s'), 'trailing below the horizontal\nin travel at 2 m/s (deg)'),
              (('swing', 'first_through_vertical_s'), 'a hanging panel let go at 30°:\ntime to the vertical (s)'),
              (('added_mass_over_fabric', 'fold'), "the air's added mass over the\nfabric's, for a 20 cm fold")]
    for k, (key, title) in enumerate(panels):
        for wn, col, off in (('moon', MOON_C, -0.2), ('earth', EARTH_C, 0.2)):
            vals = [r[wn][key] if isinstance(key, str) else r[wn][key[0]][key[1]] for r in rows]
            ax[k].barh(y + off, vals, height=0.38, color=col, label=wn)
        ax[k].set_title(title, fontsize=9)
        ax[k].grid(axis='x', lw=0.3, alpha=0.5)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(names, fontsize=7.5)
    ax[0].legend(fontsize=8, frameon=False)
    fig.suptitle('The fabric library on the Moon and on Earth (closed-form screen)', fontsize=10)
    fig.text(0.01, 0.01, 'Bending rigidity: measured unless marked (i) inferred from a measured drape, (b) a stated '
             'bracket or (c) a modelling choice; weights from measurements or suppliers (garment/sources.json).',
             fontsize=7.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out / 'screen.png', dpi=130)
    plt.close(fig)


def drape_figure(R, out, shapes, status=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rs = [r for r in R.get('drape', []) if r.get('seed', 0) == 0]
    if not rs:
        return
    status = status or {}
    mark = {'measured': '', 'inferred': ' (i)', 'bracket': ' (b)', 'chosen': ' (c)'}
    by = {}
    for r in rs:
        by.setdefault(r['fabric'], {})[r['world']] = r
    names = sorted(by, key=lambda n: by[n].get('earth', by[n].get('moon'))['drape_coefficient'])
    fig = plt.figure(figsize=(15, 7.5))
    ax = fig.add_axes([0.1, 0.1, 0.3, 0.82])
    y = np.arange(len(names))
    for wn, col, off in (('moon', MOON_C, -0.2), ('earth', EARTH_C, 0.2)):
        v = [by[n][wn]['drape_coefficient'] if wn in by[n] else np.nan for n in names]
        ax.barh(y + off, v, height=0.38, color=col, label=wn)
        for yy, n in zip(y, names):
            if wn in by[n]:
                ax.text(by[n][wn]['drape_coefficient'] + 0.01, yy + off, f"{by[n][wn]['nodes']}", fontsize=6.5,
                        va='center', color=col)
    ax.set_yticks(y)
    ax.set_yticklabels([_short(n) + mark.get(status.get(n), '') for n in names], fontsize=8)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel('drape coefficient (numbers: folds round the hem)')
    ax.legend(fontsize=8, frameon=False, loc='lower right')
    ax.set_title("Cusick's drape test, 30 cm on an 18 cm disc (solver)", fontsize=9)
    fig.text(0.01, 0.015, 'Bending rigidity measured unless marked (i) inferred from a measured drape, (b) a stated '
             'bracket or (c) a modelling choice.', fontsize=7.5)
    picks = [n for n in ('silk_chiffon', 'silk_satin_heavy', 'wool_plain', 'silk_organza')
             if (shapes / f'drape_{n}_moon_s0.npz').exists() and (shapes / f'drape_{n}_earth_s0.npz').exists()]
    for i, n in enumerate(picks):
        for j, wn in enumerate(('earth', 'moon')):
            d = np.load(shapes / f'drape_{n}_{wn}_s0.npz')
            a3 = fig.add_axes([0.42 + 0.145 * i, 0.5 - 0.44 * j, 0.15, 0.46], projection='3d')
            _mesh(a3, d['faces'], d['verts'], MOON_C if wn == 'moon' else EARTH_C, lim=0.13)
            a3.set_title(f'{_short(n)}, {wn}: {by[n][wn]["drape_coefficient"]:.2f}', fontsize=8, y=0.92)
            a3.view_init(elev=22, azim=-60)
    fig.savefig(out / 'drape.png', dpi=130)
    plt.close(fig)


def _mesh(a3, F, X, colour, lim=None):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    tri = X[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    light = np.array([0.4, -0.5, 0.75])
    light /= np.linalg.norm(light)
    shade = 0.35 + 0.65 * np.abs(n @ light)
    base = np.array(__import__('matplotlib').colors.to_rgb(colour))
    cols = np.clip(base[None] * shade[:, None] + 0.15, 0, 1)
    a3.add_collection3d(Poly3DCollection(tri, facecolors=cols, edgecolors='none', linewidths=0))
    c = X.mean(0)
    r = lim or 0.55 * (X.max(0) - X.min(0)).max()
    a3.set_xlim(c[0] - r, c[0] + r)
    a3.set_ylim(c[1] - r, c[1] + r)
    a3.set_zlim(c[2] - r, c[2] + r)
    a3.set_axis_off()


def turn_figure(G, R, out, shapes):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .materials import FABRICS
    from .screen import wind_up_index, worlds
    rs = [r for r in R.get('turn', []) if r['hem_kg_m'] == 0] + R.get('hem', [])
    if not rs:
        return
    W = worlds()
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for r in rs:
        tt = [t_ for t_ in r['turns'] if t_['rpm'] == 27.0]
        if not tt:
            continue
        mat = FABRICS[r['fabric']]
        cut = math.radians(r['cut_deg'])
        idx = wind_up_index(mat.areal_density, W[r['world']]['rho'], cut, hem_kg_m=r['hem_kg_m'],
                            factor=mat.normal_factor)
        col = MOON_C if r['world'] == 'moon' else EARTH_C
        mk = 'o' if r['cut_deg'] == 360 else '^'
        face = 'none' if r['hem_kg_m'] > 0 else col
        label = _short(r['fabric']) + (f" +{1000 * r['hem_kg_m']:.0f} g/m" if r['hem_kg_m'] > 0 else '') + \
            (' (half)' if r['cut_deg'] == 180 else '') + (', earth' if r['world'] == 'earth' else '') + \
            (' *' if r['strain']['p99_stretch'] > 0.05 else '')
        for k, key in enumerate(('flare_deg', 'lag_deg')):
            v = tt[0][key]
            ax[k].plot(idx, v, mk, color=col, mfc=face, ms=7)
            ax[k].annotate(label, (idx, v), fontsize=6.3, xytext=(4, 3), textcoords='offset points')
    ax[0].axhline(74.8, color=MOON_C, lw=0.8, ls='--')
    ax[0].text(0.02, 76, 'a full circle without air, Moon', fontsize=7, color=MOON_C, transform=ax[0].get_yaxis_transform())
    ax[0].axhline(30.0, color='k', lw=0.6, ls=':')
    ax[0].text(0.02, 31, "a half circle's smooth cone", fontsize=7, transform=ax[0].get_yaxis_transform())
    for k, (lab, title) in enumerate((('flare at 27 rpm (deg)', 'How far the skirt opens'),
                                       ("hem's lag behind the waist at 27 rpm (deg)", 'How far it winds round'))):
        ax[k].set_xscale('log')
        ax[k].set_xlabel("wind-up index: the folds' drag over the centrifugal pull (screen)")
        ax[k].set_ylabel(lab)
        ax[k].set_title(title + ' (solver; circles full, triangles half circles, open: weighted hem)', fontsize=8.5)
        ax[k].grid(lw=0.3, alpha=0.5)
    fig.text(0.01, 0.01, '* the cloth sheared beyond what a woven fabric allows (99th-percentile strain over 5%): the '
             'solver\'s shear does not lock, so the shape is not trustworthy; a lag past 180 degrees reads as negative.',
             fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out / 'turn.png', dpi=130)
    plt.close(fig)
    # the turning skirts seen from the side
    files = sorted(shapes.glob('turn_*_moon_*_rpm27.npz'), key=lambda p: ('hem0' not in p.name, p.name)) + \
        sorted(shapes.glob('turn_*_earth_*_rpm27.npz'))
    if files:
        n = len(files)
        cols = min(5, n)
        rows_ = math.ceil(n / cols)
        fig = plt.figure(figsize=(3.6 * cols, 3.5 * rows_))
        for i, p in enumerate(files):
            d = np.load(p)
            a3 = fig.add_axes([(i % cols) / cols, 1 - (i // cols + 1) / rows_ * 0.96, 1 / cols, 0.96 / rows_],
                              projection='3d')
            wn = 'moon' if '_moon_' in p.name else 'earth'
            _mesh(a3, d['faces'], d['verts'], MOON_C if wn == 'moon' else EARTH_C, lim=0.62)
            a3.view_init(elev=14, azim=-70)
            name, cutpart, hempart = p.stem.replace('turn_', '').split(f'_{wn}_')[0], *p.stem.split(f'_{wn}_')[1].split('_')[:2]
            cut = 'half circle' if cutpart == 'cut180' else 'full circle'
            hem = int(hempart.replace('hem', ''))
            a3.set_title(f"{_short(name)}, {cut}{f', {hem} g/m hem' if hem else ''}, {wn}", fontsize=8, y=0.9)
        fig.suptitle('Skirts turning steadily at 27 rpm, seen from the side (the body is not drawn)', fontsize=9)
        fig.savefig(out / 'turn_shapes.png', dpi=110)
        plt.close(fig)


def stream_figure(G, R, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rs = R.get('stream', [])
    if not rs:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    est = {r['fabric']: r for r in G['fabrics']}
    for fab in sorted({r['fabric'] for r in rs}):
        for wn, col in (('moon', MOON_C), ('earth', EARTH_C)):
            pts = sorted([r for r in rs if r['fabric'] == fab and r['world'] == wn], key=lambda r: r['U_m_s'])
            if not pts:
                continue
            U = [r['U_m_s'] for r in pts]
            ax.errorbar(U, [r['trail_deg_mean'] for r in pts], yerr=[r['trail_deg_std'] for r in pts], fmt='o-',
                        color=col, ms=4, lw=1, capsize=2, alpha=0.9)
            ax.plot(U, [est[fab][wn]['trail_deg'][f'{u:g}_m_s'] for u in U], ':', color=col, lw=1)
            if wn == 'moon':
                ax.annotate(_short(fab), (U[-1], pts[-1]['trail_deg_mean']), fontsize=6.5, xytext=(4, -2),
                            textcoords='offset points', color=col)
            else:
                ax.annotate(_short(fab), (U[0], pts[0]['trail_deg_mean']), fontsize=6.5, xytext=(-6, 0),
                            textcoords='offset points', color=col, ha='right')
    ax.plot([], [], 'o-', color=MOON_C, label='moon')
    ax.plot([], [], 'o-', color=EARTH_C, label='earth')
    ax.legend(fontsize=7.5, frameon=False)
    ax.set_xlim(0.55, 3.45)
    ax.set_xlabel('travel speed (m/s)')
    ax.set_ylabel('panel trailing below horizontal (deg)')
    ax.set_title('A 1.2 m panel trailing in travel (solver; dotted: straight-panel estimate)', fontsize=9)
    ax.grid(lw=0.3, alpha=0.5)
    fig.tight_layout()
    fig.savefig(out / 'stream.png', dpi=130)
    plt.close(fig)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    d = Path(argv[0])
    out = d / 'figures'
    out.mkdir(parents=True, exist_ok=True)
    G = json.loads(SCREEN.read_text())
    R = json.loads((d / 'trials.json').read_text()) if (d / 'trials.json').exists() else {}
    screen_figure(G, out)
    drape_figure(R, out, d / 'shapes', {r['fabric']: r.get('bending_status') for r in G['fabrics']})
    turn_figure(G, R, out, d / 'shapes')
    stream_figure(G, R, out)
    print('figures in', out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
