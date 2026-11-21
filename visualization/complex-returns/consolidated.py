"""Plot reviewed sponsor cash returns without importing the economic model."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'research/studies/complex_returns'
SOURCE = STUDY / 'results/consolidated.json'
OUTPUT = ROOT / 'research/runs/complex_returns/consolidated-returns'
SCHEMA = 'terluna.research.consolidated-returns/1'
METRIC = 'sponsor_cash_NPV_with_saving_release'
HORIZONS = (30, 100, 300, 800)
PORTFOLIOS = (
    ('raw_resource_price', 'Raw-resource price proxy'),
    ('commodity_refined', 'Commodity refined'),
    ('specialty', 'Specialty metals'),
    ('premium_lunar', 'Premium lunar products'),
    ('propellant_logistics', 'Propellant and logistics'),
    ('compute_digital', 'Computing and digital'),
    ('services', 'Worldwide services'),
    ('diversified', 'Diversified'),
    ('prior_manufacturing_mix', 'Earlier manufacturing mix'),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_product(path):
    """Require the current producer and every declared dependency to match."""
    product = json.loads(path.read_text())
    if product['schema'] != SCHEMA:
        raise ValueError(f"Unexpected consolidated-return schema: {product['schema']}")
    producer = STUDY / 'consolidated.py'
    if digest(producer) != product['producer_sha256']:
        raise ValueError('Consolidated-return producer hash mismatch')
    for name, expected in product['dependencies'].items():
        if digest(STUDY / name) != expected:
            raise ValueError(f'Consolidated-return dependency hash mismatch: {name}')
    return product


def display_data(product):
    """Select saved values and convert their display units, with no new return calculation."""
    rows = {}
    for row in product['financial']:
        if row['mode'] == 'latest_U_package':
            key = (row['strategy'], row['market'], row['horizon'])
            if key in rows:
                raise ValueError(f'Duplicate financial result: {key}')
            rows[key] = row
    displayed = []

    def value(name, market, horizon, rate):
        row = rows[name, market, horizon]
        raw = row['discounted'][rate][METRIC]
        if not math.isfinite(raw):
            raise ValueError('Nonfinite sponsor cash NPV')
        continuing = row['status'] == 'continuing_conditional_path'
        displayed.append(dict(strategy=name, market=market, horizon=horizon,
                              real_discount_rate=rate, status=row['status'],
                              saved_cash_NPV_USD2041=raw,
                              plotted_cash_NPV_trillion_USD2041=raw / 1e12 if continuing else None))
        return raw / 1e12 if continuing else np.nan

    heat = np.array([[value(name, 'breakthrough', year, '0.03')
                      for year in HORIZONS] for name, _ in PORTFOLIOS])
    middle = np.array([value(name, 'global', 100, '0.03') for name, _ in PORTFOLIOS])
    strong = np.array([value(name, 'breakthrough', 100, '0.03') for name, _ in PORTFOLIOS])
    higher_rate = np.array([value(name, 'breakthrough', 100, '0.07') for name, _ in PORTFOLIOS])
    return heat, middle, strong, higher_rate, displayed


def render(source=SOURCE, output=OUTPUT):
    product = read_product(source)
    settings = json.loads((STUDY / 'consolidated_inputs.json').read_text())
    heat, middle, strong, higher_rate, displayed = display_data(product)
    identity = dict(product_sha256=digest(source), producer_sha256=product['producer_sha256'],
                    renderer_sha256=digest(Path(__file__)), dependencies=product['dependencies'])
    metadata = dict(schema='terluna.visualization.consolidated-returns/1',
                    source=str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else str(source),
                    metric=METRIC, mode='latest_U_package', units='trillion constant USD2041',
                    evidence=product['evidence'], reading_rule=product['reading_rule'],
                    **identity)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'terluna-consolidated-returns-v1',
                         'svg.fonttype': 'none', 'axes.unicode_minus': True})
    teal, orange, navy = '#16766e', '#b15b3b', '#253d54'
    fig = plt.figure(figsize=(14, 8.2), facecolor='white')
    a = fig.add_axes((.225, .245, .265, .54))
    b = fig.add_axes((.565, .245, .415, .54))
    y = np.arange(len(PORTFOLIOS))
    labels = [label for _, label in PORTFOLIOS]
    maximum = max(1., math.ceil(float(np.nanmax(np.abs(heat))) * 2) / 2)
    cmap = LinearSegmentedColormap.from_list('return', [orange, '#fbfaf7', teal])
    cmap.set_bad('#e4e7e8')
    picture = a.imshow(heat, cmap=cmap, norm=Normalize(-maximum, maximum), aspect='auto')
    a.set(xticks=range(len(HORIZONS)), xticklabels=HORIZONS, yticks=y, yticklabels=labels,
          xlabel='Elapsed years from first investment')
    a.tick_params(axis='both', length=0, pad=9)
    a.set_xticks(np.arange(-.5, len(HORIZONS), 1), minor=True)
    a.set_yticks(np.arange(-.5, len(PORTFOLIOS), 1), minor=True)
    a.grid(which='minor', color='white', linewidth=1.5)
    a.tick_params(which='minor', bottom=False, left=False)
    for row in range(heat.shape[0]):
        for column in range(heat.shape[1]):
            v = heat[row, column]
            color = 'white' if math.isfinite(v) and abs(v) > maximum * .60 else '#172a36'
            number = f'{v:+.3f}' if math.isfinite(v) and 0 < abs(v) < .005 else f'{v:+.2f}'
            a.text(column, row, number if math.isfinite(v) else 'Failed',
                   ha='center', va='center', fontsize=10, color=color)
    for spine in a.spines.values():
        spine.set_visible(False)
    a.set_title('A  Strong demand across horizons', loc='left', fontsize=12, fontweight='bold', pad=35)
    a.text(0., 1.025, '3% real discount rate', transform=a.transAxes, fontsize=10, color='#4a5862')

    b.barh(y-.20, middle, height=.31, color=orange, edgecolor='none', zorder=3)
    b.barh(y+.17, strong, height=.31, color=teal, edgecolor='none', zorder=3)
    b.scatter(higher_rate, y+.17, color=navy, marker='D', s=23, zorder=4)
    b.set(ylim=(len(PORTFOLIOS)-.5, -.5), yticks=y, yticklabels=[],
          xlabel='Cash NPV, trillion constant USD2041')
    b.tick_params(axis='y', length=0)
    b.spines['left'].set_visible(False)
    b.spines['bottom'].set_color('#aeb7bd')
    b.axvline(0., color='#526573', linewidth=.9, zorder=2)
    b.grid(axis='x', color='#e1e5e7', linewidth=.7, zorder=0)
    b.set_xlim(min(-3., float(np.nanmin(middle))-.35), max(5., float(np.nanmax(strong))+.55))
    for i, (left, right) in enumerate(zip(middle, strong)):
        if math.isfinite(left):
            b.text(left-.10, i-.20, f'{left:+.2f}', ha='right', va='center', fontsize=9, color='#643c2f')
        if math.isfinite(right):
            b.text(right+.10, i+.17, f'{right:+.2f}', ha='left', va='center', fontsize=9, color='#14564f')
    b.set_title('B  100-year demand and discount sensitivity', loc='left', fontsize=12,
                fontweight='bold', pad=35)
    handles = [Patch(facecolor=orange, label='Middle demand · 3%'),
               Patch(facecolor=teal, label='Strong demand · 3%'),
               Line2D([], [], color=navy, marker='D', linestyle='none', markersize=4,
                      label='Strong demand · 7%')]
    b.legend(handles=handles, loc='lower left', bbox_to_anchor=(-.015, 1.002),
             ncol=3, frameon=False, fontsize=8.7, columnspacing=1.2, handlelength=1.1)
    for series, offset in ((middle, -.20), (strong, .17), (higher_rate, .17)):
        for i, v in enumerate(series):
            if not math.isfinite(v):
                b.text(0., i+offset, 'Failed', ha='center', va='center', fontsize=8)

    color_axes = fig.add_axes((.225, .155, .265, .018))
    bar = fig.colorbar(picture, cax=color_axes, orientation='horizontal')
    bar.outline.set_visible(False)
    bar.set_label('Cash NPV, trillion constant USD2041', fontsize=9.5, labelpad=4)
    bar.ax.tick_params(labelsize=9, length=2)
    fig.text(.035, .945, 'Sponsor cash returns for nine alternative export portfolios',
             fontsize=19, fontweight='bold', color='#172a36')
    fig.text(.035, .910,
             f"U package: {settings['latest_capital_supplier_U_share']:.0%} outside capital-supplier acceptance; "
             f"{settings['local_U_rate_reduction']*10000:g} basis-point local loan-rate relief",
             fontsize=11, color='#4a5862')
    fig.text(.035, .879,
             'Initial seed and monetary reserves deducted; terminal book claims excluded.',
             fontsize=10.5, color='#4a5862')
    fig.text(.035, .075,
             'Annual interest savings are released to outside saving assets; the bank loses the same income. '
             'These are sponsor returns, not global welfare.', fontsize=9.3, color='#33434e')
    fig.text(.035, .049,
             'Conditional scenarios, not forecasts. Each alternative includes worldwide service receipts; '
             'small-module engineering and future demand remain assumptions.', fontsize=9.3, color='#33434e')
    fig.text(.035, .019,
             f"SHA-256  product {identity['product_sha256'][:16]}  ·  producer {identity['producer_sha256'][:16]}"
             f"  ·  renderer {identity['renderer_sha256'][:16]}  |  Full identities and plotted values accompany the figure.",
             fontsize=7.6, color='#61717a')
    output.parent.mkdir(parents=True, exist_ok=True)
    description = json.dumps(metadata, sort_keys=True)
    png, svg = output.with_suffix('.png'), output.with_suffix('.svg')
    fig.savefig(png, dpi=200, metadata={'Description': description})
    fig.savefig(svg, metadata={'Date': None, 'Description': description})
    plt.close(fig)
    manifest = dict(metadata, displayed_values=displayed,
                    artifacts={png.name: digest(png), svg.name: digest(svg)})
    output.with_suffix('.manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    return png, svg


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--product', type=Path, default=SOURCE)
    parser.add_argument('--output-stem', type=Path, default=OUTPUT)
    args = parser.parse_args()
    for artifact in render(args.product.resolve(), args.output_stem.resolve()):
        print(artifact)
