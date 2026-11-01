"""Render frozen condensation-growth controls without importing domain physics.

    OPENBLAS_NUM_THREADS=1 python visualization/snow/growth.py

The mass CDF is a discrete sum at stored particle masses; no within-bin shape is
invented. Generated images and the manifest belong in ignored generated/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, PercentFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PRODUCTS = {
    'first_order': ('atmosphere/snow/results/growth_transport.json',
                    'terluna.atmosphere.snow-growth-transport/1'),
    'limited': ('atmosphere/snow/results/limited_growth.json',
                'terluna.atmosphere.snow-limited-growth/1'),
}
GRIDS = (96, 144, 216)
COLORS = {'first_order': '#bf5b1b', 'limited': '#2367aa', 'characteristic': '#172b3a'}
INK, MUTED, GRID = '#172b3a', '#52646d', '#dce3e6'
VIEW_MG = (.02, 8.)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_products(root=ROOT):
    products, identities = {}, {}
    for key, (relative, schema) in PRODUCTS.items():
        raw = (root / relative).read_bytes()
        product = json.loads(raw)
        if product.get('schema') != schema:
            raise ValueError(f'{relative}: unsupported schema')
        if not product.get('evidence') or not product.get('reading_rule'):
            raise ValueError(f'{relative}: missing evidence or reading rule')
        products[key] = product
        identities[relative] = hashlib.sha256(raw).hexdigest()
    reference = PRODUCTS['first_order'][0]
    if products['limited']['producer']['inputs'].get(reference) != identities[reference]:
        raise ValueError('Limited-growth product does not bind the first-order product being plotted')
    verified = {}
    for product in products.values():
        for group in ('files', 'inputs'):
            for relative, expected in product['producer'][group].items():
                actual = verified.setdefault(relative, digest(root / relative))
                if actual != expected:
                    raise ValueError(f'Upstream {group} identity mismatch: {relative}')
    if products['first_order']['settings'] != products['limited']['settings']:
        raise ValueError('Products use different physical settings')
    return products, identities, verified


def one(rows, identity):
    matches = [row for row in rows if row['id'] == identity]
    if len(matches) != 1:
        raise ValueError(f'Expected one {identity}, found {len(matches)}')
    return matches[0]


def check(actual, expected, description):
    if not np.isclose(actual, expected, rtol=2e-10, atol=0.):
        raise ValueError(f'{description}: reconstructed {actual} differs from stored {expected}')


def curve(masses, numbers, snapshot, identity):
    mass, number = np.asarray(masses, float), np.asarray(numbers, float)
    if (mass.ndim != 1 or number.shape != mass.shape or not np.isfinite(mass).all()
            or not np.isfinite(number).all() or np.any(mass <= 0) or np.any(number < 0)):
        raise ValueError(f'{identity}: invalid particle masses or numbers')
    total_number = number.sum()
    total_mass = np.dot(number, mass)
    if total_number <= 0 or total_mass <= 0:
        raise ValueError(f'{identity}: empty distribution')
    mean = total_mass / total_number
    variance = np.dot(number, (mass - mean)**2) / total_number
    moments = dict(number_m3=float(total_number), mass_kg_m3=float(total_mass),
                   mean_mass_kg=float(mean), variance_mass_kg2=float(variance),
                   relative_mass_variance=float(variance / mean**2))
    for key, value in moments.items():
        check(value, snapshot[key], f'{identity}/{key}')
    order = np.argsort(mass, kind='stable')
    mass, number = mass[order], number[order]
    weights = mass * number
    cdf = np.cumsum(weights) / total_mass
    check(cdf[-1], 1., f'{identity}/CDF normalization')
    if np.any(np.diff(cdf) < 0):
        raise ValueError(f'{identity}: nonmonotone CDF')
    in_mg = mass * 1e6
    return dict(id=identity, particle_mass_mg=in_mg.tolist(), number_m3=number.tolist(),
                cumulative_mass_fraction=cdf.tolist(), moments=moments,
                mass_fraction_below_view=float(weights[in_mg < VIEW_MG[0]].sum()/total_mass),
                mass_fraction_above_view=float(weights[in_mg > VIEW_MG[1]].sum()/total_mass))


def figure_data(products):
    first, limited = products['first_order'], products['limited']
    duration = first['settings']['duration_s']
    physical = one(first['physical_characteristics'], 'q128_rtol1e-11')
    physical_curve = curve(physical['final_masses_kg'], physical['number_weights_m3'],
                           physical['snapshots'][-1], physical['id'])
    check(physical['snapshots'][-1]['time_s'], duration, 'physical endpoint time')
    panels = []
    grid_bounds = None
    for bins in GRIDS:
        identity = f'm{bins}_dt30'
        a, b = one(first['runs'], identity), one(limited['runs'], identity)
        if (a['mass_bins'] != bins or b['mass_bins'] != bins
                or a['max_step_s'] != 30. or b['max_step_s'] != 30.
                or a['mass_pivots_kg'] != b['mass_pivots_kg']):
            raise ValueError(f'{identity}: inconsistent grid or time-step limit')
        pivots = np.asarray(a['mass_pivots_kg'])
        if len(pivots) != bins or np.any(np.diff(pivots) <= 0):
            raise ValueError(f'{identity}: invalid pivot grid')
        bounds = [float(pivots[0]), float(pivots[-1])]
        if grid_bounds is None:
            grid_bounds = bounds
        elif bounds != grid_bounds:
            raise ValueError('Grid refinement changes the represented mass domain')
        ref = one(first['grid_matched_characteristics'], f'projected_seed_m{bins}')
        final_ref = ref['snapshots'][-1]
        check(final_ref['time_s'], duration, f'{identity}/reference endpoint time')
        curves = {'characteristic': curve(ref['final_masses_kg'], ref['number_weights_m3'],
                                          final_ref, ref['id'])}
        ledgers = {}
        for key, run in (('first_order', a), ('limited', b)):
            final = run['snapshots'][-1]
            check(final['time_s'], duration, f'{identity}/{key} endpoint time')
            for moment in ('number_m3', 'mass_kg_m3', 'relative_mass_variance'):
                check(final['characteristic'][moment], final_ref[moment],
                      f'{identity}/{key} characteristic {moment}')
            curves[key] = curve(pivots, run['final_numbers_m3'], final, f'{key}/{identity}')
            ledgers[key] = final['water_ledger']
            if final['water_ledger']['upper_exit'] / final['mass_kg_m3'] > 1e-8:
                raise ValueError(f'{identity}/{key}: upper exit needs an explicit figure treatment')
        panels.append(dict(mass_bins=bins, run_id=identity, curves=curves,
                           water_ledgers=ledgers,
                           once_projected_reference_relative_variance=
                           a['snapshots'][-1]['once_projected_characteristic']['relative_mass_variance']))
    support = [min(physical_curve['particle_mass_mg']), max(physical_curve['particle_mass_mg'])]
    max_hidden = max(c['mass_fraction_below_view'] + c['mass_fraction_above_view']
                     for p in panels for c in p['curves'].values())
    if max_hidden >= 1e-4:
        raise ValueError('Common display window excludes at least 0.01% of a curve; revise the view')
    return dict(settings=first['settings'], panels=panels, physical_reference=physical_curve,
                physical_support_mg=support, represented_mass_domain_kg=grid_bounds,
                displayed_mass_window_mg=list(VIEW_MG), maximum_mass_fraction_outside_view=max_hidden,
                interpretation={
                    'cdf': 'CDF of ice mass remaining in the box: sum(m*n) up to each stored particle mass, divided by sum(m*n) over the entire array. Curves are discrete steps, including the quadrature reference; no within-bin shape is assumed.',
                    'variance': 'Number-weighted variance of particle mass divided by squared number-weighted mean mass; dimensionless. All values are reconstructed from the plotted arrays and checked against stored moments.',
                    'reference': 'Each black curve independently grows the same projected seed population used by its grid; characteristic masses remain continuous and are not projected again for this plot.',
                    'support': 'The shared shaded band spans the independently grown original continuous seed support, including its zero-weight endpoint tracers. Projecting the initial seeds onto each grid slightly changes its own reference support.',
                    'view': 'Identical logarithmic mass axes in every CDF panel. All tails stay in the normalization and manifest; less than 0.01% of each curve lies outside the common display window.',
                    'duration': '132.589 h is the closed-box simulation duration selected to stress numerical growth transport. It is not a flake residence time.',
                    'scope': 'Fixed-air condensation only: no collisions, settling or new seed injection. These numerical controls do not predict a lunar snowstorm or validate ice physics empirically.'})


def render(data, output, metadata):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.labelcolor': INK, 'text.color': INK,
                         'svg.hashsalt': 'terluna-snow-growth-v1'})
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.5))
    fig.subplots_adjust(left=.09, right=.97, bottom=.17, top=.79, hspace=.43, wspace=.25)
    fig.text(.09, .95, 'Numerical spreading during vapour growth', fontsize=21, weight='bold')
    settings = data['settings']
    fig.text(.09, .911, f'Closed condensation-only box  ·  {settings["duration_s"]/3600:.3f} h simulation stress  ·  '
             f'{settings["temperature_c"]:g} °C, {settings["pressure_pa"]/1000:.2f} kPa',
             fontsize=11, color=MUTED)
    handles = [Line2D([], [], color=COLORS[k], lw=2.3,
                       ls='--' if k == 'characteristic' else '-', label=label)
               for k, label in (('first_order', 'First-order'), ('limited', 'Limited'),
                                ('characteristic', 'Independent characteristics'))]
    handles.append(Patch(facecolor='#cad4d8', alpha=.5, label='Original seed support after growth'))
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.084, .886),
               frameon=False, ncol=2, fontsize=10, columnspacing=2.2)
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#aab8bf')
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.grid(axis='y', color=GRID, linewidth=.7)
        ax.set_axisbelow(True)
    for letter, ax, panel in zip('abc', axes.flat, data['panels']):
        ax.axvspan(*data['physical_support_mg'], color='#cad4d8', alpha=.5, linewidth=0)
        for key in ('first_order', 'limited', 'characteristic'):
            row = panel['curves'][key]
            mass, cdf = np.asarray(row['particle_mass_mg']), np.asarray(row['cumulative_mass_fraction'])
            ax.step(np.r_[mass[0]/2, mass, mass[-1]*2], np.r_[0., cdf, cdf[-1]], where='post',
                    color=COLORS[key], linewidth=1.8, ls='--' if key == 'characteristic' else '-')
        ax.set_xscale('log')
        ax.set_xlim(VIEW_MG)
        ax.set_ylim(-.015, 1.025)
        ax.set_xticks([.02, .05, .1, .2, .5, 1, 2, 5])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
        ax.yaxis.set_major_formatter(PercentFormatter(1.))
        ax.set_xlabel('Particle mass (mg; logarithmic axis)')
        ax.set_ylabel('Cumulative share of ice mass')
        ax.set_title(f'{letter}  {panel["mass_bins"]} mass pivots', loc='left', pad=11, weight='bold', fontsize=12)
    ax = axes[1, 1]
    for key in ('first_order', 'limited', 'characteristic'):
        values = [p['curves'][key]['moments']['relative_mass_variance'] for p in data['panels']]
        ax.plot(GRIDS, values, color=COLORS[key], marker='o', lw=2,
                ls='--' if key == 'characteristic' else '-')
        ax.annotate(f'{values[-1]:.4f}', (GRIDS[-1], values[-1]), xytext=(7, 0),
                    textcoords='offset points', fontsize=9, color=COLORS[key], va='center')
    physical_variance = data['physical_reference']['moments']['relative_mass_variance']
    ax.axhline(physical_variance, color='#8b9ca5', linewidth=1, zorder=0)
    ax.text(.04, .18, 'Thin grey line: original continuous seeds', transform=ax.transAxes,
            color=MUTED, fontsize=8.5)
    ax.set_yscale('log')
    ax.set_ylim(.018, 1.)
    ax.set_xlim(88, 244)
    ax.set_xticks(GRIDS)
    ax.set_yticks([.02, .05, .1, .2, .5, 1.])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.set_xlabel('Number of mass pivots')
    ax.set_ylabel('Relative mass variance: variance / mean²\n(number weighted; logarithmic axis)')
    ax.set_title('d  Distribution width with grid refinement', loc='left', pad=11, weight='bold', fontsize=12)
    fig.text(.09, .106, 'Black curves grow each grid’s initial projected seeds independently. '
             'All runs use a 30 s time-step cap.', fontsize=9.3, color=MUTED)
    fig.text(.09, .083, 'CDFs include all remaining ice in their normalization; less than 0.01% lies outside the displayed mass window.',
             fontsize=9.3, color=MUTED)
    fig.text(.09, .049, 'The duration is a numerical stress test, not a flake residence time. '
             'No collisions or settling; not a lunar snowstorm forecast.', fontsize=9.3, color=INK)
    output.mkdir(parents=True, exist_ok=True)
    paths = [output / 'growth_comparison.svg', output / 'growth_comparison.png']
    fig.savefig(paths[0], metadata={'Date': None, 'Title': 'Numerical spreading during vapour growth',
                                  'Description': json.dumps(metadata, sort_keys=True)})
    fig.savefig(paths[1], dpi=180, metadata={'Description': json.dumps(metadata, sort_keys=True)})
    plt.close(fig)
    return {p.name: digest(p) for p in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'generated')
    args = parser.parse_args()
    products, inputs, verified = read_products()
    data = figure_data(products)
    identity = dict(renderer={'visualization/snow/growth.py': digest(__file__)}, inputs=inputs)
    outputs = render(data, args.output, identity)
    manifest = dict(schema='terluna.visualization.snow-growth/1', **identity,
                    verified_upstream_files_and_inputs=verified,
                    environment={'matplotlib': matplotlib.__version__, 'numpy': np.__version__},
                    upstream={k: {field: p[field] for field in ('schema', 'producer', 'evidence', 'reading_rule')}
                              for k, p in products.items()},
                    plotted=data, outputs=outputs)
    path = args.output / 'growth_manifest.json'
    path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'manifest': str(path), 'outputs': list(outputs),
                      'maximum_mass_fraction_outside_view': data['maximum_mass_fraction_outside_view']}))


if __name__ == '__main__':
    main()
