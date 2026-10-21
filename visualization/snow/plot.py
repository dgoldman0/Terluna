"""Plot versioned snow products without importing or rerunning their physics.

    OPENBLAS_NUM_THREADS=1 python visualization/snow/plot.py

Generated figures and the manifest retain hashes and exact plotted values. The
cumulative curves normalize stored bin moments; no within-bin size distribution
is invented. Direct normal and diffuse horizontal light are separate series.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, FuncFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / 'generated'
PRODUCTS = {
    'population': ('atmosphere/snow/results/population.json', 'terluna.atmosphere.snow-population/1'),
    'turbulence': ('atmosphere/snow/results/turbulence.json', 'terluna.atmosphere.snow-turbulence/1'),
    'optics': ('atmosphere/snow/results/optics.json', 'terluna.atmosphere.snow-optics/1'),
}
GEOMETRY = 'bf95_dmax_dendrites_area'
BOXES = (
    ('box_0e_small', 'Reference box'),
    ('box_0e_elec_corrected', 'Electrified box'),
)
GROUPS = (
    ('all_winter_air', 'All cold air'),
    ('snow_ge_1e-05_kg_kg', 'Snow-bearing'),
    ('snow_and_liquid_cloud', 'Snow + droplets'),
)
COLORS = ('#2367aa', '#df7224', '#23866e')
INK, MUTED, GRID = '#172b3a', '#52646d', '#dce3e6'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_products():
    products, inputs = {}, {}
    for name, (relative, schema) in PRODUCTS.items():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f'Required {name} product is missing: {relative}. Generate the domain product first.')
        raw = path.read_bytes()
        product = json.loads(raw)
        if product.get('schema') != schema:
            raise ValueError(f'{relative}: expected {schema}, received {product.get("schema")}')
        products[name] = product
        inputs[relative] = hashlib.sha256(raw).hexdigest()
    pop_path = PRODUCTS['population'][0]
    if products['optics'].get('inputs', {}).get(pop_path) != inputs[pop_path]:
        raise ValueError('Optics does not bind the population being plotted; regenerate the optics product.')
    return products, inputs


def one(rows, predicate, description):
    selected = [row for row in rows if predicate(row)]
    if len(selected) != 1:
        raise ValueError(f'Expected one {description}; found {len(selected)}')
    return selected[0]


def vector(value, name, nonnegative=True):
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or not np.isfinite(array).all():
        raise ValueError(f'{name}: finite one-dimensional data required')
    if nonnegative and np.any(array < 0):
        raise ValueError(f'{name}: negative values')
    return array


def check_total(actual, expected, name):
    if not np.isclose(actual, expected, rtol=1e-9, atol=1e-18):
        raise ValueError(f'{name}: stored total {expected} differs from plotted sum {actual}')


def cumulative(distribution, weight, total_key):
    raw_edges = distribution['size_edges_m']
    open_tail = distribution.get('open_ended_last_bin', False)
    if open_tail:
        if not raw_edges or raw_edges[-1] is not None:
            raise ValueError('Open-ended distribution must have one null final edge')
        edges = vector(raw_edges[:-1], 'finite diameter edges')
    else:
        edges = vector(raw_edges, 'diameter edges')
    weights = vector(weight, total_key)
    if len(raw_edges) != len(weights) + 1 or np.any(np.diff(edges) <= 0) or edges[0] <= 0:
        raise ValueError('Increasing positive bin edges must bracket every weight')
    total = float(weights.sum())
    check_total(total, distribution['totals'][total_key], total_key)
    if total <= 0:
        raise ValueError(f'Cannot normalize an empty {total_key} distribution')
    finite_weights = weights[:-1] if open_tail else weights
    overflow = float(weights[-1] / total) if open_tail else 0.
    fractions = np.r_[0., np.cumsum(finite_weights) / total]
    if abs(fractions[-1] + overflow - 1.) > 1e-12 or np.any(np.diff(fractions) < -1e-14):
        raise ValueError('Cumulative finite bins and overflow do not sum to one')
    occupied = np.flatnonzero(finite_weights > 0)
    low = edges[occupied[0]] if occupied.size else edges[-2]
    high = edges[-1] if overflow > 0 or not occupied.size else edges[occupied[-1]+1]
    return dict(diameter_edges_cm=(edges * 100).tolist(), cumulative_fraction=fractions.tolist(),
                bin_weights=weights.tolist(), total=total,
                open_ended_last_bin=bool(open_tail), overflow_fraction=overflow,
                overflow_lower_dimension_cm=float(edges[-1]*100) if open_tail else None,
                occupied_diameter_cm=[float(low * 100), float(high * 100)])


def figure_data(products):
    pop, turb, optics = (products[key] for key in ('population', 'turbulence', 'optics'))
    curves, selected = [], {}
    for rate in (.3, 1., 3.):
        identity = f'moon_{rate:g}_100_{GEOMETRY}_s0'
        run = one(pop['runs'], lambda row: row['id'] == identity, identity)
        if (run['world'] != 'moon' or run['seed'] != 0 or run['geometry'] != GEOMETRY
                or run['crystals_per_litre'] != 100. or run['snowfall_set_mm_h'] != rate):
            raise ValueError(f'{identity}: case identity and parameters disagree')
        selected[rate] = run
        level = run['freezing_level']
        curves.append(dict(run_id=identity, imposed_condensation_mm_h=rate,
                           arriving_water_mm_h=run['summary']['snowfall_mm_h'],
                           **cumulative(level, level['mass_flux_kg_m2_s'], 'mass_flux_kg_m2_s')))
    heavy = selected[3.]
    level = heavy['freezing_level']
    number = vector(level['number_m3'], 'number concentration')
    area = vector(level['projected_area_m2'], 'mean projected area')
    if number.shape != area.shape:
        raise ValueError('Number and projected-area bins differ')
    area_weight = number * area
    optical = one(optics['distributions'], lambda row: row['run_id'] == heavy['id'] and row['level'] == 'freezing_level',
                  'heavy-case freezing-level optics')
    # The optical product specifies one size-independent Q_ext. Its extinction
    # fractions therefore equal the already computed projected-area fractions.
    # This is normalization of stored moments, not a new scattering model.
    efficiency = float(optics['assumptions']['q_ext'])
    check_total(float(area_weight.sum()) * efficiency, optical['extinction_m_inverse'], 'extinction total')
    weights = dict(
        number=cumulative(level, number, 'number_m3'),
        ice_mass=cumulative(level, level['mass_kg_m3'], 'mass_kg_m3'),
        extinction=cumulative(level, area_weight, 'projected_area_m_inv'),
    )
    bars = []
    for case_id, label in BOXES:
        case = one(turb['cases'], lambda row: row['case'] == case_id, case_id)
        records = case['sampling']['records']
        first, last = records[0]['model_day'], records[-1]['model_day']
        if last - first < 29.:
            raise ValueError(f'{case_id}: full-cycle comparison requires at least 29 days of sampled time')
        conditions = case['conditioning']
        if conditions['temperature_c'] != [-20., 0.] or conditions['height_m'] != [20000., 60000.]:
            raise ValueError(f'{case_id}: unexpected winter temperature or height selection')
        values = []
        for key, group_label in GROUPS:
            row = case['summaries'][key]
            fraction = row['fractions']['abs_w_gt_0.6_m_s']
            if fraction is None or not 0 <= fraction <= 1 or row['sample_count'] < 1:
                raise ValueError(f'{case_id}/{key}: missing conditioned wind samples')
            values.append(dict(group=key, label=group_label, fraction=fraction, sample_count=row['sample_count']))
        bars.append(dict(case=case_id, label=label, model_day_range=[first, last],
                         snapshots=case['sampling']['snapshots'], conditioning=conditions, values=values))
    light = sorted((row for row in optics['lighting']['rows']
                    if row['source'] == 'shielded_sun' and row['altitude_m'] == 35000.),
                   key=lambda row: row['source_elevation_deg'])
    if len(light) < 3:
        raise ValueError('Need at least three shielded-solar elevations at 35 km')
    angles = vector([row['source_elevation_deg'] for row in light], 'solar elevations', nonnegative=False)
    if np.any(np.diff(angles) <= 0):
        raise ValueError('Repeated or unordered source elevations')
    direct = vector([row['direct_normal_xyz'][1] for row in light], 'direct normal lux')
    diffuse = vector([row['diffuse_down_horizontal_xyz'][1] for row in light], 'diffuse downward horizontal lux')
    return dict(arriving_mass_curves=curves, heavy_run_id=heavy['id'], concentration_curves=weights,
                extinction_efficiency=efficiency, optical_total_m_inverse=optical['extinction_m_inverse'],
                wind_bars=bars,
                incident_light=dict(altitude_m=35000., source_elevation_deg=angles.tolist(),
                                    direct_normal_lux=direct.tolist(), diffuse_down_horizontal_lux=diffuse.tolist(),
                                    zero_direct_elevations_deg=angles[direct == 0].tolist()),
                interpretation=dict(cdf='Cumulative normalized finite-bin moments at upper bin edges; open-ended overflow is retained separately, with no invented finite upper edge.',
                    population='Moon, seed 0, corrected maximum dimension and representative dendrite area; 100/L flux normalization.',
                    numerical_limitation='Raw population outputs inherit inconsistent fractional-multiplicity rounding and fragment cutoffs, plus batched sticking variance. Extreme unbounded tails are not physical flake-size predictions. See atmosphere/snow/numerical_audit.md.',
                    weighting='Panel a uses arriving mass flux. Panel b uses instantaneous concentration, not arriving flux.',
                    resolved_wind='Equal conditioned cell-times in two 29.5-day archived box intervals; not flake-scale turbulence.',
                    light='Separate receiving-plane conventions: normal direct illuminance and downward horizontal diffuse illuminance.',
                    zero_light='Exact zero direct values are retained in the manifest and omitted from the logarithmic line.'))


def style_axis(ax, title, subtitle):
    ax.set_title(title, loc='left', fontsize=13, weight='bold', color=INK, pad=37)
    ax.text(0, 1.035, subtitle, transform=ax.transAxes, color=MUTED, fontsize=9, va='bottom', linespacing=1.4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_color('#aab8bf')
    ax.tick_params(labelsize=9, colors=MUTED)
    ax.grid(axis='y', color=GRID, linewidth=.65)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)


def cdf_axis(ax, curves):
    ax.set_xscale('log')
    ax.set_ylim(-.015, 1.03)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.set_xlabel('Maximum aggregate dimension (cm)', fontsize=10)
    ax.set_ylabel('Cumulative share', fontsize=10)
    low = min(curve['occupied_diameter_cm'][0] for curve in curves)
    high = max(curve['occupied_diameter_cm'][1] for curve in curves)
    ax.set_xlim(low / 1.2, high * 1.2)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f'{value:g}'))
    ax.legend(loc='lower right', frameon=False, fontsize=9)


def tail_labels(ax, curves):
    if not any(curve['open_ended_last_bin'] for _, curve in curves):
        return
    lower = next(curve['overflow_lower_dimension_cm']/100 for _, curve in curves if curve['open_ended_last_bin'])
    lines = [f'Raw model tail ≥ {lower:g} m (not shown); no imposed size cap',
             ' · '.join(f'{name}: {100*curve["overflow_fraction"]:.3g}%' for name, curve in curves)]
    ax.text(0, -.24, '\n'.join(lines), transform=ax.transAxes, va='top', ha='left', fontsize=8.1,
            color=MUTED, linespacing=1.4)


def render(data, provenance, output):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.hashsalt': 'terluna-snow-overview',
                         'axes.labelsize': 10, 'savefig.facecolor': '#fcfcfb', 'figure.facecolor': '#fcfcfb'})
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.subplots_adjust(left=.08, right=.97, top=.84, bottom=.18, hspace=.72, wspace=.28)
    fig.suptitle('Aerial winter: snow, motion and light', x=.08, y=.973, ha='left', fontsize=22, color=INK, weight='bold')
    fig.text(.08, .926, 'Conditional population calculations, archived weather and an independently solved clear sky',
             color=MUTED, fontsize=11)
    ax = axes[0, 0]
    style_axis(ax, 'a  Raw model sizes by supplied snowfall',
               'Arriving mass; numerical size sensitivity unresolved\n100/L input normalization · corrected-dimension, dendrite-area scenario')
    curves = data['arriving_mass_curves']
    for color, curve in zip(COLORS, curves):
        ax.step(curve['diameter_edges_cm'], curve['cumulative_fraction'], where='post', color=color, lw=2,
                label=f'{curve["imposed_condensation_mm_h"]:g} mm/h')
    cdf_axis(ax, curves)
    tail_labels(ax, [(f'{curve["imposed_condensation_mm_h"]:g} mm/h', curve) for curve in curves])
    ax = axes[0, 1]
    style_axis(ax, 'b  Raw model count, mass and optical area',
               'Heavy case: 3 mm/h imposed condensation\nInstantaneous concentrations at the freezing level; grey Qext = 2')
    curves = data['concentration_curves']
    for (key, label), color in zip((('number', 'Number'), ('ice_mass', 'Ice mass'), ('extinction', 'Extinction')), COLORS):
        curve = curves[key]
        ax.step(curve['diameter_edges_cm'], curve['cumulative_fraction'], where='post', color=color, lw=2, label=label)
    cdf_axis(ax, list(curves.values()))
    tail_labels(ax, [('Number', curves['number']), ('Ice mass', curves['ice_mass']), ('Extinction', curves['extinction'])])
    ax = axes[1, 0]
    style_axis(ax, 'c  Snow-bearing air is less often quiet',
               'Share of selected samples with |vertical wind| > 0.6 m/s\n−20 to 0 °C, 20–60 km; equal cell-time weights')
    x, width = np.arange(len(GROUPS)), .34
    for i, case in enumerate(data['wind_bars']):
        heights = [row['fraction'] for row in case['values']]
        bars = ax.bar(x + (i-.5)*width, heights, width, color=COLORS[i], label=case['label'])
        for bar, value in zip(bars, heights):
            ax.text(bar.get_x()+bar.get_width()/2, value+.018, f'{100*value:.1f}%', ha='center', va='bottom',
                    fontsize=9, color=INK)
    ax.set_xticks(x, [label for _, label in GROUPS])
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.set_ylim(0, 1.10)
    ax.set_ylabel('Within each conditioned sample', fontsize=10)
    ax.legend(loc='upper left', frameon=False, fontsize=9)
    ax.text(0, -.22, 'Snow and droplet thresholds: each ≥ 0.01 g/kg.\nResolved bulk motion; centimetre-scale turbulence is unresolved.',
            transform=ax.transAxes, color=MUTED, fontsize=8.5, va='top', linespacing=1.4)
    ax = axes[1, 1]
    style_axis(ax, 'd  Incident illumination at 35 km',
               'Shielded solar source; molecular clear sky before snow\nDirect normal and diffuse horizontal values are separate')
    light = data['incident_light']
    angle = np.asarray(light['source_elevation_deg'])
    for key, label, color in [('direct_normal_lux', 'Direct normal', COLORS[1]),
                               ('diffuse_down_horizontal_lux', 'Diffuse downward horizontal', COLORS[0])]:
        values = np.asarray(light[key])
        ax.plot(angle, np.where(values > 0, values, np.nan), '-o', ms=4, lw=2, color=color, label=label)
    ax.set_yscale('log')
    ax.set_xlim(float(angle.min())-2, float(angle.max())+2)
    ax.set_xlabel('Geometric solar elevation (degrees)', fontsize=10)
    ax.set_ylabel('Illuminance (lux)', fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc='lower right')
    zeros = ', '.join(f'{value:g}°' for value in light['zero_direct_elevations_deg'])
    ax.text(0, -.22, f'Exact direct-beam zeros omitted on the log axis: {zeros or "none"}.\nIncident field is not the light inside an optically thick storm.',
            transform=ax.transAxes, color=MUTED, fontsize=8.5, va='top', linespacing=1.4)
    fig.text(.08, .036, 'Uncapped population outputs contain extreme tails; unresolved numerical sensitivity prevents interpreting these as physical flake-size predictions.\n'
             'Persistent-cloud populations, weather boxes and clear-sky illumination are distinct calculations, not a coupled snowstorm forecast.',
             color=MUTED, fontsize=9.3, linespacing=1.45)
    encoded = json.dumps(provenance, sort_keys=True, separators=(',', ':'))
    png, svg = output/'overview.png', output/'overview.svg'
    fig.savefig(png, dpi=180, metadata={'Title': 'Aerial winter: snow, motion and light',
                                      'Description': encoded, 'Software': 'Terluna visualization/snow/plot.py'})
    fig.savefig(svg, metadata={'Title': 'Aerial winter: snow, motion and light', 'Description': encoded,
                              'Creator': 'Terluna visualization/snow/plot.py', 'Date': None})
    plt.close(fig)
    return [png, svg]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    products, inputs = read_products()
    data = figure_data(products)
    source = {'visualization/snow/plot.py': digest(__file__)}
    provenance = dict(schema='terluna.visualization.snow-overview/1', inputs=inputs, source=source)
    args.output.mkdir(parents=True, exist_ok=True)
    outputs = render(data, provenance, args.output)
    manifest = dict(provenance, plotted_values=data,
                    evidence={key: product.get('evidence') for key, product in products.items()},
                    output_sha256={path.name: digest(path) for path in outputs},
                    checks=['Expected product schemas', 'Optics binds the plotted population hash',
                            'CDF moments reproduce stored totals', 'Optical area total reproduces extinction',
                            'Finite CDFs monotonic; finite and overflow fractions sum to one', 'Both weather boxes span at least 29 days'])
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {args.output / "overview.png"}')
    print(f'Wrote {args.output / "overview.svg"}')
    print(f'Wrote {args.output / "manifest.json"}')


if __name__ == '__main__':
    main()
