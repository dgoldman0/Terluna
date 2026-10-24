"""Render the early/later growth product and its bound predecessor."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PRODUCT = ROOT/'atmosphere/snow/results/mcsnow_growth_limits.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render():
    data = json.loads(PRODUCT.read_text())
    if data['schema'] != 'terluna.atmosphere.mcsnow-growth-limits/1':
        raise ValueError('Unsupported early/later growth product')
    if not data.get('evidence') or not data.get('reading_rule'):
        raise ValueError('Missing evidence or reading rule')
    expected = {**data['producer'], **data['inputs'],
                data['previous_product_path']: data['previous_product_sha256']}
    for name, sha in expected.items():
        if digest(ROOT/name) != sha:
            raise ValueError('Input identity mismatch: '+name)
    previous = json.loads((ROOT/data['previous_product_path']).read_text())
    selected = [f for f in data['searches'] if f['name'] in data['selected_searches']]
    fig, axes = plt.subplots(4, 3, figsize=(11, 11), sharex=True)
    plotted = []

    def draw(ax, ages, values, scale, label, color, style, identity):
        values = np.asarray(values)
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError('Invalid plotted values')
        line, = ax.plot(ages, values*scale, style, color=color, label=label)
        if not np.array_equal(line.get_xdata(), ages):
            raise ValueError('Plot ages changed')
        if not np.allclose(line.get_ydata()/scale, values, rtol=3e-15, atol=0):
            raise ValueError('Plot unit conversion failed')
        plotted.append(dict(identity=identity, ages_min=list(ages),
                            native_values=values.tolist(), scale=scale))

    labels = ['Mass (mg)', 'Width (mm)', 'Polar dimension (mm)']
    for i, (temperature, initialization) in enumerate(
            [(t, init) for t in data['training']['temperatures_C'] for init in ('sphere', 'axes')]):
        diameter = 0. if temperature == -14.4 else 40.
        eligible = [f for f in selected if f['temperature_C'] == temperature
                    and f['family'] == 'axis' and f['seed_diameter_um'] == diameter
                    and f['specification']['early_axis_multiplier_bounds'] == data['default_early_bounds']
                    and f['specification'].get('ventilation_max', 2.) == 2.]
        defaults = [f for f in eligible if not f['specification'].get('fixed')]
        if len(defaults) != 1:
            raise ValueError('Expected one default-bound fit per temperature')
        before = next(f for f in previous['temperature_searches']
                      if f['name'] == previous['selected_temperature_searches'][str(temperature)])
        curves = [(before, 'Previous early limit fixed', '#cc8c27'),
                  (defaults[0], 'Independent early and later limits', '#087f8c')]
        constrained = [f for f in eligible if set(f['specification'].get('fixed',{})) == {'0','5'}]
        if len(constrained) > 1:
            raise ValueError('Ambiguous constrained comparison')
        if constrained:
            curves.append((constrained[0], 'Original transport and mean ventilation', '#8553a3'))
        if temperature == -5.3:
            smooth = [f for f in selected if f['temperature_C'] == temperature
                      and f['family'] == 'axis' and f['seed_diameter_um'] == 0.
                      and f['specification']['early_axis_multiplier_bounds'] == data['default_early_bounds']
                      and f['specification'].get('ventilation_max', 2.) == 2.
                      and not f['specification'].get('fixed')]
            if len(smooth) > 1:
                raise ValueError('Ambiguous warm comparison without compact cutoff')
            if smooth:
                curves.append((smooth[0], 'Warm refit without compact cutoff', '#b53d59'))
        for j, scale in enumerate((1e6, 1e3, 1e3)):
            ax = axes[i, j]
            for fit, label, color in curves:
                rows = [r for r in fit['assessment']['rows'] if r['initialization'] == initialization]
                draw(ax, [r['age_min'] for r in rows], [r['values'][j] for r in rows],
                     scale, label, color, '.-', [temperature, initialization, j, fit['name']])
            for k, segment in enumerate(previous['laboratory_plot_segments'][str(temperature)][data['quantities'][j]]):
                draw(ax, segment['age_min'], segment['values'], scale,
                     'Laboratory fitted history' if k == 0 else None, '#222222', '--',
                     [temperature, initialization, j, 'laboratory', k])
            ax.set_yscale('log')
            ax.grid(alpha=.18)
            if i == 0:
                ax.set_title(labels[j])
            if j == 0:
                ax.set_ylabel(f'{temperature:g} °C\n'+('4 µm seed' if initialization == 'sphere' else '3-minute restart'))
            if i == 3:
                ax.set_xlabel('Age (minutes)')
    legend = {}
    for ax in (axes[0, 0], axes[2, 0]):
        handles, labels = ax.get_legend_handles_labels()
        legend.update(zip(labels, handles))
    legend_order = ['Previous early limit fixed', 'Independent early and later limits',
                    'Warm refit without compact cutoff', 'Original transport and mean ventilation',
                    'Laboratory fitted history']
    legend_labels = [name for name in legend_order if name in legend]
    fig.legend([legend[name] for name in legend_labels], legend_labels,
               loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(.5, .969))
    fig.suptitle('Independent early and later growth: Earth diagnostic comparisons', y=.995)
    fig.text(.5, .012, 'Fit ages: 5, 10, 20 min. Other measurements excluded from the score. Separate fits at each temperature.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, .90))
    out = HERE/'generated'
    out.mkdir(parents=True, exist_ok=True)
    identity = dict(product_sha256=digest(PRODUCT), renderer_sha256=digest(Path(__file__)), evidence=data['evidence'])
    png, pdf = out/'growth_limits.png', out/'growth_limits.pdf'
    fig.savefig(png, dpi=160, metadata={'Description': json.dumps(identity, sort_keys=True)})
    fig.savefig(pdf, metadata={'Title': 'Early and later Earth crystal growth',
                              'Subject': json.dumps(identity, sort_keys=True), 'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    manifest = dict(schema='terluna.visualization.snow-growth-limits/1', **identity,
                    inputs=expected, reading_rule=data['reading_rule'], plotted_values=plotted,
                    outputs={p.name: digest(p) for p in (png, pdf)})
    (out/'growth_limits_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(curves=len(plotted), outputs=manifest['outputs']), indent=2))


if __name__ == '__main__':
    render()
