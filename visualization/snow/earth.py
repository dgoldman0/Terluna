"""Compare saved author-rerun Earth profiles with hash-verified raw archives.

    OPENBLAS_NUM_THREADS=1 python visualization/snow/earth.py --archive-dir PATH

This renderer imports no domain model and performs no physical calculation.
It reconstructs the product's profile differences, then draws the same 20 m bins.
External archive bytes stay external; generated images/manifests are ignored.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import scipy
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.welss-author-comparison/1'
CASES = (
    ('growth', 'welss_author_growth.json', 'ob_hp0.nc', False, 'Deposition + settling'),
    ('aggregation', 'welss_author_aggregation.json', 'ob_agg0.nc', True,
     'Deposition + settling + aggregation'),
)
MOMENTS = ('Nd', 'Md', 'Fn', 'Fm')
SCALES = dict(Nd=.001, Md=1000., Fn=.001, Fm=1000.)
LABELS = dict(Nd='Number concentration, $N_d$\n(L⁻¹)',
              Md='Ice mass concentration, $M_d$\n(g m⁻³)',
              Fn='Number velocity moment, $F_n$\n(10³ m⁻² s⁻¹)',
              Fm='Mass velocity moment, $F_m$\n(g m⁻² s⁻¹)')
INK, MUTED, GRID = '#172b3a', '#52646d', '#dce3e6'
ARCHIVE, RERUN = '#667681', '#206ca7'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encoded_identity(value):
    """Exact serialization used for the producer's embedded run/build records."""
    return hashlib.sha256((json.dumps(value, indent=2, allow_nan=False)+'\n').encode()).hexdigest()


def verify(actual, expected, description):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(description+': metric keys differ')
        for key in expected:
            verify(actual[key], expected[key], description+'/'+key)
    elif expected is None:
        if actual is not None:
            raise ValueError(description+': expected missing region')
    elif not np.isclose(actual, expected, rtol=2e-12, atol=1e-14):
        raise ValueError(f'{description}: {actual} differs from saved {expected}')


def profile_metrics(model, reference, height):
    """Reconstruct the saved comparison on its exact valid-bin masks."""
    valid = np.isfinite(model) & np.isfinite(reference) & (reference > 0) & (height < 3000.)
    result = {}
    for region, mask in (('all_populated', valid), ('below_nucleation', valid & (height < 2700.)),
                         ('nucleation_belt', valid & (height >= 2700.))):
        if not np.any(mask):
            result[region] = None
            continue
        delta, ref = model[mask]-reference[mask], reference[mask]
        fractional = delta/ref
        result[region] = dict(count=int(mask.sum()),
            relative_L1=float(np.abs(delta).sum()/ref.sum()),
            relative_bias=float(delta.sum()/ref.sum()),
            maximum_absolute_relative_difference=float(np.max(np.abs(fractional))),
            maximum_difference_height_m=float(height[mask][np.argmax(np.abs(fractional))]))
    return result


def read_archive(raw):
    """Decode total moments and two aggregate-only diagnostics with actual missing masks."""
    with netcdf_file(io.BytesIO(raw), mmap=False) as source:
        arrays = {}
        for name in ('height', *MOMENTS, 'Nd_unr', 'Fm_unr'):
            variable = source.variables[name]
            values = np.asarray(variable[:], float).copy()
            for attr in ('_FillValue', 'missing_value'):
                if hasattr(variable, attr):
                    values[values == getattr(variable, attr)] = np.nan
            arrays[name] = values
    return arrays


def prepare(archive_dir, root=ROOT):
    rows, products, inputs, verified = [], {}, {}, {}
    for key, filename, member, aggregation, title in CASES:
        path = root/'atmosphere/snow/results'/filename
        raw = path.read_bytes()
        product = json.loads(raw)
        inputs[str(path.relative_to(root))] = hashlib.sha256(raw).hexdigest()
        if product.get('schema') != SCHEMA or not product.get('evidence') or not product.get('reading_rule'):
            raise ValueError(filename+': missing supported schema/evidence/reading rule')
        for name, expected in product['producer'].items():
            if Path(name).name != name:
                raise ValueError('unexpected producer path: '+name)
            source = root/'atmosphere/snow'/name
            actual = digest(source)
            if actual != expected:
                raise ValueError('producer identity differs: '+name)
            verified[str(source.relative_to(root))] = actual
        run = product['run']
        if run['producer'] != product['producer'] or encoded_identity(run) != product['run_record_sha256']:
            raise ValueError(filename+': embedded run identity differs')
        if run['aggregation'] is not aggregation or run['returncode'] != 0:
            raise ValueError(filename+': case or execution state differs')
        build = run['build']
        if encoded_identity(build) != run['build_record_sha256']:
            raise ValueError(filename+': embedded build identity differs')
        for name, expected in build['producer'].items():
            if product['producer'].get(name) != expected:
                raise ValueError(filename+': build producer differs')
        source_manifest = json.loads((root/'atmosphere/snow/welss_aggregation_sources.json').read_text())
        if build['revision'] != source_manifest['revision'] or build['source_files'] != {
                item['path']: item['sha256'] for item in source_manifest['files']}:
            raise ValueError(filename+': pinned build source list differs')
        reference = product['reference']
        if reference['member'] != member:
            raise ValueError(filename+': wrong archive member')
        archive_bytes = (archive_dir/member).read_bytes()
        archive_hash = hashlib.sha256(archive_bytes).hexdigest()
        if archive_hash != reference['sha256'] or archive_hash != source_manifest['archive']['members'][member]:
            raise ValueError(member+': archive identity differs')
        inputs['external_archive/'+member] = archive_hash
        archive = read_archive(archive_bytes)
        height = np.asarray(product['height_m'], float)
        if not np.array_equal(height, archive['height']) or not np.array_equal(height, np.arange(251)*20.):
            raise ValueError(filename+': expected common 20 m bins')
        mode = product['profiles']['last5h_31']
        if mode['snapshot_count'] != 31 or mode['snapshot_times_s'] != list(range(18000,36001,600)):
            raise ValueError(filename+': expected 31 inclusive last-five-hour snapshots')
        plotted = height < 3000.
        moments = {}
        for name in MOMENTS:
            model = np.asarray(mode['raw_profiles'][name], float)
            if model.shape != height.shape or np.any(~np.isfinite(model)) or np.any(model < 0):
                raise ValueError(filename+'/'+name+': invalid raw profile')
            ref = archive[name]
            metrics = profile_metrics(model, ref, height)
            verify(metrics, mode['raw_comparison'][name], filename+'/'+name)
            moments[name] = dict(author_raw_si=model[plotted].tolist(),
                archive_raw_si=[None if not np.isfinite(x) else float(x) for x in ref[plotted]],
                display_factor=SCALES[name], label=LABELS[name], reconstructed_comparison=metrics)
        categories = {}
        if aggregation:
            for name in ('Nd_unr', 'Fm_unr'):
                model = np.asarray(mode['raw_profiles'][name], float)
                metrics = profile_metrics(model, archive[name], height)
                verify(metrics, mode['raw_comparison'][name], filename+'/'+name)
                categories[name] = metrics
        rows.append(dict(id=key, title=title, archive=member, height_lower_edge_m=height[plotted].tolist(),
                         bin_width_m=20., moments=moments, unplotted_aggregate_category_comparison=categories))
        products[key] = product
    a, b = products['growth']['run'], products['aggregation']['run']
    if a['build']['executable_sha256'] != b['build']['executable_sha256']:
        raise ValueError('the rows are not from the same author executable')
    return rows, products, inputs, verified


def step_profile(axis, values, lower, color, **style):
    """Horizontal coordinate is a cell average, vertical edges are explicit."""
    axis.plot(np.repeat(values, 2), np.column_stack((lower, lower+20.)).ravel()/1000.,
              color=color, **style)


def render(rows, output, identity):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'text.color': INK, 'axes.labelcolor': INK, 'svg.hashsalt': 'terluna-earth-profiles-v1'})
    fig, axes = plt.subplots(2, 4, figsize=(14.2, 9.2), sharey=True)
    fig.subplots_adjust(left=.065, right=.98, bottom=.21, top=.78, hspace=.45, wspace=.24)
    fig.text(.065, .95, 'Earth snowfall totals: declared author rerun and archive', fontsize=20, weight='bold')
    fig.text(.065, .912, 'Raw 20 m bins · 31 snapshots from model hours 5–10 inclusive · no vertical smoothing',
             fontsize=11, color=MUTED)
    handles = [Line2D([], [], color=ARCHIVE, lw=1.3, marker='.', markersize=5, label='2023 archived author profile'),
               Line2D([], [], color=RERUN, lw=1.8, label='Pinned public-source author rerun'),
               Patch(facecolor='#e8edf1', label='Prescribed nucleation layer')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.059, .888), frameon=False, ncol=3,
               fontsize=10, columnspacing=2.0)
    for col, name in enumerate(MOMENTS):
        maximum = max(np.nanmax(np.asarray(row['moments'][name][kind], float))*SCALES[name]
                      for row in rows for kind in ('author_raw_si', 'archive_raw_si'))
        for r, row in enumerate(rows):
            axis = axes[r, col]; moment = row['moments'][name]
            z = np.asarray(row['height_lower_edge_m'])
            archive = np.asarray(moment['archive_raw_si'], float)*SCALES[name]
            model = np.asarray(moment['author_raw_si'])*SCALES[name]
            axis.axhspan(2.7, 3., color='#e8edf1', zorder=0)
            step_profile(axis, archive, z, ARCHIVE, lw=1.1, alpha=.85)
            axis.plot(archive, (z+10.)/1000., '.', color=ARCHIVE, ms=2.5, alpha=.7)
            step_profile(axis, model, z, RERUN, lw=1.6)
            axis.set_xlim(0, maximum*1.08)
            axis.set_ylim(0, 3.)
            axis.xaxis.set_major_locator(MaxNLocator(4))
            axis.set_yticks(np.arange(0,3.01,.5))
            axis.set_xlabel(LABELS[name], fontsize=9.5)
            if col == 0:
                axis.set_ylabel('Height (km)')
            axis.spines[['top','right']].set_visible(False)
            axis.spines[['left','bottom']].set_color('#aab8bf')
            axis.tick_params(labelsize=9, colors=MUTED)
            axis.grid(axis='y', color=GRID, lw=.6)
            axis.set_axisbelow(True)
            error=moment['reconstructed_comparison']['all_populated']['relative_L1']*100
            tx, ty, align = (.03, .80, 'left') if name == 'Fn' else (.97, .80 if name in ('Md', 'Fm') else .055, 'right')
            axis.text(tx, ty, f'L1 difference: {error:.2f}%', transform=axis.transAxes,
                      ha=align, va='bottom', fontsize=9, color=RERUN,
                      bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=2))
    for r, row in enumerate(rows):
        pos=axes[r,0].get_position()
        fig.text(pos.x0, pos.y1+.018, f'{row["title"]}  ·  {row["archive"][:-3]}',
                 fontsize=12, weight='bold')
    fig.text(.065, .105, 'L1 = Σ |rerun − archive| / Σ archive over populated bins below 3 km; each column uses the same horizontal scale in both rows.',
             fontsize=9.3, color=MUTED)
    fig.text(.065, .078, 'Fn and Fm are cell-average Nv and Nmv velocity moments, not measured boundary crossings.',
             fontsize=9.3, color=MUTED)
    categories = rows[1]['unplotted_aggregate_category_comparison']
    nd = 100*categories['Nd_unr']['all_populated']['relative_L1']
    fm = 100*categories['Fm_unr']['all_populated']['relative_L1']
    fig.text(.065, .053, f'Totals shown. Aggregate-only L1 differences: Nd {nd:.1f}%, Fm {fm:.1f}%; full category comparison: welss_aggregation.md.',
             fontsize=9.3, color=MUTED)
    fig.text(.065, .027, 'Archive producer unidentified; this is an author-solver comparison, not independent population validation. One seed; convergence remains untested.',
             fontsize=9.3)
    output.mkdir(parents=True, exist_ok=True)
    paths=[output/'earth_comparison.svg', output/'earth_comparison.png']
    description=json.dumps(identity, sort_keys=True)
    fig.savefig(paths[0], metadata={'Date': None, 'Title': 'Declared-source Earth snowfall profile comparison',
                                  'Description': description})
    fig.savefig(paths[1], dpi=180, metadata={'Description': description})
    plt.close(fig)
    return {p.name:digest(p) for p in paths}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE/'generated')
    args=parser.parse_args()
    rows, products, inputs, verified=prepare(args.archive_dir)
    identity=dict(renderer={'visualization/snow/earth.py':digest(__file__)}, inputs=inputs)
    outputs=render(rows, args.output, identity)
    manifest=dict(schema='terluna.visualization.snow-earth/1', **identity,
        verified_producer_files=verified,
        environment=dict(matplotlib=matplotlib.__version__, numpy=np.__version__, scipy=scipy.__version__),
        upstream={key:{name:p[name] for name in ('schema','producer','evidence','reading_rule','caveats','reference')}
                  for key,p in products.items()},
        plotted=rows,
        display_rule='Raw cell means are drawn as vertical constant-value segments on their actual 20 m bins. Archive missing values remain gaps; no smoothing or curve interpolation. All shown bins lie below 3000 m. Column axes match between the two cases.',
        category_scope='Panels show total populations. The manifest separately reconstructs aggregate-only Nd and Fm differences, which are substantially larger; see atmosphere/snow/welss_aggregation.md for all category comparisons.',
        source_scope='The input products verify the original run files. This renderer verifies embedded run/build records, repository producer hashes and the independently read archive bytes; it does not require the external executable or raw run directory.',
        outputs=outputs)
    manifest_path=args.output/'earth_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(manifest=str(manifest_path), outputs=outputs,
        relative_L1={row['id']:{k:row['moments'][k]['reconstructed_comparison']['all_populated']['relative_L1']
                               for k in MOMENTS} for row in rows}), indent=2))


if __name__ == '__main__':
    main()
