"""Storm lives in the electrified box's windows written every ten minutes (cm1_run.py's
box_0e_elec_corrected_storms_first and _second, the main run's busiest days run again from its restarts): each storm's
core followed from output to output, with its updraft, graupel and hail, charge, rain and flashes.

    python3 -m climate.crm.storm_lives box_0e_elec_corrected_storms_first box_0e_elec_corrected_storms_second
    python3 -m climate.crm.storm_lives supercell_elec        # Earth's benchmark storm, the same rules

writes ../results/crm/storms_<case>.json.

A core is the columns whose graupel and hail path reaches CORE_KG_M2, joined across sides and corners and, in a box with
periodic edges (the lunar boxes; the benchmark's are open), across them. A core continues the track of the core at the
output before that it overlaps most, when that core overlaps it most in turn; a core that overlaps none continues the
nearest earlier core within MATCH_KM whose track has not gone on. Where cores merge, the track with the largest overlap
carries on and the others end, merged into it; where one splits, the largest part carries on and the others begin new
tracks, split from it. Each flash belongs to the core, at the output nearest its time, that holds the column it starts
in, or failing that a column beside it.

The updraft, the cloud ice and snow and the charge count over the core's columns and those beside them, since a storm's
updraft stands at the edge of its precipitation and its charge spreads; the charge counts where the net density passes
0.1 nC/m3, as in elec_analysis.structure. The charging zone is -30 to -5 C (elec_analysis.CHARGING_ZONE_C), where
Saunders and Peck's law separates most of the charge. The summary is repeated with cores at each of
CORE_SENSITIVITY's paths, to show how far the threshold moves it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy import stats as sstats

from climate.crm import elec_analysis as ea
from climate.crm import ring_analysis as ra
from climate.crm.mixed_phase_analysis import significant, stats

SCHEMA = 'terluna.climate.crm-storm-lives/1'
CORE_KG_M2 = 1.0           # graupel and hail path of a core's columns
CORE_SENSITIVITY = (0.5, 2.0)
MATCH_KM = 20.0            # the farthest a core overlapping none may lie from an earlier core and continue its track
UPDRAFT_M_S = 5.0          # the updraft volume counts air rising at least this fast
CONDENSATE_KG_KG = 1.0e-5  # a core's top: the highest level with this much condensate
CHARGE_C_M3 = 0.1e-9       # as elec_analysis.structure
CHARGED_C = 10.0           # a core counts as charged once it holds this much charge of either sign
EIGHT = np.ones((3, 3), bool)
# the quantities that predict a flash rate, set against each core's flashes at each output
PREDICTORS = ('gh_zone_kg', 'gh_kg', 'ice_snow_zone_kg', 'ice_product_kg2', 'updraft_km3', 'w_max_m_s', 'top_km',
              'area_km2')
EVIDENCE = ('CM1 r22.0 with WRF-ELEC (climate/crm/cm1_elec.py) at the case\'s outputs. The lunar windows run the main '
            'run\'s busiest days again from its restarts with output every ten minutes; CM1 on several threads does '
            'not repeat itself bit for bit, so a window is a fresh realization of its days, whose storms part from the '
            'main run\'s within hours, and its flash counts stand beside the main run\'s without adding to them. '
            'Cores, tracks and the flashes\' assignment follow the rules in the module\'s description; times resolve '
            'to the outputs.')
READING_RULE = ('Hours count from the window\'s start (start_s). A track is one core followed through the outputs; '
                'whole marks the tracks that begin after the window\'s first output and end before its last and '
                'neither split from another track nor merged into one: storms whose own lives the window holds '
                'entire. Each output stands for the output_s around it, so a life is its outputs times output_s, '
                'and flashes per output are those nearest it in time. lead_h is from a new core\'s first output to '
                'its first flash, charge_lead_h from its first output holding 10 C of either sign, each counted '
                'from half an output before; the lags are between the outputs where its updraft, its graupel and '
                'hail in the charging zone, its flashes and the rain beneath it peak. The correlations are '
                'Spearman\'s, over every output of every track that flashed, of its flashes against each quantity '
                'there. sensitivity repeats the summary with cores at other graupel and hail paths.')


def periodic_labels(mask: np.ndarray, periodic: bool = True) -> np.ndarray:
    """Labels 1..n of the connected regions of a (ny, nx) mask, sides and corners joining, and across the edges when
    they are periodic."""
    labels, n = ndimage.label(mask, structure=EIGHT)
    if n == 0 or not periodic:
        return labels
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    ny, nx = mask.shape
    pairs = [(labels[j, 0], labels[(j + d) % ny, nx - 1]) for j in range(ny) for d in (-1, 0, 1)]
    pairs += [(labels[0, i], labels[ny - 1, (i + d) % nx]) for i in range(nx) for d in (-1, 0, 1)]
    for a, b in pairs:
        if a and b:
            ra_, rb = find(a), find(b)
            if ra_ != rb:
                parent[max(ra_, rb)] = min(ra_, rb)
    roots = np.array([find(a) for a in range(n + 1)])
    compact = np.unique(roots, return_inverse=True)[1].reshape(-1)      # the background's root, 0, stays 0
    return compact[labels]


def widen(labels: np.ndarray, periodic: bool = True) -> np.ndarray:
    """The labels with each unlabelled column beside a core (sides and corners, and across periodic edges) given that
    core's label, the larger where two cores meet."""
    near = ndimage.maximum_filter(labels, size=3, mode='wrap' if periodic else 'constant')
    return np.where(labels > 0, labels, near)


def periodic_centroid(jj, ii, weights, ny: int, nx: int, periodic: bool = True) -> tuple:
    """The weighted centre (row, column) of columns, on a periodic grid from the mean of their angles."""
    if not periodic:
        return float(np.average(jj, weights=weights)), float(np.average(ii, weights=weights))
    out = []
    for idx, n in ((jj, ny), (ii, nx)):
        a = 2.0 * np.pi * np.asarray(idx, float) / n
        c = np.arctan2(np.sum(weights * np.sin(a)), np.sum(weights * np.cos(a)))
        out.append(float((c * n / (2.0 * np.pi)) % n))
    return tuple(out)


def periodic_distance(a: tuple, b: tuple, ny: int, nx: int, periodic: bool = True) -> float:
    """Distance in columns between two (row, column) points, the shorter way round on a periodic grid."""
    dj, di = abs(a[0] - b[0]), abs(a[1] - b[1])
    if periodic:
        dj, di = min(dj, ny - dj), min(di, nx - di)
    return float(np.hypot(dj, di))


def periodic_case(case: Path) -> bool:
    """Whether a case's box has periodic edges on all four sides (CM1's boundary condition 1)."""
    from climate.crm.cm1_run import namelist_value
    text = (case / 'namelist.input').read_text()
    return all(namelist_value(text, 'param2', side) == '1' for side in ('wbc', 'ebc', 'sbc', 'nbc'))


def column_fields(snap: dict, s: dict, zh: np.ndarray, area_m2: float) -> dict:
    """Per column (flattened as the snapshot is) what the cores sum or take the largest of. snap is the raw snapshot,
    s elec_analysis.snapshot_charge's reading of it."""
    vol = s['volume']
    gh = s['graupel']
    tz = s['t'] - 273.15
    zone = (tz >= ea.CHARGING_ZONE_C[0]) & (tz <= ea.CHARGING_ZONE_C[1])
    ice_snow = (snap['qi'] + snap['qs']) * s['rho']
    condensate = sum(snap[q] for q in ('qc', 'qr', 'qi', 'qs', 'qg', 'qhl'))
    has = condensate >= CONDENSATE_KG_KG
    top = has.shape[0] - 1 - np.argmax(has[::-1], axis=0)
    net, w = s['net'], s['w']
    return dict(path_kg_m2=(gh * vol).sum(axis=0) / area_m2, gh_kg=(gh * vol).sum(axis=0),
                gh_zone_kg=(gh * zone * vol).sum(axis=0), ice_snow_zone_kg=(ice_snow * zone * vol).sum(axis=0),
                w_max_m_s=w.max(axis=0), updraft_m3=((w >= UPDRAFT_M_S) * vol).sum(axis=0),
                top_km=np.where(has.any(axis=0), zh[top], 0.0), rain_mm_h=snap['prate'] * 3600.0,
                positive_c=(np.where(net > CHARGE_C_M3, net, 0.0) * vol).sum(axis=0),
                negative_c=(np.where(net < -CHARGE_C_M3, net, 0.0) * vol).sum(axis=0))


def core_table(fields: dict, labels: np.ndarray, dx: float, dy: float, periodic: bool = True) -> list:
    """Each core's size, place and contents at one output (fields from column_fields, labels from periodic_labels)."""
    ny, nx = labels.shape
    flat, near_flat = labels.reshape(-1), widen(labels, periodic).reshape(-1)
    out = []
    for c in range(1, int(labels.max()) + 1):
        inside, near = flat == c, near_flat == c
        jj, ii = np.divmod(np.flatnonzero(inside), nx)
        cj, ci = periodic_centroid(jj, ii, fields['path_kg_m2'][inside], ny, nx, periodic)
        gh_zone, ice_zone = float(fields['gh_zone_kg'][inside].sum()), float(fields['ice_snow_zone_kg'][near].sum())
        out.append(dict(
            columns=int(inside.sum()), area_km2=float(inside.sum() * dx * dy / 1.0e6), row=cj, column=ci,
            x_km=(ci + 0.5) * dx / 1000.0, y_km=(cj + 0.5) * dy / 1000.0,
            gh_kg=float(fields['gh_kg'][inside].sum()), gh_zone_kg=gh_zone, ice_snow_zone_kg=ice_zone,
            ice_product_kg2=gh_zone * ice_zone, w_max_m_s=float(fields['w_max_m_s'][near].max()),
            updraft_km3=float(fields['updraft_m3'][near].sum() / 1.0e9), top_km=float(fields['top_km'][inside].max()),
            rain_mm_h=float(fields['rain_mm_h'][inside].max()), positive_c=float(fields['positive_c'][near].sum()),
            negative_c=float(fields['negative_c'][near].sum())))
    return out


def link(prev: np.ndarray, cur: np.ndarray, prev_cores: list, cur_cores: list, match_columns: float,
         periodic: bool = True) -> tuple:
    """How the cores at one output (cur) follow those at the output before (prev): for each current core label, the
    earlier label whose track it continues, or the earlier label it split from; for each earlier label whose track
    ends in a merger, the current label it merged into."""
    ny, nx = cur.shape
    n_prev, n_cur = int(prev.max()), int(cur.max())
    overlap = np.zeros((n_prev + 1, n_cur + 1), int)
    np.add.at(overlap, (prev.reshape(-1), cur.reshape(-1)), 1)
    overlap[0, :] = 0
    overlap[:, 0] = 0
    best_cur, best_prev = overlap.argmax(axis=1), overlap.argmax(axis=0)
    continues, split_from, merged_into = {}, {}, {}
    for b in range(1, n_cur + 1):
        a = int(best_prev[b])
        if overlap[a, b] == 0:
            continue
        if best_cur[a] == b:
            continues[b] = a
        else:
            split_from[b] = a
    for a in range(1, n_prev + 1):
        b = int(best_cur[a])
        if overlap[a, b] > 0 and continues.get(b) != a:
            merged_into[a] = b
    free = [a for a in range(1, n_prev + 1) if a not in continues.values() and a not in merged_into]
    for b in range(1, n_cur + 1):
        if b in continues or b in split_from or not free:
            continue
        here = (cur_cores[b - 1]['row'], cur_cores[b - 1]['column'])
        d, a = min((periodic_distance((prev_cores[a - 1]['row'], prev_cores[a - 1]['column']), here, ny, nx,
                                      periodic), a) for a in free)
        if d <= match_columns:
            continues[b] = a
            free.remove(a)
    return continues, split_from, merged_into


def follow(label_maps: list, core_lists: list, match_columns: float, periodic: bool = True) -> list:
    """Tracks through the outputs: each a list of (output, core label) with the track it split from and the one it
    merged into, as numbers in the returned list."""
    tracks, current = [], {}
    for k, (labels, cores) in enumerate(zip(label_maps, core_lists)):
        if k == 0:
            continues, split_from, merged_into = {}, {}, {}
        else:
            continues, split_from, merged_into = link(label_maps[k - 1], labels, core_lists[k - 1], cores,
                                                      match_columns, periodic)
        nxt = {}
        for b in range(1, len(cores) + 1):
            if b in continues:
                t = current[continues[b]]
                tracks[t]['outputs'].append(k)
                tracks[t]['labels'].append(b)
            else:
                tracks.append(dict(outputs=[k], labels=[b], merged_into=None,
                                   split_from=current[split_from[b]] if b in split_from else None))
                t = len(tracks) - 1
            nxt[b] = t
        for a, b in merged_into.items():
            tracks[current[a]]['merged_into'] = nxt[b]
        current = nxt
    return tracks


def assign_flashes(f: dict, out_times, label_maps: list, dx: float, dy: float, periodic: bool = True) -> tuple:
    """For each flash, the output nearest its time and the core there holding its starting column, or failing that a
    column beside it (0 when none does)."""
    t = np.asarray(out_times, float)
    if not f['time_s'].size:
        return np.zeros(0, int), np.zeros(0, int)
    k = np.abs(f['time_s'][:, None] - t[None, :]).argmin(axis=1)
    ny, nx = label_maps[0].shape
    i = (f['x_m'] // dx).astype(int) % nx
    j = (f['y_m'] // dy).astype(int) % ny
    wide = {}
    core = []
    for kk, jj, ii in zip(k, j, i):
        c = int(label_maps[kk][jj, ii])
        if not c:
            if kk not in wide:
                wide[kk] = widen(label_maps[kk], periodic)
            c = int(wide[kk][jj, ii])
        core.append(c)
    return k, np.array(core, int)


def track_record(track: dict, core_lists: list, hours: np.ndarray, flashes: dict, n_out: int, output_h: float) -> dict:
    """One track's series and the moments of its life (hours from the window's start)."""
    ks = track['outputs']
    cores = [core_lists[k][b - 1] for k, b in zip(ks, track['labels'])]
    keys = ('area_km2', 'x_km', 'y_km', 'gh_kg', 'gh_zone_kg', 'ice_snow_zone_kg', 'w_max_m_s', 'updraft_km3', 'top_km',
            'rain_mm_h', 'positive_c', 'negative_c')
    series = {key: [c[key] for c in cores] for key in keys}
    series['hour'] = [float(hours[k]) for k in ks]
    series['flashes'] = [int(flashes['count'].get((k, b), 0)) for k, b in zip(ks, track['labels'])]
    series['ground_strikes'] = [int(flashes['ground'].get((k, b), 0)) for k, b in zip(ks, track['labels'])]
    times = sorted(t for (k, b) in zip(ks, track['labels']) for t in flashes['times'].get((k, b), []))
    peak = lambda key: (float(np.max(series[key])), series['hour'][int(np.argmax(series[key]))])
    rec = dict(first_h=series['hour'][0], last_h=series['hour'][-1], life_h=len(ks) * output_h,
               whole=bool(ks[0] > 0 and ks[-1] < n_out - 1 and track['split_from'] is None
                          and track['merged_into'] is None), split_from=track['split_from'],
               merged_into=track['merged_into'], flashes=int(sum(series['flashes'])),
               ground_strikes=int(sum(series['ground_strikes'])))
    for key in ('w_max_m_s', 'gh_zone_kg', 'gh_kg', 'top_km', 'area_km2', 'positive_c', 'rain_mm_h'):
        rec[f'peak_{key}'], rec[f'peak_{key}_h'] = peak(key)
    rec['peak_negative_c'] = float(np.min(series['negative_c']))
    if times:
        rec.update(first_flash_h=times[0], last_flash_h=times[-1], flash_hours=times[-1] - times[0],
                   peak_flashes_per_min=float(np.max(series['flashes']) / (output_h * 60.0)),
                   peak_flashes_h=series['hour'][int(np.argmax(series['flashes']))],
                   after_last_flash_h=series['hour'][-1] + 0.5 * output_h - times[-1])
        if track['split_from'] is None and ks[0] > 0:
            rec['lead_h'] = times[0] - (series['hour'][0] - 0.5 * output_h)
            charged = [h for h, p, n in zip(series['hour'], series['positive_c'], series['negative_c'])
                       if max(p, -n) >= CHARGED_C]
            if charged:
                rec['charge_lead_h'] = times[0] - (charged[0] - 0.5 * output_h)
        rec['graupel_after_updraft_h'] = round(rec['peak_gh_zone_kg_h'] - rec['peak_w_max_m_s_h'], 6)
        rec['flashes_after_graupel_h'] = round(rec['peak_flashes_h'] - rec['peak_gh_zone_kg_h'], 6)
        rec['rain_after_flashes_h'] = round(rec['peak_rain_mm_h_h'] - rec['peak_flashes_h'], 6)
    rec['series'] = series
    return rec


def spread(values) -> dict:
    """mixed_phase_analysis.stats with the smallest value as well."""
    out = stats(values)
    if out['n']:
        out['min'] = float(np.min(np.asarray(values, float)))
    return out


def summarize(tracks: list, flashes_total: int, flashes_assigned: int) -> dict:
    """The window's tracks together: how long storms live, how soon and how long they flash, what peaks before what,
    and which quantity follows the flashes best."""
    flashing = [t for t in tracks if t['flashes']]
    whole = [t for t in tracks if t['whole']]
    pick = lambda group, key: [t[key] for t in group if key in t]
    out = dict(tracks=len(tracks), whole_tracks=len(whole), flashing_tracks=len(flashing),
               flashes=flashes_total, flashes_in_cores=flashes_assigned,
               ground_strikes_in_cores=int(sum(t['ground_strikes'] for t in tracks)),
               life_h=dict(flashing=spread(pick([t for t in whole if t['flashes']], 'life_h')),
                           without_flashes=spread(pick([t for t in whole if not t['flashes']], 'life_h'))),
               lead_h=spread(pick(flashing, 'lead_h')), charge_lead_h=spread(pick(flashing, 'charge_lead_h')),
               flash_hours=spread(pick(flashing, 'flash_hours')),
               after_last_flash_h=spread(pick([t for t in flashing if t['whole']], 'after_last_flash_h')),
               flashes_per_track=spread(pick(flashing, 'flashes')),
               peak_flashes_per_min=spread(pick(flashing, 'peak_flashes_per_min')),
               graupel_after_updraft_h=spread(pick(flashing, 'graupel_after_updraft_h')),
               flashes_after_graupel_h=spread(pick(flashing, 'flashes_after_graupel_h')),
               rain_after_flashes_h=spread(pick(flashing, 'rain_after_flashes_h')))
    out['peaks'] = {group: {key: spread(pick(members, f'peak_{key}'))
                            for key in ('w_max_m_s', 'gh_zone_kg', 'top_km', 'area_km2', 'positive_c')}
                    for group, members in (('flashing', flashing),
                                           ('without_flashes', [t for t in tracks if not t['flashes']]))}
    y, x = [], {key: [] for key in PREDICTORS}
    for t in flashing:
        s = t['series']
        y += s['flashes']
        for key in PREDICTORS:
            x[key] += ([g * c for g, c in zip(s['gh_zone_kg'], s['ice_snow_zone_kg'])] if key == 'ice_product_kg2'
                       else s[key])
    corr = {}
    if len(y) > 2 and np.ptp(y) > 0:
        for key in PREDICTORS:
            if np.ptp(x[key]) > 0.0:
                corr[key] = float(sstats.spearmanr(x[key], y).statistic)
        corr['outputs'] = len(y)
    out['flashes_spearman'] = corr
    return out


def tracked(fields_list: list, f: dict, out_times: np.ndarray, start: float, grid: tuple, output_h: float,
            threshold: float, periodic: bool) -> dict:
    """Cores at one threshold through every output, their tracks and the flashes each holds."""
    ny, nx, dx, dy = grid
    label_maps = [periodic_labels(fl['path_kg_m2'].reshape(ny, nx) >= threshold, periodic) for fl in fields_list]
    core_lists = [core_table(fl, labels, dx, dy, periodic) for fl, labels in zip(fields_list, label_maps)]
    k_of, core_of = assign_flashes(f, out_times, label_maps, dx, dy, periodic)
    ground = np.isin(f['kind'], (2, 3))
    flashes = dict(count={}, ground={}, times={})
    for k, c, g, t in zip(k_of, core_of, ground, f['time_s']):
        flashes['count'][(k, c)] = flashes['count'].get((k, c), 0) + 1
        flashes['ground'][(k, c)] = flashes['ground'].get((k, c), 0) + int(g)
        flashes['times'].setdefault((k, c), []).append(float((t - start) / 3600.0))
    hours = (out_times - start) / 3600.0
    tracks = follow(label_maps, core_lists, MATCH_KM * 1000.0 / dx, periodic)
    records = [track_record(t, core_lists, hours, flashes, len(out_times), output_h) for t in tracks]
    return dict(label_maps=label_maps, k_of=k_of, core_of=core_of, records=records,
                summary=summarize(records, int(f['time_s'].size), int(np.sum(core_of > 0))))


def analyse(name: str) -> dict:
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    grid = record['grid']
    nx, ny = grid['nx'], grid.get('ny', grid['nx'])
    dx, dy = grid['dx_m'], grid.get('dy_m', grid['dx_m'])
    periodic = periodic_case(case)
    zh = ea.heights_km(case / 'cm1out_s.ctl')
    zw = (np.loadtxt(case / 'input_grid_z') if (case / 'input_grid_z').exists()      # evenly spaced levels have none
          else np.concatenate([[0.0], 1000.0 * (zh[:-1] + zh[1:]) / 2.0,
                               [1000.0 * (2 * zh[-1] - (zh[-1] + zh[-2]) / 2.0)]]))
    progress = json.loads((case / 'progress.json').read_text())
    start = float(progress.get('start_s', 0.0))
    times = ea.output_times(case)
    outputs = sorted(n for n in (int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat')) if n in times)
    output_s = record['configuration']['output_s'] * record.get('time_scale', 1.0)
    fields_list, domain = [], []
    for n in outputs:
        snap = ra.read_snapshot(case, n)
        fields = column_fields(snap, ea.snapshot_charge(case, n, dx, dy, zw, snap=snap), zh, dx * dy)
        fields_list.append(fields)
        domain.append(dict(hour=(times[n] - start) / 3600.0, w_max_m_s=float(fields['w_max_m_s'].max()),
                           gh_kg=float(fields['gh_kg'].sum()), gh_zone_kg=float(fields['gh_zone_kg'].sum()),
                           positive_c=float(fields['positive_c'].sum()), negative_c=float(fields['negative_c'].sum())))
    out_times = np.array([times[n] for n in outputs])
    f = ea.read_log(case / 'terluna_flashes.txt', ea.FLASH_COLUMNS)
    f = {k: v[f['time_s'] > start] for k, v in f.items()}
    at = lambda threshold: tracked(fields_list, f, out_times, start, (ny, nx, dx, dy), output_s / 3600.0, threshold,
                                   periodic)
    base = at(CORE_KG_M2)
    for k, row in enumerate(domain):
        here = base['k_of'] == k
        row.update(cores=int(base['label_maps'][k].max()), flashes=int(here.sum()),
                   flashing_cores=int(np.unique(base['core_of'][here & (base['core_of'] > 0)]).size))
    return dict(
        schema=SCHEMA, case=name, purpose=record['configuration'].get('purpose'),
        producer=dict(domain='climate', files={'crm/storm_lives.py': hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest()[:16]}),
        evidence=EVIDENCE, reading_rule=READING_RULE, start_s=start, output_s=output_s, periodic=periodic,
        settings=dict(core_kg_m2=CORE_KG_M2, match_km=MATCH_KM, updraft_m_s=UPDRAFT_M_S,
                      condensate_kg_kg=CONDENSATE_KG_KG, charge_c_m3=CHARGE_C_M3,
                      charging_zone_c=list(ea.CHARGING_ZONE_C)),
        executables=sorted({seg.get('executable') for seg in progress['segments']}),
        summary=base['summary'],
        sensitivity={f'core_{threshold:g}_kg_m2': at(threshold)['summary'] for threshold in CORE_SENSITIVITY},
        domain=domain, tracks=base['records'])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('cases', nargs='+')
    args = parser.parse_args(argv)
    for name in args.cases:
        out = significant(analyse(name))
        path = ra.RESULTS / f'storms_{name}.json'
        path.write_text(json.dumps(out, indent=1) + '\n')
        s = out['summary']
        print(f"{name}: {s['tracks']} tracks ({s['whole_tracks']} whole), {s['flashing_tracks']} flashing; "
              f"{s['flashes_in_cores']} of {s['flashes']} flashes in cores -> {path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
