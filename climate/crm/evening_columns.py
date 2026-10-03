"""Export second-cycle cloud histories across the corrected regional experiments.

python -m climate.crm.evening_columns --output research/runs/optical_comfort/evening
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from climate.crm.cloud_columns import dry_density, read_cm1_constants, sha256, SPECIES
from climate.crm.ring_analysis import case_geometry, hour_angle, read_snapshot
from climate.crm.cm1_run import RUNS
from shared.constants import SYNODIC_MONTH_DAYS

SCHEMA = 'terluna.climate.regional-evening-clouds/1'
CASES = ('ring_equator', 'ring_70_45e', 'ring_70_135e', 'box_highland_own_height')
THRESHOLD = 1e-5


def snapshot_numbers(record):
    dt = record['configuration']['output_s'] / 86400
    return np.arange(int(np.ceil(SYNODIC_MONTH_DAYS/dt))+1,
                     int(np.floor(min(record['configuration']['days'], 2*SYNODIC_MONTH_DAYS)/dt))+2)


def layer_edges(centres, ground, top):
    """Bound terrain-following mass cells by adjacent centre midpoints."""
    return np.vstack((ground[None, :], (centres[1:]+centres[:-1])/2,
                      np.full((1, centres.shape[1]), top)))


def column_summary(fields, density, centres, edges):
    cloud = fields['qc'] + fields['qi'] >= THRESHOLD
    dz = np.diff(edges, axis=0)
    if np.any(dz <= 0):
        raise ValueError('Nonpositive layer thickness')
    top = np.max(np.where(cloud, centres, 0), axis=0)
    base = np.min(np.where(cloud, centres, np.inf), axis=0)
    base[~cloud.any(axis=0)] = 0
    out = dict(cloud_top_m=top, cloud_base_m=base)
    for species in SPECIES:
        out[f'{species}_path_kg_m2'] = (density*fields[species]*dz).sum(axis=0)
    for floor in (1000, 20000, 40000, 60000):
        overlap = np.maximum(0, edges[1:]-np.maximum(edges[:-1], floor))
        out[f'cloud_path_above_{floor}_kg_m2'] = (density*(fields['qc']+fields['qi'])*overlap).sum(axis=0)
    out['low_cloud'] = np.any(cloud & (centres-edges[0] < 2000), axis=0)
    return out


def export_case(name, output):
    case = RUNS/name
    record = json.loads((case/'case.json').read_text())
    if record['reference']['run'] != 'A28_dim5_moon' or record.get('time_scale', 1.) != 1.:
        raise ValueError('Corrected full-length lunar forcing required')
    constants, sources = read_cm1_constants(record)
    ring = record['configuration'].get('kind') == 'tilted'
    if ring:
        geo = case_geometry(case)
        lat, lon = [np.array(record['path'][k]) for k in ('lat', 'lon')]
        land = geo['land']
        centres = np.broadcast_to(geo['zh'][:, None], (len(geo['zh']), len(lat)))
        edges = np.broadcast_to(geo['zw'][:, None], (len(geo['zw']), len(lat)))
    rows, times, snapshots, raw = {}, [], [], []
    numbers = snapshot_numbers(record)
    for count, n in enumerate(numbers):
        d = read_snapshot(case, int(n))
        rho, _ = dry_density(d['prs'], d['th'], d['qv'], constants)
        for species in SPECIES:
            if np.any(d[species] < 0) or not np.isfinite(d[species]).all():
                raise ValueError(f'Invalid {species} in {name}/{n}')
        day = (n-1)*record['configuration']['output_s']/86400
        if not ring:
            centres = d['zhval'].astype(float)
            ground = d['zs'].astype(float)
            edges = layer_edges(centres, ground, record['grid']['ztop_m'])
            lat = np.full(centres.shape[1], record['site']['lat_deg'])
            lon = np.full(len(lat), record['site']['lon_deg'])
            land = d['xland'] == 1
            h = np.full(len(lat), (record['configuration']['start_hour_angle_deg']+lon[0]+
                                  360*day/SYNODIC_MONTH_DAYS+180) % 360-180)
        else:
            h = hour_angle(geo, day*86400, SYNODIC_MONTH_DAYS*86400)
        row = column_summary(d, rho, centres, edges)
        row.update(hour_angle_deg=h, sun_elevation_deg=np.degrees(np.arcsin(np.cos(np.radians(lat))*np.cos(np.radians(h)))),
                   rain_mm_h=d['prate']*3600, temperature_2m_k=d['t2'])
        for key, value in row.items():
            rows.setdefault(key, []).append(np.asarray(value, dtype='f4'))
        times.append(day); snapshots.append(int(n))
        path = case/f'cm1out_t{n:06d}_s.dat'
        raw.append(dict(snapshot=int(n), time_day=day, file=path.name, sha256=sha256(path)))
        if count % 40 == 0:
            print(name, count+1, '/', len(numbers), flush=True)
    metadata = dict(schema=SCHEMA, case=name, gcm_run=record['reference']['run'],
                    producer={str(Path(__file__).name):sha256(__file__),
                              'cloud_columns.py':sha256(Path(__file__).with_name('cloud_columns.py')),
                              'ring_analysis.py':sha256(Path(__file__).with_name('ring_analysis.py'))},
                    case_sha256=sha256(case/'case.json'), source_hashes=sources,
                    span_days=[times[0],times[-1]], snapshots=len(times), grid=record['grid'],
                    path=record.get('path'), site=record.get('site'), terrain=record.get('terrain'),
                    constants=constants, droplet_number_cm3=record['configuration']['droplets_cm3'],
                    geometry='Flat great-circle ring at each column coordinate' if ring else
                             'Terrain-following three-dimensional box; solar clock and map marker use the forcing site',
                    reading_rule='Instantaneous three-hour samples within the second complete lunar cycle. '
                                 'Cloud top/base use qc+qi >= 1e-5 kg/kg. Paths retain all saved positive condensate. '
                                 'Ring heights are above the simulated sea-level surface; box heights are ASL. '
                                 'The box layer interfaces are reconstructed from saved mass-level heights.',
                    units=dict(heights='m', water_path='kg m^-2', angles='degrees', time='Earth days'),
                    raw_inputs=raw)
    output.mkdir(parents=True, exist_ok=True)
    path = output/f'{name}_columns.npz'
    np.savez_compressed(path, metadata=json.dumps(metadata), time_day=times, snapshot=snapshots,
                        latitude_deg=lat, longitude_deg=lon, land=land, ground_m=edges[0],
                        **{k:np.stack(v) for k,v in rows.items()})
    print(path, flush=True)
    return dict(case=name, path=str(path), sha256=sha256(path), schema=SCHEMA,
                span_days=metadata['span_days'], snapshots=len(times), columns=len(lat),
                source_case_sha256=metadata['case_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=CASES)
    args = parser.parse_args()
    products = [export_case(name, args.output) for name in args.cases]
    (args.output/'regional_inputs.json').write_text(json.dumps(dict(schema=SCHEMA, products=products), indent=2)+'\n')


if __name__ == '__main__':
    main()
