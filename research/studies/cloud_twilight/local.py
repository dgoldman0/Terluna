"""Dusk cloud composition and liquid-extinction screens from local CM1 output.

OPENBLAS_NUM_THREADS=1 python -m research.studies.cloud_twilight.local
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from climate.crm.cloud_columns import SCHEMA as INPUT_SCHEMA, SPECIES, export, sha256, water_path
from illumination.cloud_light.condensate import liquid_optical_depth
from illumination.cloud_light.model import photopic_weights, solved_moon
from .run import HERE, ROOT, path_row

SCHEMA = 'terluna.research.cloud-twilight-local/1'
RUNS = ROOT / 'research/runs/optical_comfort'
SOURCE = ROOT / 'climate/crm/runs/ring_equator'
INPUT = RUNS / 'dusk_columns.npz'
THRESHOLD = 1e-5
BANDS = ((0., 20000.), (20000., 40000.), (40000., 60000.), (60000., 150000.))


def statistics(values):
    values = np.asarray(values)
    if values.size == 0:
        return dict(count=0, mean=None, p50=None, p90=None, maximum=None)
    return dict(count=int(values.size), mean=float(values.mean()), p50=float(np.quantile(values, .5)),
                p90=float(np.quantile(values, .9)), maximum=float(values.max()))


def build():
    export(SOURCE, INPUT)
    with np.load(INPUT, allow_pickle=False) as a:
        metadata = json.loads(str(a['metadata']))
        if metadata['schema'] != INPUT_SCHEMA:
            raise ValueError('Unsupported cloud column schema')
        edges = a['height_edges_m']
        height = (edges[1:] + edges[:-1]) / 2
        rho = a['density_kg_m3'].astype(float)
        species = {key: a[key].astype(float) for key in SPECIES}
        hour, time, snapshot, column = (a[key].copy() for key in ('hour_angle_deg', 'time_day', 'snapshot', 'column_index'))
        cloud = species['qc'] + species['qi'] >= THRESHOLD
        any_cloud = cloud.any(axis=1)
        top = np.max(np.where(cloud, height, 0), axis=1)
        paths = {key: water_path(rho, values, edges) for key, values in species.items()}
        bands, band_paths = [], []
        for lo, hi in BANDS:
            values = {key: water_path(rho, q, edges, lo, hi) for key, q in species.items()}
            band_paths.append(values)
            present = cloud[:, (height >= lo) & (height < hi)].any(axis=1)
            bands.append(dict(height_km=[lo / 1000, hi / 1000], cloud_column_fraction=float(present.mean()),
                              species_water_path_kg_m2={key: statistics(v) for key, v in values.items()},
                              mean_frozen_mass_fraction=float(sum(values[k].sum() for k in ('qi', 'qs', 'qg')) /
                                                              max(sum(v.sum() for v in values.values()), 1e-300))))
        liquid = []
        for radius in (5., 10., 20.):
            tau = liquid_optical_depth(paths['qc'], radius * 1e-6, metadata['constants']['liquid_density_kg_m3'])
            liquid.append(dict(prescribed_effective_radius_um=radius, optical_depth_all_dusk_columns=statistics(tau),
                               column_fraction_tau_ge_1=float(np.mean(tau >= 1)),
                               column_fraction_tau_ge_10=float(np.mean(tau >= 10)),
                               optical_depth_cloudy_columns=statistics(tau[any_cloud])))
        # Select by stated mass metrics, before evaluating illumination.
        picks = [(int(np.argmax(paths['qc'])), 'maximum full-column cloud-liquid path'),
                 (int(np.argmax(paths['qi'])), 'maximum full-column cloud-ice path')]
        for (lo, hi), values in zip(BANDS[1:], band_paths[1:]):
            if values['qc'].max() + values['qi'].max() > 0:
                picks.append((int(np.argmax(values['qc'] + values['qi'])), f'maximum cloud path at {lo/1000:g}–{hi/1000:g} km'))
        model, _ = solved_moon()
        weights = photopic_weights(True)
        examples = []
        for i, selection in picks:
            if not any_cloud[i]:
                continue
            examples.append(dict(selection=selection, snapshot=int(snapshot[i]), time_day=float(time[i]),
                                 column_index=int(column[i]), local_hour_angle_deg=float(hour[i]),
                                 cloud_top_km=float(top[i] / 1000),
                                 species_path_kg_m2={key: float(v[i]) for key, v in paths.items()},
                                 cloud_top_molecular_illumination=path_row(model, weights, top[i] / 1000, 0., float(hour[i] - 90))))
        local_time = []
        for lo in (90., 95., 100., 105.):
            mask = (hour >= lo) & (hour < lo + 5)
            local_time.append(dict(hour_angle_deg=[lo, lo + 5], columns=int(mask.sum()),
                                   cloud_fraction=float(any_cloud[mask].mean()),
                                   cloud_top_km_cloudy=statistics(top[mask & any_cloud] / 1000),
                                   species_path_kg_m2={key: statistics(v[mask]) for key, v in paths.items()}))
        closure = max(float(np.max(abs(sum(b[key] for b in band_paths) - paths[key]))) for key in SPECIES)
    manifest = RUNS / 'dusk_source_manifest.json'
    manifest.write_text(json.dumps(metadata, indent=2, allow_nan=False) + '\n')
    source_metadata = {key: value for key, value in metadata.items() if key != 'raw_inputs'}
    files = [Path(__file__), ROOT / 'climate/crm/cloud_columns.py', ROOT / 'illumination/cloud_light/condensate.py',
             ROOT / 'illumination/cloud_light/model.py', HERE / 'run.py', HERE / 'local_sources.json']
    return dict(schema=SCHEMA, producer=dict(lane='research', runner=str(Path(__file__).relative_to(ROOT)),
                files={str(p.relative_to(ROOT)): sha256(p) for p in files},
                input=dict(file=str(INPUT.relative_to(ROOT)), sha256=sha256(INPUT), schema=INPUT_SCHEMA),
                raw_input_manifest=dict(file=str(manifest.relative_to(ROOT)), sha256=sha256(manifest))),
                evidence='Post-spin-up local CM1 dusk columns with separate cloud liquid, rain, cloud ice, snow and graupel. '
                         'Liquid optical depth uses prescribed effective radii and geometric-optics extinction. '
                         'Molecular illumination at selected cloud tops excludes intervening cloud extinction.',
                reading_rule='Statistics sample the equatorial ring at local hour angle [90,110), over the saved '
                             'second-cycle snapshots. Fractions weight each sampled 6-km column equally; they describe '
                             'this ring experiment. Clear columns contribute zero condensate. Cloudy means qc+qi >= '
                             '1e-5 kg/kg in any layer. Height bands split partial layers by overlap. Effective radii '
                             '5, 10 and 20 micrometres are sensitivity choices; the raw fixed droplet count is retained '
                             'for later size-distribution work. Selected mass maxima illustrate mechanisms. '
                             'Appearance, self-shadowing, ice scattering and cloud lifetimes require further transport.',
                units=dict(water_path='kg m^-2', optical_depth='dimensionless', height='km',
                           occurrence='fraction of sampled columns'),
                source=source_metadata, cloud_column_fraction=float(any_cloud.mean()),
                cloud_top_km_cloudy=statistics(top[any_cloud] / 1000),
                species_path_kg_m2={key: statistics(v) for key, v in paths.items()},
                height_bands=bands, local_time_bins=local_time, liquid_extinction_scenarios=liquid,
                illumination_examples=examples,
                illumination_source='Same current solved molecular column and binned ASTM solar source as cloud_twilight.json; '
                                    'these incident cloud-top values retain the earlier spectral approximation.',
                checks=dict(**metadata['checks'], maximum_height_band_mass_closure_kg_m2=closure))


def main():
    product = build()
    path = HERE / 'results/local_clouds.json'
    path.write_text(json.dumps(product, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(source=product['source']['sampled_span_days'],
                          cloud_fraction=product['cloud_column_fraction'], checks=product['checks']), indent=2))


if __name__ == '__main__':
    main()
