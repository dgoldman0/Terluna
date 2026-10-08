"""The nitrogen the lunar lightning fixes, from the main run's flashes and Earth's measured yields.

    python -m research.studies.atmospheric_electricity.nitrogen

writes results/nitrogen_oxides.json beside this file.

The model's own nitrogen oxides (lightmsz's yield per channel point) count the points its channels take on the grid, so
they run about fifty times low on Earth's benchmark storm and lower still on the lunar 6-km grid. The estimate here
takes the main run's flashes (box_0e_elec_corrected) and three ways Earth's yields are expressed, which part for flashes
as large as these: no published yield covers flashes that release tens to hundreds of gigajoules and move a hundred
coulombs or more (literature note lightning_nox.md on the data drive):
- per flash: Earth's mean, 250 mol of NOx a flash (150-350 likely, 33-664 the full range; [S53], [S54]), the same in
  cloud and to ground ([S54], [S59]);
- per joule of the energy a flash releases: 9-10e16 molecules a joule ([S55], [S56]), 1-50e16 the range ([S53]);
- per metre of channel: Wang et al.'s laboratory yield at the channel's pressure, 0.34e21 + 1.30e16 p[Pa] molecules a
  metre ([S33] as fitted in [S58]), 0.3-13e21 the observed range ([S53], [S54]), on a channel as long as the flash's
  height span plus the square root of its area (a ground strike's from its lowest point in cloud down to the ground).
  The flash log draws no branches, so this is a lower bound.
Each is taken in the design air's oxygen, 0.175 against Earth's 0.209, at 0.84-0.93 of Earth's yield (derived in the
note from equilibrium; no measurement exists), and set against Earth's 5 Tg N a year (3-7 [S54]; 2-8 [S53]) over its
surface, 9.8 kg N per km2 a year, from 44 flashes a second ([S53]).

The product also records the flashes' sizes beside the Earth benchmark storm's and how that storm's flashes change from its
1-km grid to its 2-km one (the electrified CM1's products in climate/results/crm/): coarsening the grid makes its flashes
fewer and larger, so the box's 6-km grid bears on the routes in opposite ways.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from climate.crm import elec_analysis as ea
from climate.crm import ring_analysis as ra
from shared.constants import AVOGADRO, EARTH_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.research.atmospheric-electricity-nitrogen/1'
RUN = 'box_0e_elec_corrected'
COLUMN = ROOT / 'atmosphere' / 'electricity' / 'results' / 'conductivity_moon.json'
CRM_RESULTS = ROOT / 'climate' / 'results' / 'crm'
BENCHMARK = ('supercell_elec', 'supercell_elec_2km')                # Earth's benchmark storm on its 1-km and 2-km grids
N_G_MOL = 14.007
PER_FLASH_MOL = dict(central=250.0, likely=(150.0, 350.0), full=(33.0, 664.0))
PER_JOULE = dict(central=9.5e16, published=(9.0e16, 10.0e16), full=(1.0e16, 50.0e16))     # molecules per joule
PER_METRE = dict(a=0.34e21, b=1.30e16, observed=(0.3e21, 13.0e21))                       # molecules per metre
OXYGEN = dict(central=0.93, range=(0.84, 0.93))
EARTH_TG_N_YR = dict(central=5.0, range=(3.0, 7.0))                  # [S54]; [S53] gives 2-8
EARTH_FLASHES_PER_S = 44.0                                           # [S53]
EARTH_KM2 = 4.0 * np.pi * EARTH_RADIUS ** 2 / 1.0e6
EARTH_KG_N_KM2_YR = dict(central=EARTH_TG_N_YR['central'] * 1.0e9 / EARTH_KM2,
                         range=tuple(v * 1.0e9 / EARTH_KM2 for v in EARTH_TG_N_YR['range']))


def per_metre(p_pa):
    """Wang et al.'s yield (molecules a metre) at pressure p_pa, as Davis et al. fit it."""
    return PER_METRE['a'] + PER_METRE['b'] * np.asarray(p_pa, float)


def channel_metres(z_low_m, z_high_m, area_m2, ground):
    """A flash's channel length: its height span (a ground strike's down to the ground) plus the square root of its
    area, the same proxy the thunder takes."""
    low = np.where(ground, 0.0, z_low_m)
    return np.maximum(np.asarray(z_high_m) - low, 1000.0) + np.sqrt(np.asarray(area_m2))


def estimate(energy_j, z_low_m, z_high_m, area_m2, ground, column_z_m, column_p_pa, area_km2, days):
    """Nitrogen fixed (mol of NOx) by the flashes under each way of counting, and per km2 a year over the box."""
    years = days / 365.25
    length = channel_metres(z_low_m, z_high_m, area_m2, ground)
    low = np.where(ground, 0.0, z_low_m)
    p_mid = np.interp(0.5 * (low + np.asarray(z_high_m)), column_z_m, column_p_pa)
    n = np.asarray(energy_j).size
    routes = {
        'per_flash': dict(central=n * PER_FLASH_MOL['central'],
                          range=[n * PER_FLASH_MOL['full'][0], n * PER_FLASH_MOL['full'][1]]),
        'per_joule': dict(central=float(np.sum(energy_j)) * PER_JOULE['central'] / AVOGADRO,
                          range=[float(np.sum(energy_j)) * v / AVOGADRO for v in PER_JOULE['full']]),
        'per_metre': dict(central=float(np.sum(length * per_metre(p_mid))) / AVOGADRO,
                          range=[float(np.sum(length)) * v / AVOGADRO for v in PER_METRE['observed']]),
    }
    out = {}
    for name, r in routes.items():
        per_km2 = lambda mol: mol / area_km2 / years
        kg_n = lambda mol: per_km2(mol) * N_G_MOL / 1000.0
        out[name] = dict(
            total_mol=r['central'], total_mol_range=r['range'],
            per_flash_mol=r['central'] / n,
            kg_n_per_km2_yr=kg_n(r['central']) * OXYGEN['central'],
            kg_n_per_km2_yr_range=[kg_n(r['range'][0]) * OXYGEN['range'][0], kg_n(r['range'][1]) * OXYGEN['range'][1]],
            share_of_earth=kg_n(r['central']) * OXYGEN['central'] / EARTH_KG_N_KM2_YR['central'])
    out['channel_km'] = dict(median=float(np.median(length)) / 1000.0, total=float(np.sum(length)) / 1000.0)
    rate = n / area_km2 / years
    out['flashes_per_km2_yr'] = dict(box=rate, earth=EARTH_FLASHES_PER_S * 365.25 * 86400.0 / EARTH_KM2)
    out['flashes_per_km2_yr']['box_over_earth'] = rate / out['flashes_per_km2_yr']['earth']
    return out


def flash_sizes(run: str = RUN, benchmark=BENCHMARK) -> dict:
    """The run's flashes beside the Earth benchmark storm's, from the electrified CM1's products: median energy and
    charge of each kind, and how the benchmark's flash count, median energy and total energy change from its 1-km grid
    to its 2-km one."""
    lightning = lambda name: json.loads((CRM_RESULTS / f'elec_{name}.json').read_text())['lightning']
    run_l = lightning(run)
    out = {kind: dict(flashes=run_l[kind]['count'],
                      energy_j={k: run_l[kind]['energy_dissipated_j'][k] for k in ('p10', 'median', 'p90', 'max')},
                      charge_c_median=max(run_l[kind]['positive_c']['median'], run_l[kind]['negative_c']['median']))
           for kind in ('in_cloud', 'negative_to_ground')}
    earth = {}
    for name in benchmark:
        b = lightning(name)['in_cloud']
        earth[name] = dict(flashes=b['count'], median_energy_j=b['energy_dissipated_j']['median'],
                           total_energy_j=b['count'] * b['energy_dissipated_j']['mean'])
    fine, coarse = (earth[name] for name in benchmark)
    earth['coarse_over_fine'] = {k: coarse[k] / fine[k] for k in fine}
    out['earth_benchmark'] = earth
    out['median_energy_over_earth_benchmark'] = {kind: out[kind]['energy_j']['median'] / fine['median_energy_j']
                                                 for kind in ('in_cloud', 'negative_to_ground')}
    return out


def main(argv=None) -> int:
    case = ra.RUNS / RUN
    record = json.loads((case / 'case.json').read_text())
    f = ea.read_log(case / 'terluna_flashes.txt', ea.FLASH_COLUMNS)
    ground = np.isin(f['kind'], (2, 3))
    energy = f['energy_before_j'] - f['energy_after_j']
    z_low = np.where(f['z_low_m'] < 0.0, f['z_m'], f['z_low_m'])
    z_high = np.where(f['z_high_m'] < 0.0, f['z_m'], f['z_high_m'])
    column = json.loads(COLUMN.read_text())['solar_minimum']
    grid = record['grid']
    area_km2 = grid['nx'] * grid['dx_m'] * grid.get('ny', grid['nx']) * grid.get('dy_m', grid['dx_m']) / 1.0e6
    days = record['configuration']['days']
    result = dict(
        schema=SCHEMA,
        producer=dict(study='research/studies/atmospheric_electricity',
                      files={'nitrogen.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}),
        run=RUN, flashes=int(ground.size), ground_strikes=int(ground.sum()), days=days, box_km2=area_km2,
        evidence=('Earth\'s measured lightning nitrogen-oxide yields, per flash, per joule of the energy a flash releases '
                  'and per metre of channel (literature note lightning_nox.md), applied to the main run\'s flashes. '
                  'None is measured for flashes that release tens to hundreds of gigajoules, and the three part by '
                  'more than a hundred times for them; the oxygen factor is a scaling argument. The flashes\' number '
                  'and energies depend on the 6-km grid (flash_sizes). The box is one equatorial lowland site.'),
        reading_rule=('total_mol is the NOx the run\'s flashes make over its two lunar days under each way of counting, '
                      'with the published range; kg_n_per_km2_yr spreads it over the box and a year and takes the '
                      'design air\'s oxygen; share_of_earth sets that against Earth\'s 9.8 kg N per km2 a year. The '
                      'per-metre route counts no branches and is a lower bound. flash_sizes sets the flashes beside '
                      'the Earth benchmark storm\'s, whose 2-km grid against its 1-km one gives fewer flashes, each '
                      'releasing more energy, and more energy in all; the box\'s coarser grid likely raises the per-'
                      'joule count and lowers the per-flash one, by amounts no finer lunar run has measured.'),
        estimate=estimate(energy, z_low, z_high, f['extent_m2'], ground, 1000.0 * np.array(column['z_km']),
                          100.0 * np.array(column['pressure_hpa']), area_km2, days),
        flash_sizes=flash_sizes(),
        yields=dict(per_flash_mol=PER_FLASH_MOL, per_joule_molecules=PER_JOULE, per_metre_molecules=PER_METRE,
                    oxygen_factor=OXYGEN, earth_tg_n_per_yr=EARTH_TG_N_YR, earth_kg_n_per_km2_yr=EARTH_KG_N_KM2_YR,
                    earth_flashes_per_s=EARTH_FLASHES_PER_S))
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'nitrogen_oxides.json').write_text(json.dumps(result, indent=1) + '\n')
    for name in ('per_flash', 'per_joule', 'per_metre'):
        e = result['estimate'][name]
        print(f"{name:10s} {e['per_flash_mol']:9.0f} mol a flash; {e['kg_n_per_km2_yr']:.3g} kg N per km2 a year "
              f"({e['kg_n_per_km2_yr_range'][0]:.2g}-{e['kg_n_per_km2_yr_range'][1]:.2g}), {e['share_of_earth']:.3g} of Earth's")
    r, f = result['estimate']['flashes_per_km2_yr'], result['flash_sizes']
    print(f"flashes {r['box']:.3g} per km2 a year, {r['box_over_earth']:.3g} of Earth's; median energy "
          f"{f['median_energy_over_earth_benchmark']} times the benchmark's; benchmark 2 km over 1 km "
          f"{f['earth_benchmark']['coarse_over_fine']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
