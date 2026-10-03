"""Stage 1 of the atmospheric-electricity study: can the CM1 storms build a breakdown field against the air's leakage?

    python -m research.studies.atmospheric_electricity.stage1

writes results/stage1.json beside this file.

A storm column separates charge at the generator current density J, the charging rate integrated up its charging zone
(climate/crm/mixed_phase_analysis.py, by charging law). The field between the separated charges grows as
dE/dt = (J - sigma E)/eps0, so it can reach the lightning-initiation threshold E_th only where J exceeds the current the
air conducts at that field, J_crit = sigma E_th; the time to reach it is tau ln(J / (J - J_crit)) with tau = eps0/sigma.
E_th is the relativistic-runaway threshold, 284 kV/m scaled by the air's number density against 2.69e25 per m3
(Dwyer 2007); sigma is the conductivity inside cloud at the charging zone's height and cloud water
(atmosphere/electricity/conductivity.py). Earth's reference uses the same chain at 6 km in a mid-latitude storm.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

from atmosphere.electricity import conductivity as cd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.research.atmospheric-electricity-stage1/1'
RUNS = ('box_0e', 'box_0e_small', 'box_0e_small_earth_fall', 'ring_equator', 'ring_70_45e', 'ring_70_135e')
E_TH_STANDARD = 284.0e3                       # V/m at n0 (Dwyer 2007)
N0 = 2.69e25                                  # m-3
AIR_MOLAR_MASS = 0.02896
AVOGADRO = 6.02214076e23
EARTH_ZONE_KM = 6.0                           # Earth's main charging zone, about -10 to -20 C at mid-latitudes
EARTH_CLOUD = dict(lwc_g_m3=1.0, droplets_m3=3.0e8)
EARTH_STORM_CURRENT_NA_M2 = (20.0, 1000.0)    # charging rates of 10-300 pC m-3 s-1 over 2-4 km (Mansell et al. 2005;
                                              # Kuhlman et al. 2006), as integrated here


def threshold_field(rho_kg_m3):
    """The runaway-breakdown threshold (V/m) at air density rho."""
    return E_TH_STANDARD * np.asarray(rho_kg_m3, float) * AVOGADRO / AIR_MOLAR_MASS / N0


def zone_mean(bands: list) -> dict:
    """The charging zone's properties weighted by each temperature band's depth."""
    rows = [b for b in bands if b.get('depth_m', 0.0) > 0.0]
    w = np.array([b['depth_m'] for b in rows])
    mean = lambda key: float(np.sum(w * np.array([b[key] for b in rows])) / w.sum())
    return dict(height_km=mean('height_km'), pressure_hpa=mean('pressure_hpa'), air_density_kg_m3=mean('air_density_kg_m3'),
                cloud_water_g_m3=mean('cloud_water_g_m3'))


def column_at(column: dict, pressure_hpa: float, key: str, lwc_g_m3: float | None = None) -> float:
    """A column quantity at a pressure: the clear-air or in-cloud conductivity, the latter interpolated in cloud water."""
    p = np.array(column['pressure_hpa'])[::-1]
    at = lambda values: float(np.exp(np.interp(pressure_hpa, p, np.log(np.array(values)[::-1]))))
    if key == 'cloud':
        lwcs = sorted(column['cloud_s_m'], key=lambda k: float(k.split('_')[0]))
        x = np.log([float(k.split('_')[0]) for k in lwcs])
        y = np.log([at(column['cloud_s_m'][k]) for k in lwcs])
        return float(np.exp(np.interp(np.log(lwc_g_m3), x, y)))
    return at(column['clear_air_s_m'][key])


def exceeding(histogram: dict, threshold_na_m2: float) -> float:
    """The share of storm columns whose current exceeds a threshold, from the log histogram (linear within a bin)."""
    edges = np.array(histogram['edges_log10_na_m2'])
    counts = np.array(histogram['counts'], float)
    if counts.sum() == 0:
        return 0.0
    x = np.log10(threshold_na_m2)
    above = counts[edges[:-1] >= x].sum()
    k = np.searchsorted(edges, x) - 1
    if 0 <= k < len(counts):
        above += counts[k] * (edges[k + 1] - x) / (edges[k + 1] - edges[k])
    return float(above / counts.sum())


def time_to_breakdown_h(j_na_m2: float, j_crit_na_m2: float, tau_s: float):
    if j_na_m2 is None or j_na_m2 <= j_crit_na_m2:
        return None
    return float(tau_s * np.log(j_na_m2 / (j_na_m2 - j_crit_na_m2)) / 3600.0)


def earth_reference() -> dict:
    """Earth's in-cloud conductivity, threshold field and critical current at 6 km in a mid-latitude storm."""
    z = np.array([EARTH_ZONE_KM])
    p, t = cd.us_standard_atmosphere(z)
    mu = cd.mobility_air(cd.MOBILITY_STANDARD[0], p, t)
    radius = (3.0 * EARTH_CLOUD['lwc_g_m3'] * 1e-3 / (4.0 * np.pi * 1000.0 * EARTH_CLOUD['droplets_m3'])) ** (1.0 / 3.0)
    sink = cd.droplet_attachment(t, mu, EARTH_CLOUD['droplets_m3'], radius)
    q, sigma_clear, _ = cd.earth_column(z, 3.0, 400.0)
    _, sigma_cloud, _ = cd.earth_column(z, 3.0, 400.0, sink)
    rho = p / (287.053 * t)
    e_th = threshold_field(rho)
    j_crit = float(sigma_cloud[0] * e_th[0] * 1e9)
    return dict(height_km=EARTH_ZONE_KM, pressure_hpa=float(p[0] / 100.0), ion_pairs_per_cm3_s=float(q[0] / 1e6),
                clear_air_s_m=float(sigma_clear[0]), cloud_s_m=float(sigma_cloud[0]),
                relaxation_in_cloud_min=float(cd.relaxation_time(sigma_cloud[0]) / 60.0),
                threshold_kv_m=float(e_th[0] / 1e3), critical_current_na_m2=j_crit,
                storm_current_na_m2=list(EARTH_STORM_CURRENT_NA_M2),
                time_to_breakdown_min=[float(cd.EPSILON_0 * e_th[0] / (j * 1e-9) / 60.0) for j in EARTH_STORM_CURRENT_NA_M2],
                cloud=EARTH_CLOUD)


def assess(solar: str = 'solar_minimum') -> dict:
    column = json.loads((ROOT / 'atmosphere' / 'electricity' / 'results' / 'conductivity_moon.json').read_text())[solar]
    out = {}
    for run in RUNS:
        mp = json.loads((ROOT / 'climate' / 'results' / 'crm' / f'mixed_phase_{run}.json').read_text())
        zone = zone_mean(mp['zone']['bands'])
        sigma_cloud = column_at(column, zone['pressure_hpa'], 'cloud', zone['cloud_water_g_m3'])
        sigma_clear = column_at(column, zone['pressure_hpa'], '0_per_cm3')
        e_th = float(threshold_field(zone['air_density_kg_m3']))
        j_crit = sigma_cloud * e_th * 1e9
        tau = float(cd.relaxation_time(sigma_cloud))
        laws = {}
        for law, regimes in mp['charging'].items():
            laws[law] = {}
            for regime, c in regimes.items():
                cc = c['column_current_na_m2']
                share = exceeding(c['column_current_histogram'], j_crit)
                laws[law][regime] = dict(
                    share_of_storm_columns_above_critical=share,
                    share_of_all_columns_above_critical=share * mp['zone']['share_of_columns'],
                    current_p90_na_m2=cc.get('p90'), current_max_na_m2=cc.get('max'),
                    hours_to_breakdown_p90=time_to_breakdown_h(cc.get('p90'), j_crit, tau),
                    hours_to_breakdown_max=time_to_breakdown_h(cc.get('max'), j_crit, tau))
        out[run] = dict(zone=zone, cloud_s_m=sigma_cloud, clear_air_s_m=sigma_clear,
                        relaxation_in_cloud_h=tau / 3600.0, threshold_kv_m=e_th / 1e3, critical_current_na_m2=j_crit,
                        zone_share_of_columns=mp['zone']['share_of_columns'], laws=laws)
    return out


def main(argv=None) -> int:
    result = dict(schema=SCHEMA, solar_minimum=assess('solar_minimum'), solar_maximum=assess('solar_maximum'),
                  earth_reference=earth_reference())
    result['producer'] = dict(domain='research', files={f'studies/atmospheric_electricity/{Path(__file__).name}':
                                                        hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'stage1.json').write_text(json.dumps(result, indent=1) + '\n')
    e = result['earth_reference']
    print(f"Earth at {e['height_km']} km: q {e['ion_pairs_per_cm3_s']:.3g}/cm3/s, cloud sigma {e['cloud_s_m']:.2e} S/m, "
          f"E_th {e['threshold_kv_m']:.0f} kV/m, J_crit {e['critical_current_na_m2']:.3g} nA/m2")
    for run, r in result['solar_minimum'].items():
        print(f"{run}: zone {r['zone']['height_km']:.1f} km, cloud sigma {r['cloud_s_m']:.2e} S/m, tau {r['relaxation_in_cloud_h']:.2g} h, "
              f"E_th {r['threshold_kv_m']:.0f} kV/m, J_crit {r['critical_current_na_m2']:.3g} nA/m2")
        for law, regimes in r['laws'].items():
            x = regimes['run']
            print(f"    {law:16s} storm columns above: {x['share_of_storm_columns_above_critical']:.3f}; hours to breakdown "
                  f"p90 {x['hours_to_breakdown_p90']} max {x['hours_to_breakdown_max']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
