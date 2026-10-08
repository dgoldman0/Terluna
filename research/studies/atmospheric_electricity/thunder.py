"""Thunder from the lunar box's flashes, heard at the ground (stage 2 of the atmospheric-electricity study).

    python -m research.studies.atmospheric_electricity.thunder

writes results/thunder.json beside this file.

The flashes are stage 2's main run's (box_0e_elec_corrected, two lunar days under the lunar rules). Each releases the
electrostatic energy its log gives, a share of which becomes sound (Holmes et al. 1971: 0.18 %, with a tenth and ten
times that as bounds), spread evenly over its channel from its lowest to its highest point (at its starting point where
the log gives no channel), as point sources a kilometre apart; a ground strike's channel runs on from its lowest point
in cloud down to the ground, which the log does not draw, and takes its share of the energy there. The spectrum peaks
at
Few's frequency for the energy per metre of a channel as long as its height span plus the square root of its area.
The air is the box's own, averaged over the box and over its three-hourly outputs at the hours that flashed:
temperature, pressure, vapour and wind by height, in the design air's composition. atmosphere/electricity/thunder.py
traces the rays in eight directions and gathers the sound at the ground. Earth's benchmark flash (supercell_elec's
median: 1.2 GJ over 4.75-10.75 km and 27 km2), in cloud and as a ground strike, in the US standard atmosphere checks
the method against the 15-25 km over which thunder is heard on Earth.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from atmosphere.electricity import thunder as th
from climate.crm import elec_analysis as ea
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.research.atmospheric-electricity-thunder/1'
RUNS = ra.RUNS
GCM_PROGRESS = ROOT / 'climate' / 'gcm' / 'runs' / 'A28_dim5_moon' / 'progress.json'
FLASH_SETS = (('box_0e_elec_corrected', 0.0, None),)
EFFICIENCY = (th.HOLMES_EFFICIENCY, 0.1 * th.HOLMES_EFFICIENCY, 10.0 * th.HOLMES_EFFICIENCY)
AZIMUTHS = tuple(range(0, 360, 45))
DZ = 250.0
SOURCE_STEP = 1000.0
N_RAYS = 7201
X_EDGES = np.arange(0.0, 400001.0, 2000.0)
QUIET_DBA = (30.0, 45.0)                    # a quiet countryside night, an ordinary daytime background
EARTH_FLASH = dict(energy_j=1.2e9, z_low_m=4750.0, z_high_m=10750.0, extent_m2=27.0e6)
EARTH_GROUND_FLASH = dict(EARTH_FLASH, z_low_m=0.0)       # the same flash, its channel run down to the ground


def flashes() -> dict:
    """The flashes: energy released (J), channel bottom and top (m), extent (m2), start height (m) and whether they
    strike the ground, whose channel then reaches it."""
    rows = []
    for case, t0, t1 in FLASH_SETS:
        f = ea.read_log(RUNS / case / 'terluna_flashes.txt', ea.FLASH_COLUMNS)
        sel = (f['time_s'] >= t0) & ((f['time_s'] <= t1) if t1 else True)
        for k in np.nonzero(sel)[0]:
            rows.append((case, f['time_s'][k], f['energy_before_j'][k] - f['energy_after_j'][k], f['z_low_m'][k],
                         f['z_high_m'][k], f['extent_m2'][k], f['z_m'][k], f['kind'][k] in (2, 3)))
    cases, t, e, lo, hi, ext, z0, ground = zip(*rows)
    lo, hi, z0, ground = np.array(lo), np.array(hi), np.array(z0), np.array(ground)
    unrecorded = lo < 0.0                                     # a flash whose channel the log does not give (-1)
    lo = np.where(unrecorded, z0, lo)
    hi = np.where(unrecorded, z0, hi)
    lo = np.where(ground, 0.0, lo)                            # a ground strike's channel reaches the ground
    return dict(case=np.array(cases), time_s=np.array(t), energy_j=np.array(e), z_low_m=lo, z_high_m=hi,
                extent_m2=np.array(ext), z_start_m=z0, unrecorded=unrecorded, ground=ground)


def air_composition() -> dict:
    gases = json.loads(GCM_PROGRESS.read_text())['configuration']['gases_bar']
    total = sum(gases.values())
    x = {k: v / total for k, v in gases.items()}
    molar = (x['pN2'] * 28.0134 + x['pO2'] * 31.9988 + x['pAr'] * 39.948 + x['pCO2'] * 44.0095) / 1000.0
    return dict(x_n2=x['pN2'], x_o2=x['pO2'], x_ar=x['pAr'], x_co2=x['pCO2'], molar_mass=molar)


def box_air(f: dict, air: dict) -> tuple:
    """The box's air at the hours that flashed, averaged over the box and those outputs, on edges DZ apart."""
    profiles = []
    used = []
    for case in np.unique(f['case']):
        output_s = json.loads((RUNS / case / 'case.json').read_text())['configuration']['output_s']
        for n in sorted({int(round(t / output_s)) + 1 for t in f['time_s'][f['case'] == case]}):
            if not (RUNS / case / f'cm1out_t{n:06d}_s.dat').exists():
                continue
            d = ra.read_snapshot(RUNS / case, n)
            t = d['th'] * (d['prs'] / 1.0e5) ** (287.04 / 1005.7)
            profiles.append([a.mean(axis=1) for a in (t, d['prs'], d['qv'], d['uinterp'], d['vinterp'])])
            used.append(f'{case}:{n}')
    zw = np.loadtxt(RUNS / f['case'][0] / 'input_grid_z')
    zh = 0.5 * (zw[1:] + zw[:-1])
    mean = np.mean(np.array(profiles), axis=0)
    spread = np.std(np.array(profiles), axis=0)
    z = np.arange(0.0, zw[-1] - DZ + 1.0, DZ)
    t, p, qv, u, v = (np.interp(z, zh, a) for a in mean)
    p = np.exp(np.interp(z, zh, np.log(mean[1])))
    atm = th.Atmosphere(z, t, p, qv, u, v, air['molar_mass'], air['x_o2'], air['x_n2'])
    return atm, used, dict(t_k_std_max=float(spread[0].max()), u_m_s_std_max=float(spread[3].max()))


def earth_air() -> th.Atmosphere:
    """The US standard atmosphere to 80 km, still, with 10 g/kg of vapour at the ground falling off over 2.5 km."""
    layers = ((0.0, -6.5e-3), (11000.0, 0.0), (20000.0, 1.0e-3), (32000.0, 2.8e-3), (47000.0, 0.0),
              (51000.0, -2.8e-3), (71000.0, -2.0e-3))
    z = np.arange(0.0, 80001.0, DZ)
    t = np.full(z.size, 288.15)
    for (zl, lapse), (zn, _) in zip(layers, list(layers[1:]) + [(np.inf, 0.0)]):
        t += lapse * np.clip(np.minimum(z, zn) - zl, 0.0, None)
    g, rd = 9.80665, 287.05
    p = np.empty(z.size); p[0] = 101325.0
    for i in range(1, z.size):
        p[i] = p[i - 1] * np.exp(-g * DZ / (rd * 0.5 * (t[i] + t[i - 1])))
    qv = 0.010 * np.exp(-z / 2500.0)
    return th.Atmosphere(z, t, p, qv, np.zeros(z.size), np.zeros(z.size))


def sources(atm: th.Atmosphere, z_lo: float, z_hi: float) -> np.ndarray:
    k = np.arange(np.round(z_lo / SOURCE_STEP), np.round(z_hi / SOURCE_STEP) + 1) * SOURCE_STEP
    return k[(k > 0) & (k < atm.z[-1])]


def peak_frequency(atm: th.Atmosphere, energy, z_lo, z_hi, extent) -> float:
    zm = 0.5 * (z_lo + z_hi)
    length = max(z_hi - z_lo, SOURCE_STEP) + np.sqrt(extent)
    return float(th.few_peak_hz(energy / length, np.interp(zm, atm.z, atm.p), np.interp(zm, atm.z, atm.c)))


def hear(atm, traces, energy, z_lo, z_hi, extent, efficiency):
    """Ground levels from one flash in one direction, from traces by source height."""
    zs = sources(atm, z_lo, z_hi)
    acoustic = efficiency * energy / zs.size
    return th.ground_levels([traces[float(z)] for z in zs], [acoustic] * zs.size,
                            peak_frequency(atm, energy, z_lo, z_hi, extent), atm, X_EDGES)


def summary(levels: dict, f_peak: float) -> dict:
    below = 0
    return dict(peak_hz=round(f_peak, 1),
                below=dict(first_s=float(levels['first_s'][below]), span_s=float(levels['span_s'][below]),
                           peak_db_by_band={str(b): round(float(v), 1) for b, v in zip(th.BANDS, levels['peak_db'][below])},
                           peak_dba=round(float(levels['peak_dba'][below]), 1)),
                audible_km=round(th.audible_range(levels) / 1e3, 1),
                **{f'above_{int(q)}_dba_km': round(th.audible_range(levels, q) / 1e3, 1) for q in QUIET_DBA},
                by_range=[dict(x_km=float(x / 1e3), peak_dba=round(float(a), 1), peak_db_63=round(float(p[2]), 1),
                               peak_db_31=round(float(p[1]), 1), first_s=None if np.isnan(fs) else round(float(fs), 1),
                               span_s=None if np.isnan(sp) else round(float(sp), 1))
                          for x, a, p, fs, sp in zip(levels['x_m'], levels['peak_dba'], levels['peak_db'],
                                                     levels['first_s'], levels['span_s'])
                          if x / 1e3 in (1.0, 11.0, 25.0, 51.0, 75.0, 101.0, 151.0, 201.0, 251.0, 301.0)])


def finite(o):
    """The result with silence (a level of minus infinity) and missing values as null, for strict JSON."""
    if isinstance(o, dict):
        return {k: finite(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [finite(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def main(argv=None) -> int:
    f = flashes()
    air = air_composition()
    atm, used, spread = box_air(f, air)
    q = lambda a, p: float(np.percentile(a, p))

    def kinds(sel):
        """The typical flash of a kind (its medians) and its large one (the top tenth by energy)."""
        typical = dict(energy_j=q(f['energy_j'][sel], 50), z_low_m=q(f['z_low_m'][sel], 50),
                       z_high_m=q(f['z_high_m'][sel], 50), extent_m2=q(f['extent_m2'][sel], 50))
        big = sel & (f['energy_j'] >= q(f['energy_j'][sel], 90))
        large = dict(energy_j=q(f['energy_j'][sel], 90), z_low_m=q(f['z_low_m'][big], 50),
                     z_high_m=q(f['z_high_m'][big], 50), extent_m2=q(f['extent_m2'][big], 50))
        return typical, large
    typical, large = kinds(~f['ground'])
    named = {'typical': typical, 'large': large}
    if f['ground'].any():
        named['ground_typical'], named['ground_large'] = kinds(f['ground'])
    heights = sorted({float(z) for lo, hi in zip(f['z_low_m'], f['z_high_m']) for z in sources(atm, lo, hi)})
    out = {}
    per_flash = {}
    for az in AZIMUTHS:
        traces = {z: th.trace(atm, z, azimuth_deg=az, n_rays=N_RAYS) for z in heights}
        print(f'azimuth {az}: {len(heights)} source heights traced', flush=True)
        for name, fl in named.items():
            for eff in EFFICIENCY:
                key = f'{name}_eff{eff:g}'
                lv = hear(atm, traces, fl['energy_j'], fl['z_low_m'], fl['z_high_m'], fl['extent_m2'], eff)
                out.setdefault(key, {})[str(az)] = summary(lv, peak_frequency(atm, fl['energy_j'], fl['z_low_m'],
                                                                             fl['z_high_m'], fl['extent_m2']))
        ranges = []
        for k in range(f['energy_j'].size):
            lv = hear(atm, traces, f['energy_j'][k], f['z_low_m'][k], f['z_high_m'][k], f['extent_m2'][k],
                      th.HOLMES_EFFICIENCY)
            ranges.append((th.audible_range(lv), th.audible_range(lv, QUIET_DBA[0]), th.audible_range(lv, QUIET_DBA[1]),
                           float(lv['peak_dba'][0])))
        per_flash[str(az)] = np.array(ranges)
    stack = np.array([per_flash[str(az)] for az in AZIMUTHS])           # azimuth, flash, measure
    def ranges(sel):
        mean = stack[:, sel, :].mean(axis=0)
        out = {name: {p: round(float(np.percentile(mean[:, i], p)) / 1e3, 1) for p in (10, 50, 90)}
               for i, name in enumerate(('audible_km', 'above_30_dba_km', 'above_45_dba_km'))}
        out['below_peak_dba'] = {p: round(float(np.percentile(mean[:, 3], p)), 1) for p in (10, 50, 90)}
        return out
    flash_ranges = ranges(~f['ground'])
    ground_ranges = ranges(f['ground']) if f['ground'].any() else None

    earth = earth_air()
    e_heights = [float(z) for z in sources(earth, EARTH_GROUND_FLASH['z_low_m'], EARTH_GROUND_FLASH['z_high_m'])]
    e_traces = {z: th.trace(earth, z, n_rays=N_RAYS) for z in e_heights}
    earth_runs = {}
    for name, fl in (('in_cloud', EARTH_FLASH), ('ground', EARTH_GROUND_FLASH)):
        earth_runs[name] = {str(eff): summary(hear(earth, e_traces, fl['energy_j'], fl['z_low_m'], fl['z_high_m'],
                                                   fl['extent_m2'], eff),
                                              peak_frequency(earth, fl['energy_j'], fl['z_low_m'], fl['z_high_m'],
                                                             fl['extent_m2']))
                            for eff in EFFICIENCY}
    earth_out = earth_runs['in_cloud']

    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    result = dict(
        schema=SCHEMA,
        producer=dict(study='research/studies/atmospheric_electricity', files={
            'thunder.py': digest(__file__), 'atmosphere/electricity/thunder.py': digest(th.__file__)}),
        evidence=('Geometric acoustics through the lunar box\'s mean air at the hours that flashed (horizontally '
                  'uniform, with its mean wind), ISO 9613-1 absorption, a hard flat ground, and flashes whose sound is '
                  'a fixed share of their electrostatic energy spread evenly along a vertical channel, a ground strike\'s '
                  'run down to the ground. Not a model of '
                  'thunder\'s generation: the acoustic share is uncertain by at least ten times either way, storms\' '
                  'own winds and cold pools near the ground are left out, and geometric acoustics fails near caustics '
                  'and at grazing angles.'),
        reading_rule=('Levels are sound pressure levels at a listener on the ground, in dB re 20 uPa, as the loudest '
                      'one-second average; dBA A-weighted over the audible bands (31.5 Hz and up). audible_km is the '
                      'farthest 2-km range bin where some audible band reaches the threshold of hearing (ISO 226) in '
                      'its loudest second; above_30_dba_km and above_45_dba_km where the A-weighted level reaches a '
                      'quiet countryside night and an ordinary daytime background. first_s is the delay after the flash; '
                      'span_s holds 90 % of the audible sound energy. Directions are azimuths clockwise from north.'),
        inputs=dict(flashes={case: int((f['case'] == case).sum()) for case in np.unique(f['case'])},
                    ground_strikes=int(f['ground'].sum()),
                    flashes_without_a_channel_placed_at_their_start=int(f['unrecorded'].sum()),
                    profile_outputs=used, profile_spread=spread, air=air,
                    acoustic_efficiency=dict(central=th.HOLMES_EFFICIENCY, bounds=EFFICIENCY[1:])),
        typical_flash=typical, large_flash=large, lunar=out, flash_ranges_mean_over_directions=flash_ranges,
        ground_strikes=dict(typical=named.get('ground_typical'), large=named.get('ground_large'),
                            ranges_mean_over_directions=ground_ranges),
        earth_benchmark_flash=dict(flash=EARTH_FLASH, results=earth_out),
        earth_benchmark_ground_strike=dict(flash=EARTH_GROUND_FLASH, results=earth_runs['ground']),
        lapse_k_per_km=dict(lunar_0_30_km=round(float((atm.t[0] - np.interp(30000.0, atm.z, atm.t)) / 30.0), 2)))
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'thunder.json').write_text(json.dumps(finite(result), indent=1, allow_nan=False) + '\n')
    print(json.dumps(dict(typical=out[f'typical_eff{th.HOLMES_EFFICIENCY:g}']['90'] | {'by_range': None},
                          large=out[f'large_eff{th.HOLMES_EFFICIENCY:g}']['90'] | {'by_range': None},
                          flashes=flash_ranges,
                          earth=earth_out[str(th.HOLMES_EFFICIENCY)] | {'by_range': None}), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
