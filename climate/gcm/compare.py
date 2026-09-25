"""Compare ExoPlaSim runs over chosen years: the numbers behind the README's run tables.

    climate/gcm/.venv/bin/python -m climate.gcm.compare A:30-39 A_earth_path_clouds:10-19
    python -m climate.gcm.compare --settle A_earth_path_clouds:1-14 --baseline A:30-39

Each argument is a run folder under climate/gcm/runs and an inclusive range of model years. Values are
means over those years of the 3-day output means, weighted by area with the T21 Gaussian weights.
Ranges and extremes are of 3-day means at a place, not of hourly values; a range is each year's
warmest minus coldest 3-day mean, averaged over the years. Sunlit means have the Sun's zenith-angle
cosine above 0.3 over the 3 days, dark ones below 0.01. The cloud effect needs the
clear-sky fluxes, which are snapshots at each output time (see exoplasim_run.CLEAR_SKY), and is left
out for runs without them.

--settle estimates where a branch that has not finished warming or cooling would settle. It fits a
line to each year's top-of-atmosphere imbalance against its global surface temperature (a Gregory
plot) and finds where the line reaches the baseline run's settled imbalance: PlaSim does not close
its energy budget exactly, so a settled run keeps a small imbalance of its own. It reads only the
runs' progress.json summaries.
"""
import sys
import numpy as np

FIELDS = ('ts', 'tas', 'ntr', 'rst', 'rsut', 'rlut', 'clt', 'pr', 'prw', 'sic', 'lsm', 'czen')
CLEAR = ('rstcs', 'rltcs')
ROWS = [('surface_k', 'Surface temperature (K)', 2), ('sea_surface_k', 'Sea surface (K)', 1),
        ('land_air_k', 'Air over land (K)', 1), ('equator_air_k', 'Air within 12° of the equator (K)', 1),
        ('polar_air_k', 'Air poleward of 70° (K)', 1), ('toa_net_w_m2', 'Top-of-atmosphere balance (W/m²)', 2),
        ('planetary_albedo', 'Planetary albedo', 3), ('cloud_effect_sw_w_m2', 'Cloud effect, sunlight (W/m²)', 1),
        ('cloud_effect_lw_w_m2', 'Cloud effect, infrared (W/m²)', 1), ('cloud_effect_net_w_m2', 'Cloud effect, net (W/m²)', 1),
        ('cloud_cover', 'Cloud cover', 3), ('precipitation_mm_day', 'Precipitation (mm/day)', 2),
        ('water_vapour_kg_m2', 'Water vapour (kg/m²)', 0),
        ('equator_land_day_k', 'Equatorial land air, sunlit 3-day means (K)', 1),
        ('equator_land_night_k', 'Equatorial land air, dark 3-day means (K)', 1),
        ('equator_land_range_k', 'Equatorial land, range of 3-day means at a place in a year (K)', 1),
        ('lat60_land_range_k', 'Land at 55–65°, the same range (K)', 1),
        ('coldest_air_k', 'Coldest 3-day mean anywhere (K)', 1), ('warmest_air_k', 'Warmest 3-day mean anywhere (K)', 1),
        ('sea_ice_share', 'Sea ice (share of the sea)', 3)]


def area_weights(lat, nlon):
    """Gaussian quadrature weights of the T21 rows, spread over the longitudes, summing to one."""
    _, weights = np.polynomial.legendre.leggauss(len(lat))
    w = weights[np.argsort(np.argsort(np.sin(np.radians(lat))))][:, None] * np.ones(nlon)
    return w / w.sum()


def summarise(years, lat):
    """Headline numbers from a list of years, each a dict of (time, lat, lon) fields."""
    w = area_weights(lat, years[0]['ts'].shape[-1])
    mean = lambda x, m=None: float((x * w * (1.0 if m is None else m)).sum() / (w * (1.0 if m is None else m)).sum())
    avg = lambda k: np.mean([y[k].mean(axis=0) for y in years], axis=0)        # mean map over the years
    land = years[0]['lsm'][0] > 0.5
    equator = (np.abs(lat) < 12.0)[:, None] * np.ones_like(land, dtype=float)
    pole = (np.abs(lat) > 70.0)[:, None] * np.ones_like(land, dtype=float)
    swing = np.mean([y['tas'].max(axis=0) - y['tas'].min(axis=0) for y in years], axis=0)
    lat60 = ((np.abs(lat) > 55.0) & (np.abs(lat) < 65.0))[:, None] * np.ones_like(land, dtype=float)
    tas, sun = np.concatenate([y['tas'] for y in years]), np.concatenate([y['czen'] for y in years])
    eq_land = (equator * land)[None] > 0
    out = dict(surface_k=mean(avg('ts')), sea_surface_k=mean(avg('ts'), ~land), land_air_k=mean(avg('tas'), land),
               equator_air_k=mean(avg('tas'), equator), polar_air_k=mean(avg('tas'), pole),
               toa_net_w_m2=mean(avg('ntr')),
               planetary_albedo=-mean(avg('rsut')) / mean(avg('rst') - avg('rsut')),
               cloud_cover=mean(avg('clt')), precipitation_mm_day=mean(avg('pr')) * 86400e3,
               water_vapour_kg_m2=mean(avg('prw')), equator_land_range_k=mean(swing, equator * land), lat60_land_range_k=mean(swing, lat60 * land),
               equator_land_day_k=float(tas[(sun > 0.3) & eq_land].mean()),
               equator_land_night_k=float(tas[(sun < 0.01) & eq_land].mean()),
               coldest_air_k=float(min(y['tas'].min() for y in years)),
               warmest_air_k=float(max(y['tas'].max() for y in years)),
               sea_ice_share=mean(avg('sic'), ~land))
    if all(k in years[0] for k in CLEAR):
        sw, lw = mean(avg('rst') - avg('rstcs')), mean(avg('rlut') - avg('rltcs'))
        out.update(cloud_effect_sw_w_m2=sw, cloud_effect_lw_w_m2=lw, cloud_effect_net_w_m2=sw + lw)
    return out


def settle(temperatures, imbalances, baseline):
    """Where a line fitted to yearly imbalance against yearly temperature reaches the baseline
    imbalance, and the fitted feedback (W/m2 per K; positive when warming raises the net loss)."""
    slope, intercept = np.polyfit(np.asarray(temperatures, dtype=float), np.asarray(imbalances, dtype=float), 1)
    return (baseline - intercept) / slope, -slope


def _years(item):
    import json
    from climate.gcm.exoplasim_run import RUNS
    folder, span = item.split(':')
    first, last = (int(v) for v in span.split('-'))
    years = [y for y in json.loads((RUNS / folder / 'progress.json').read_text())['years'] if first <= y['year'] <= last]
    if len(years) != last - first + 1:
        raise SystemExit(f'{folder} has {len(years)} of the years {first}-{last}')
    return years


def load(folder, first, last):
    import netCDF4
    from climate.gcm.exoplasim_run import RUNS
    years = []
    for year in range(first, last + 1):
        with netCDF4.Dataset(RUNS / folder / 'model' / f'MOST.{year:05d}.nc') as d:
            years.append({k: np.asarray(d[k][:], dtype=float) for k in FIELDS + CLEAR if k in d.variables})
            lat = np.asarray(d['lat'][:], dtype=float)
    return years, lat


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    if args and args[0] == '--settle':
        years, base = _years(args[1]), _years(args[3])
        baseline = float(np.mean([y['toa_net_w_m2'] for y in base]))
        t, feedback = settle([y['surface_k'] for y in years], [y['toa_net_w_m2'] for y in years], baseline)
        print(f'{args[1]}: settles near {t:.1f} K, where its imbalance reaches {args[3]}\'s {baseline:+.2f} W/m2; '
              f'feedback {feedback:.2f} W/m2 per K; last year {years[-1]["surface_k"]:.2f} K at {years[-1]["toa_net_w_m2"]:+.2f} W/m2')
        return 0
    runs = {}
    for item in args:
        folder, span = item.split(':')
        first, last = (int(v) for v in span.split('-'))
        runs[f'{folder}, years {first}–{last}'] = summarise(*load(folder, first, last))
    print('| | ' + ' | '.join(runs) + ' |')
    print('|---|' + '---|' * len(runs))
    for key, label, digits in ROWS:
        cells = [f'{r[key]:.{digits}f}' if key in r else '–' for r in runs.values()]
        print(f'| {label} | ' + ' | '.join(cells) + ' |')
    return 0


if __name__ == '__main__':
    sys.exit(main())
