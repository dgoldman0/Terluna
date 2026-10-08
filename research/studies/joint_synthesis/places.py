"""Places on the Open Moon through a lunar day and night, with what each line of work says about them.

    python -m research.studies.joint_synthesis.places
    # -> results/places.json

Seven places the branches already study: the summit port, the four coasts of the sea-appearance study, the highland
box of the climate work and a polar plain. For each, the products give the night's own light and the fleet's
(night_sky.json, illumination/fleet_light), the lightning (atmosphere/electricity occurrence), the haze of the
aerosol study's nearest kind of region (illumination/aerosol/haze.json), the nearest sea's monthly tide
(geography/results/tides.json) and the distance to the regional magnets' nearest cable at the turn of their pattern
that keeps the most heritage and science sites clear (magnets_map.json). The kind of region is a reading of each
place, given beside it.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-places/1'
R = HERE / 'results'
NIGHT, MAGNETS = R / 'night_sky.json', R / 'magnets_map.json'
FLEET = ROOT / 'illumination' / 'fleet_light' / 'results' / 'night_light.json'
OCCURRENCE = ROOT / 'atmosphere' / 'electricity' / 'results' / 'occurrence_A28_dim5_moon.json'
HAZE = ROOT / 'illumination' / 'aerosol' / 'results' / 'haze.json'
TIDES = ROOT / 'geography' / 'results' / 'tides.json'
ATLAS = ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz'
MAGNET_SOURCE = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'magnetic_architecture.json'
OUT = R / 'places.json'
PLACES = {                                    # lat, east lon, the occurrence product's site name, the haze region read
    'summit port': (5.4125, 201.3665, 'summit port (5.4 N, 158.6 W)', 'wet_land'),
    'Oceanus Procellarum by Russell': (32.9, 284.6, 'Oceanus Procellarum by Russell (sea appearance)', 'seas'),
    'eastern Smythii headland': (-1.5, 90.0, 'eastern Smythii headland (sea appearance)', 'seas'),
    'southern Mare Nubium': (-25.0, 345.0, 'southern Mare Nubium (sea appearance)', 'seas'),
    'South Pole-Aitken coast by Mare Ingenii': (-34.0, 163.0, 'South Pole-Aitken coast by Mare Ingenii (sea appearance)', 'seas'),
    'highland box': (-44.7, 246.1, 'highland box (climate/crm box_highland)', 'fog_desert'),
    'polar plain': (75.0, 0.0, 'polar plain (75 N, 0 E)', 'polar_dry_land'),
}
EVIDENCE = 'A compilation of the joint products at seven places; each value carries its source product\'s boundaries.'
READING_RULE = ('Light in lux on level ground at local midnight; flashes per km2 a year (central, with the bracket of '
                'scalings and sea factors); visibility in km by day in the central haze case; tides as the nearest '
                'sea\'s median monthly range; cable distances in km.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def main(argv=None) -> int:
    from research.studies.joint_synthesis.magnets_map import loop_points
    night, fleet = json.loads(NIGHT.read_text()), json.loads(FLEET.read_text())
    occ, haze, tides = json.loads(OCCURRENCE.read_text()), json.loads(HAZE.read_text()), json.loads(TIDES.read_text())
    mag, mag_src = json.loads(MAGNETS.read_text()), json.loads(MAGNET_SOURCE.read_text())
    atlas = np.load(ATLAS)
    water, lat_axis, lon_axis = atlas['water_label'], atlas['lat_deg'], atlas['lon_deg']
    LA, LO = np.meshgrid(lat_axis, lon_axis, indexing='ij')
    beta, eq = np.array(fleet['beta_deg']), np.array(fleet['diffuse_lux'])[0]
    roll80 = fleet['night_roll']['80']['diffuse_lux_at_equator'][0] / fleet['night_roll']['0']['diffuse_lux_at_equator'][0]
    out = {}
    for name, (lat, lon, site, region) in PLACES.items():
        side = 'nearside' if (lon % 360.0 < 90.0 or lon % 360.0 > 270.0) else 'far side'
        latitudes = sorted({float(k.split(', ')[1].split()[0]) for k in night['places']})
        nearest_lat = min(latitudes, key=lambda b: abs(b - abs(lat)))
        row = night['places'][f'{side}, {nearest_lat:g} deg']['diffuse_rho_0.001']
        glow = float(np.interp(lat, beta, eq))
        cases = [c['sites'][site]['flash_density'] for by in occ['cases'].values() for c in by.values()]
        central = occ['cases']['linear']['1']['sites'][site]['flash_density']
        # nearest sea with a tide record
        d = np.degrees(np.arccos(np.clip(unit(LA, LO) @ unit(lat, lon), -1, 1)))
        best = None
        for label, sea in tides['seas'].items():
            mask = water == int(label)
            if mask.any():
                dist = float(d[mask].min()) * np.pi / 180 * MOON_RADIUS / 1e3
                if best is None or dist < best[0]:
                    best = (dist, sea['names'][0], sea['monthly_range_median_m']['area_median'])
        cables = {}
        for design, info in mag['designs'].items():
            loops = next(r for r in mag_src['regional'] if r['name'] == design)['loops']
            dmin = min(float(np.degrees(np.arccos(np.clip(unit(*loop_points(l, info['fewest_sites_and_targets'])[:2]) @ unit(lat, lon), -1, 1))).min())
                       for l in loops) * np.pi / 180 * MOON_RADIUS / 1e3
            cables[design] = dict(turn_deg=info['fewest_sites_and_targets'], nearest_cable_km=round(dmin, 0),
                                  inside_0_5_mT_line=dmin < info['corridor_km'])
        h = haze['regions'][region]['cases']['central']['day']
        out[name] = dict(lat=lat, lon=lon, side=side, region_read=region,
                         night=dict(natural_midnight_lux=row['natural_midnight_lux'], fleet_midnight_lux_rho_0_1pct=round(glow, 2),
                                    fleet_midnight_lux_rolled_80=round(glow * roll80, 3),
                                    hours_below_0_1_lux_natural=row['hours_below_dark']['natural']),
                         lightning=dict(central=central, bracket=[min(cases), max(cases)]),
                         haze=dict(visibility_km=h['visibility_km'], aod=h['aod']),
                         tide=dict(nearest_sea=best[1], distance_km=round(best[0], 0), monthly_range_m=round(best[2], 2)) if best else None,
                         magnets=cables)
    product = dict(schema=SCHEMA, producer=dict(study='joint_synthesis', files={'research/studies/joint_synthesis/places.py': digest(__file__)},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (NIGHT, MAGNETS, FLEET, OCCURRENCE, HAZE, TIDES, ATLAS)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, places=out)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for name, p in out.items():
        print(f"{name:42s} {p['side']:8s} night {p['night']['natural_midnight_lux']:.3g} lux (fleet {p['night']['fleet_midnight_lux_rho_0_1pct']}, rolled {p['night']['fleet_midnight_lux_rolled_80']}) | "
              f"flashes {p['lightning']['central']:.4f} [{p['lightning']['bracket'][0]:.4f}-{p['lightning']['bracket'][1]:.4f}] | vis {p['haze']['visibility_km']} km | "
              f"tide {p['tide']['nearest_sea'] if p['tide'] else '-'} {p['tide']['monthly_range_m'] if p['tide'] else ''} m at {p['tide']['distance_km'] if p['tide'] else ''} km | "
              f"cables {[ (c['nearest_cable_km'], c['inside_0_5_mT_line']) for c in p['magnets'].values()]}")
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
