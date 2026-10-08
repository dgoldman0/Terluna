"""Choose a conditional wave age from the native upwind shoreline distance."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from research.studies.sea_appearance.nobili import ROOT, digest
from illumination.water_surface.fetch import inverse_wave_age
from illumination.water_surface.full_spectrum import full_curvature
from illumination.water_surface.short_waves import capillary_scales
from shared.constants import MOON_SURFACE_GRAVITY as gravity, MOON_RADIUS


def build(scene_path, terrain_manifest=None):
    product = json.loads(scene_path.read_text())
    wind = product['wind']
    upwind = (270-wind['toward_cartesian_deg'])%360
    profile = product['terrain']['profile']
    az = np.asarray(profile['azimuth_deg'])
    index = int(np.argmin(abs((az-upwind+180)%360-180)))
    fetch = profile['first_land_m'][index]
    sampling = 50.
    fetch_source = 'Saved native terrain horizon profile'
    if fetch is None:
        if terrain_manifest is None:
            raise ValueError('No upwind coast in the admitted 40-km terrain window; fetch is only bounded below')
        if abs(product['site_lon_lat_deg'][1]) > .001 or min(abs(upwind-90),abs(upwind-270)) > .001:
            raise ValueError('The extended fetch admission currently requires an equatorial east/west wind')
        from research.studies.sea_appearance.cloud_candidates import native_heights
        sampling = 100.
        distances = np.arange(sampling, 500000.+sampling/2, sampling)
        sign = 1 if upwind < 180 else -1
        lons = product['site_lon_lat_deg'][0]+sign*np.degrees(distances/MOON_RADIUS)
        if np.any((lons <= 3)|(lons >= 177)):
            raise ValueError('Extended fetch leaves the admitted equatorial terrain strip')
        rows = [dict(observer_lon=float(l),distance_m=float(d)) for l,d in zip(lons,distances)]
        heights = native_heights(rows,json.loads(terrain_manifest.read_text()))
        land = [r for r in heights if r['native_height_above_water_m'] >= 0]
        if not land:
            raise ValueError('No native upwind shore within 500 km; do not claim an exact fetch')
        fetch = land[0]['distance_m']
        fetch_source = dict(terrain_manifest_sha256=digest(terrain_manifest),
            method='100-m equatorial samples of admitted native LOLA rows above the GRAIL geoid and atlas water level',
            first_land_lon_deg=land[0]['observer_lon'],first_land_height_m=land[0]['native_height_above_water_m'])
    u, ust = wind['speed_10m_m_s'], wind['friction_velocity_m_s']
    age = inverse_wave_age(fetch, u, gravity)
    kp = gravity*age**2/u**2
    km, _ = capillary_scales(gravity)
    k = np.geomspace(kp/30, 15*km, 12000)
    long, short = full_curvature(k, u, ust, age, gravity)
    variance = float(np.trapezoid((long+short)/k**2, np.log(k)))
    return dict(schema='terluna.research.coastal-waves/1', scene_sha256=digest(scene_path),
        selected_inverse_wave_age=age, upwind_azimuth_deg=upwind,
        sampled_azimuth_deg=float(az[index]), fetch_m=fetch, fetch_sampling_m=sampling, fetch_source=fetch_source,
        significant_height_m=4*np.sqrt(variance), peak_wavelength_m=2*np.pi/kp,
        model='Elfouhaily et al. (1997), eq. 37 with X=g*fetch/U10^2; full ECKV spectrum with corrected eq. 41.',
        source_url='https://archimer.ifremer.fr/doc/00091/20226/17877.pdf',
        evidence='Native shoreline fetch and saved experimental wind, evaluated with an Earth empirical development law at lunar gravity.',
        limits=['Steady wind over a straight unobstructed fetch is conditional; no duration/history, inherited swell, refraction or directional-fetch integration.',
                'The flat-ring wind does not resolve local terrain circulation. No validated lunar sea-state prediction is claimed.'],
        producer={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__), ROOT/'illumination/water_surface/fetch.py', ROOT/'illumination/water_surface/full_spectrum.py',ROOT/'research/studies/sea_appearance/cloud_candidates.py',ROOT/'geography/topography.py',ROOT/'geography/results/atlas.json',ROOT/'shared/constants.json']})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scene', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--terrain-manifest', type=Path)
    a = p.parse_args()
    record = build(a.scene, a.terrain_manifest)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
