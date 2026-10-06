"""Screen the saved equatorial cloud experiment for two kinds of sea view.

Scores select geometry for optical follow-up, not a promise of attractive clouds.
The regional experiment and water mask do not resolve local coastal weather.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from climate.crm.ring_analysis import grads_layout, case_geometry, hour_angle
from climate.crm.cloud_columns import dry_density, read_cm1_constants
from illumination.cloud_light.microphysics import extinction
from geography.lunar_ephemeris import earth_and_sun
from geography.topography import Grid, geoid
from research.studies.sea_appearance.lighting import geometry
from research.studies.sea_appearance.nobili import ROOT, checked, digest
from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS


def matching_date(target_solar_lon, approximate_jd):
    def residual(jd):
        s = earth_and_sun(np.array([jd]))['sun'][0]
        return (np.degrees(np.arctan2(s[1], s[0]))-target_solar_lon+180)%360-180
    return float(brentq(residual, approximate_jd-.4, approximate_jd+.4))


def native_heights(rows, manifest):
    """Equatorial native-pixel interpolation; do not extrapolate admitted rows."""
    valid = [r for r in rows if 3 < r['observer_lon'] < 177]
    if not valid:
        return []
    lon = np.array([r['observer_lon'] for r in valid])
    datum = geoid(Grid(np.array([0.]), lon, 256))[0]
    level = json.loads((ROOT/'geography/results/atlas.json').read_text())['sea_level_m']
    heights = []
    for side in ['north', 'south']:
        m = manifest['terrain'][side]
        source = checked(ROOT/m['local_image'], m['sha256'])
        a = np.memmap(source, mode='r', dtype='<i2', shape=(m['row_count'], m['samples']))
        row = a[-1] if side == 'north' else a[0]
        x = lon*manifest['terrain']['pixels_per_degree']-.5
        i = np.floor(x).astype(int)
        heights.append(((1-(x-i))*row[i]+(x-i)*row[i+1])*m['scale_m'])
    h = (heights[0]+heights[1])/2-datum-level
    return [dict(r, native_height_above_water_m=float(z), offshore_admitted=bool(z < 0))
            for r, z in zip(valid, h)]


def scan(case, mode, first=238, last=473, max_overhead_tau=None, native_manifest=None, solar_range=(-15., -2.)):
    geo = case_geometry(case)
    nx, layout = grads_layout(case/'cm1out_s.ctl')
    lon = np.array(geo['record']['path']['lon'])
    offsets = {}
    offset = 0
    for name, nz in layout:
        offsets[name] = (offset, nz)
        offset += nx*nz
    constants, _ = read_cm1_constants(geo['record'])
    native = None
    if native_manifest is not None:
        rows = [dict(observer_lon=float(v), observer_column=i) for i, v in enumerate(lon)]
        native = {r['observer_column']:r['native_height_above_water_m'] for r in native_heights(rows, native_manifest) if r['offshore_admitted']}
    candidates, sources = [], {}
    for n in range(first, last+1):
        day = (n-1)*.125
        ha = hour_angle(geo, day*86400, SYNODIC_MONTH_DAYS*86400)
        target = (lon[381]-ha[381])%360
        jd = matching_date(target, 2463769.879982721+(n-367)*.125)
        solar = 90-abs(ha)
        use = (~geo['land']) & ((solar < -45) if mode == 'earth' else ((solar > solar_range[0])&(solar < solar_range[1])))
        if mode == 'earth':
            use &= (lon >= 15)&(lon <= 80)
            if geometry(np.array([jd]), 45., 0)['earth_phase_angle'][0] > 110:
                continue
        if native is not None:
            use &= np.array([i in native for i in range(nx)])
        observers = np.flatnonzero(use)
        if not len(observers):
            continue
        source = case/f'cm1out_t{n:06d}_s.dat'
        sources[source.name] = digest(source)
        raw = np.memmap(source, dtype='<f4', mode='r')
        def field(name):
            start, nz = offsets[name]
            return raw[start:start+nz*nx].reshape(nz, nx)
        tau = None
        if max_overhead_tau is not None:
            rho, _ = dry_density(field('prs'), field('th'), field('qv'), constants)
            fields = {k:field(k) for k in ['qc','qr','qi','qs','qg','nci','ncr','ncs','ncg']}
            parts, _ = extinction(fields, rho, geo['record']['configuration']['droplets_cm3'])
            tau = (sum(parts.values())*np.diff(geo['zw'])[:,None]).sum(axis=0)
            observers = observers[tau[observers] <= max_overhead_tau]
        cloud = field('qc')+field('qi') >= 1e-5
        tops = np.where(cloud, geo['zh'][:, None], -1).max(axis=0)
        for i in observers:
            orbit = geometry(np.array([jd]), lon[i], 0)
            earth_el = orbit['earth_elevation'][0]
            if mode == 'earth' and not 12 <= earth_el <= 58:
                continue
            if mode == 'sun' and np.any(cloud[:, i]):
                continue
            step = np.arange(5, 41)
            cols = (i-step)%nx
            distance = step*geo['dx']
            theta = distance/MOON_RADIUS
            top = tops[cols]
            el = np.degrees(np.arctan2((MOON_RADIUS+top)*np.cos(theta)-(MOON_RADIUS+2),
                                      (MOON_RADIUS+top)*np.sin(theta)))
            if mode == 'earth':
                allowed = (top > 20000)&(top < 90000)&(el > 10)&(el < 55)&(abs(el-earth_el) > 5)&(abs(el-earth_el) < 35)
            else:
                allowed = (top > 30000)&(top < 90000)&(el > 15)&(el < 60)
            for j in np.flatnonzero(allowed):
                col = cols[j]
                k = np.flatnonzero(cloud[:, col])
                thickness = (geo['zw'][k+1]-geo['zw'][k]).sum()
                score = float((top[j]/1000)**.6*(thickness/1000)**.4)
                candidates.append(dict(score=score, snapshot=n, jd_tt=jd, overhead_tau=float(tau[i]) if tau is not None else None, native_height_above_water_m=native.get(int(i)) if native is not None else None,
                    observer_column=int(i), observer_lon=float(lon[i]), cloud_column=int(col), cloud_lon=float(lon[col]),
                    top_km=float(top[j]/1000), cloud_depth_sum_km=float(thickness/1000), range_km=float(distance[j]/1000),
                    cloud_elevation_deg=float(el[j]), earth_elevation_deg=float(earth_el),
                    earth_azimuth_deg=float(orbit['earth_azimuth'][0]), earth_phase_deg=float(orbit['earth_phase_angle'][0]),
                    sun_elevation_deg=float(orbit['sun_elevation'][0])))
    ranked = sorted(candidates, key=lambda r: r['score'], reverse=True)
    selected = []
    for row in ranked:
        if any(abs(row['observer_column']-v['observer_column']) < 5 and abs(row['snapshot']-v['snapshot']) < 3 for v in selected):
            continue
        selected.append(row)
        if len(selected) == 30:
            break
    return dict(mode=mode, candidate_pairs=len(candidates), selected=selected, snapshot_hashes=sources, optical_admission=dict(max_overhead_tau=max_overhead_tau, native_water_only=native is not None, solar_range=list(solar_range)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', type=Path, required=True)
    p.add_argument('--terrain-manifest', type=Path, default=Path(__file__).with_name('nobili_inputs.json'))
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--mode', choices=['earth', 'sun'], required=True)
    p.add_argument('--max-overhead-tau', type=float)
    p.add_argument('--native-water-only', action='store_true')
    p.add_argument('--solar-range', type=float, nargs=2, default=[-15., -2.])
    a = p.parse_args()
    manifest = json.loads(a.terrain_manifest.read_text())
    record = scan(a.case, a.mode, max_overhead_tau=a.max_overhead_tau, native_manifest=manifest if a.native_water_only else None, solar_range=a.solar_range)
    record.update(schema='terluna.research.cloud-candidates/1',
        native_terrain_checks=native_heights(record['selected'], manifest),
        evidence='Geometry screen of one experimental equatorial ring, second lunar cycle. No global climate optimum, 3-D cloud shape or dated weather forecast.',
        reading_rule='Cloud selection uses qc+qi >= 1e-5 kg/kg; tops are cell centres, not scale heights. Optical follow-up includes all five condensates. Clear-overhead selection is not zero optical depth. Native-water admission and finite-source optical transmission are separate checks.',
        producer_sha256=digest(__file__), case_sha256=digest(a.case/'case.json'), terrain_manifest_sha256=digest(a.terrain_manifest))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(dict(mode=a.mode, pairs=record['candidate_pairs'], native_water=sum(r['offshore_admitted'] for r in record['native_terrain_checks']))))


if __name__ == '__main__':
    main()
