"""The shield's regional magnets on the Moon's map: where their cables would run, for each turn of the pattern.

    python -m research.studies.joint_synthesis.magnets_map
    # -> results/magnets_map.json

The shield branch's regional magnetic architecture (research/studies/solar_shield_array/results/
magnetic_architecture.json) places four current loops on the surface, two about each pole, but leaves their
longitudes free: its frame puts the north pair at 0 and 180 degrees and the south pair at 90 and 270. This study turns
the whole pattern through 0-85 degrees in 5-degree steps and records, for the four 500 km loops and the four 1,000 km
loops, how much cable crosses the seas of the 28% atlas, which heritage sites (research/studies/conservation/sites.json)
and science targets lie within the cable's 0.5 mT line, and how much of each loop's cap is water. The field near a
cable follows the straight-wire law B = mu0 I / (2 pi d), which holds within a fraction of the loop's radius; 0.5 mT
is the field the guidance for implanted medical devices keeps people below. Which turn to choose is a siting decision
for people; this records the physical situation it would weigh.
"""
from __future__ import annotations
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import MOON_RADIUS, VACUUM_PERMEABILITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-magnets-map/1'
MAGNETS = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'magnetic_architecture.json'
ATLAS = ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz'
SITES = ROOT / 'research' / 'studies' / 'conservation' / 'sites.json'
TARGETS = ROOT / 'research' / 'studies' / 'conservation' / 'results' / 'science_targets.csv'
OUT = HERE / 'results' / 'magnets_map.json'
CONFIGS = ('regional_four_500', 'regional_four_1000')
FIELD_T = 0.5e-3
OFFSETS_DEG = np.arange(0.0, 90.0, 5.0)
EVIDENCE = ('Geometry of the shield branch\'s screened loop designs laid on the geography atlas, with the straight-wire '
            'field near the cables. No route has been engineered: terrain, crossings and the cables\' depth are open.')
READING_RULE = ('For each loop design and pattern turn: cable_km and cable_over_water_share; corridor_km is the 0.5 mT '
                'line\'s distance from the cable; sites_in_corridor and targets_in_corridor name what lies within it; '
                'cap_water_share is the water inside the loops. Longitudes are east, latitudes north.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def loop_points(loop, turn_deg, n=720):
    """Points (lat, lon) along a loop on the surface, its centre direction turned about the polar axis."""
    c = np.array(loop['centre']); nrm = np.array(loop['normal']); nrm /= np.linalg.norm(nrm)
    t = np.radians(turn_deg)
    rot = np.array([[np.cos(t), -np.sin(t), 0.0], [np.sin(t), np.cos(t), 0.0], [0.0, 0.0, 1.0]])
    c, nrm = rot @ c, rot @ nrm
    a = np.cross(nrm, [0.0, 0.0, 1.0]); a = a / np.linalg.norm(a) if np.linalg.norm(a) > 1e-9 else np.array([1.0, 0.0, 0.0])
    b = np.cross(nrm, a)
    s = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    p = c[None, :] + loop['radius'] * (np.cos(s)[:, None] * a + np.sin(s)[:, None] * b)
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    return np.degrees(np.arcsin(p[:, 2])), np.degrees(np.arctan2(p[:, 1], p[:, 0])) % 360.0, c / np.linalg.norm(c)


def main(argv=None) -> int:
    magnets = json.loads(MAGNETS.read_text())
    designs = {r['name']: r for r in magnets['regional'] if r['name'] in CONFIGS}
    atlas = np.load(ATLAS)
    water = atlas['water_label'] > 0
    lat_axis, lon_axis = atlas['lat_deg'], atlas['lon_deg']
    def wet(lat, lon):
        i = np.clip(np.round((lat_axis[0] - lat) / (lat_axis[0] - lat_axis[1])).astype(int), 0, lat_axis.size - 1)
        j = np.clip(np.round((lon - lon_axis[0]) / (lon_axis[1] - lon_axis[0])).astype(int) % lon_axis.size, 0, lon_axis.size - 1)
        return water[i, j]
    sites = [(s['name'], s['lat'], s['lon'] % 360.0) for s in json.loads(SITES.read_text())['sites']]
    targets = [(r['name'], float(r['lat']), float(r['lon']) % 360.0) for r in csv.DictReader(open(TARGETS))]
    r_km = MOON_RADIUS / 1e3
    out = {}
    for name, design in designs.items():
        current = float(design['loops'][0]['current'])
        corridor_km = VACUUM_PERMEABILITY * current / (2.0 * np.pi * FIELD_T) / 1e3
        rows = []
        for turn in OFFSETS_DEG:
            cable, over_water, inside, in_site, in_target, cap_cells, cap_wet = 0.0, 0.0, [], set(), set(), 0, 0
            for loop in design['loops']:
                lat, lon, centre = loop_points(loop, turn)
                seg = unit(lat, lon)
                step = np.linalg.norm(np.diff(np.vstack([seg, seg[:1]]), axis=0), axis=1) * r_km
                w = wet(lat, lon)
                cable += step.sum(); over_water += (step * w).sum()
                for label, items, bucket in (('site', sites, in_site), ('target', targets, in_target)):
                    for nm, la, lo in items:
                        d = np.min(np.arccos(np.clip(seg @ unit(la, lo), -1, 1))) * r_km
                        if d < corridor_km:
                            bucket.add(nm)
                ang = np.arcsin(loop['radius'] / MOON_RADIUS)
                LA, LO = np.meshgrid(lat_axis, lon_axis, indexing='ij')
                inside_cap = np.arccos(np.clip(unit(LA, LO) @ centre, -1, 1)) < ang
                cap_cells += int(inside_cap.sum()); cap_wet += int((inside_cap & water).sum())
            rows.append(dict(turn_deg=float(turn), cable_km=round(cable, 0), cable_over_water_share=round(over_water / cable, 3),
                             sites_in_corridor=sorted(in_site), targets_in_corridor=sorted(in_target),
                             cap_water_share=round(cap_wet / max(cap_cells, 1), 3)))
        fewest = min(rows, key=lambda r: (len(r['sites_in_corridor']), len(r['targets_in_corridor']), r['cable_over_water_share']))
        out[name] = dict(loop_radius_km=design['radius_m'] / 1e3, current_a=current, corridor_km=round(corridor_km, 1),
                         turns=rows, fewest_sites_and_targets=fewest['turn_deg'])
    product = dict(schema=SCHEMA, producer=dict(study='joint_synthesis', files={'research/studies/joint_synthesis/magnets_map.py': digest(__file__)},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (MAGNETS, ATLAS, SITES, TARGETS)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, field_t=FIELD_T, designs=out)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for name, d in out.items():
        print(f"{name}: corridor {d['corridor_km']} km either side; fewest sites/targets at a {d['fewest_sites_and_targets']:g} deg turn")
        for r in d['turns'][::3]:
            print(f"   turn {r['turn_deg']:4.0f}: cable {r['cable_km']:6.0f} km, {r['cable_over_water_share']:.0%} over water, "
                  f"caps {r['cap_water_share']:.0%} water, sites {len(r['sites_in_corridor'])} {r['sites_in_corridor'][:4]}, targets {len(r['targets_in_corridor'])}")
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
