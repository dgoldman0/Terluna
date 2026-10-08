"""Lightning over the whole Moon: the electrified box's flash rate carried to every GCM cell by the GCM's convective
rain.

    python -m atmosphere.electricity.occurrence
    # -> atmosphere/electricity/results/occurrence_A28_dim5_moon.json

The electrified CM1 box (climate/crm, box_0e_elec_corrected) ran two lunar days on 0 N, 0 E, forced by the GCM's
column there, and flashed 530 times. The box gives flashes per unit of the GCM's convective rain at its own site,
the two T21 cells on the equator at 0 E (climate/results/gcm/convection_A28_dim5_moon.json). Every other cell
takes that rate per unit rain. Three scalings bracket how lightning follows convection, since the 3-day mean rain
hides how deep the storms grow: linear in the mean rain, as the mean of its 1.5 power, and as the mean of its
square, each fixed to the box at its site. Sea cells take the land rate per unit rain, or a tenth of it, which
brackets the land-sea contrast of Earth's storms. Ground strikes keep the box's share. The flash sizes are the
box's everywhere, so a consumer scales any per-flash quantity by the flash ratio.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.lightning-occurrence/1'
CONVECTION = ROOT / 'climate' / 'results' / 'gcm' / 'convection_A28_dim5_moon.json'
BOX = ROOT / 'climate' / 'results' / 'crm' / 'elec_box_0e_elec_corrected.json'
BOX_WEATHER = ROOT / 'climate' / 'results' / 'crm' / 'box_box_0e_elec_corrected.json'
OUT = HERE / 'results' / 'occurrence_A28_dim5_moon.json'
SCALINGS = {'linear': 'prc_mean', 'power_1_5': 'prc_p15_mean', 'square': 'prc_p2_mean'}
SEA_FACTORS = (1.0, 0.1)
SECONDS_PER_YEAR = 365.25 * 86400.0
SITES = {                                                # lat, east lon; where each is studied
    'storm box (climate/crm box_0e)': (0.0, 0.0),
    'summit port (5.4 N, 158.6 W)': (5.4125, 201.3665),
    'Oceanus Procellarum by Russell (sea appearance)': (32.9, 284.6),
    'eastern Smythii headland (sea appearance)': (-1.5, 90.0),
    'southern Mare Nubium (sea appearance)': (-25.0, 345.0),
    'South Pole-Aitken coast by Mare Ingenii (sea appearance)': (-34.0, 163.0),
    'highland box (climate/crm box_highland)': (-44.7, 246.1),
    'polar plain (75 N, 0 E)': (75.0, 0.0),
}
EVIDENCE = ('One electrified CM1 box (6-km cells, two lunar days, one equatorial lowland site) calibrated against '
            'the GCM\'s convective rain at its own site and carried to the whole Moon by that rain. The scalings '
            'and the sea factor bracket the translation from 3-day mean rain to lightning; storms the GCM\'s '
            'convection scheme misses, terrain and the seas\' own storm structure are outside it.')
READING_RULE = ('Flash densities in flashes per km2 per year; maps run [lat][lon] on the GCM grid. cases[scaling]'
                '[sea_factor] holds the Moon\'s mean, its land and sea means, its flashes per second and the site '
                'values; ratio_to_box multiplies any per-flash or per-area quantity of the box (charge, energy, '
                'fixed nitrogen) to the place. central is the linear scaling with sea storms at the land rate.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def bilinear(lat_deg, lon_deg, field, lat, lon):
    """Value of a [lat][lon] field at a point, between the grid's cells (longitude periodic)."""
    lat_axis, lon_axis = np.asarray(lat), np.asarray(lon)
    order = np.argsort(lat_axis)
    la, f = lat_axis[order], np.asarray(field)[order]
    lat_deg = float(np.clip(lat_deg, la[0], la[-1]))
    i = int(np.clip(np.searchsorted(la, lat_deg) - 1, 0, la.size - 2))
    t = (lat_deg - la[i]) / (la[i + 1] - la[i])
    step = lon_axis[1] - lon_axis[0]
    x = (lon_deg % 360.0) / step
    j0 = int(np.floor(x)) % lon_axis.size
    j1, s = (j0 + 1) % lon_axis.size, x - np.floor(x)
    row = lambda k: (1 - s) * f[k, j0] + s * f[k, j1]
    return float((1 - t) * row(i) + t * row(i + 1))


def occurrence() -> dict:
    conv = json.loads(CONVECTION.read_text())
    box = json.loads(BOX.read_text())
    weather = json.loads(BOX_WEATHER.read_text())
    lat, lon = np.asarray(conv['lat_deg']), np.asarray(conv['lon_deg'])
    share = np.asarray(conv['cell_area_share'])
    land = np.asarray(conv['land_share']) > 0.5
    flashes = box['lightning']['all']['count']
    ground = box['lightning']['negative_to_ground']['count']
    days = (box['end_s'] - box['start_s']) / 86400.0
    cells = weather['gcm']['land']['cells']                      # the T21 cells the box was forced from
    box_km2 = (64 * 6.01124e3 / 1e3) ** 2                        # 64 x 64 columns of 6,011.24 m (case.json grid)
    rate0 = flashes / box_km2 / (days * 86400.0 / SECONDS_PER_YEAR)
    moon_km2 = 4.0 * np.pi * (MOON_RADIUS / 1e3) ** 2
    cases = {}
    for name, key in SCALINGS.items():
        field = np.asarray(conv[key])
        at_box = np.mean([bilinear(c[0], c[1], field, lat, lon) for c in cells])
        cases[name] = {}
        for sea in SEA_FACTORS:
            ratio = field / at_box * np.where(land, 1.0, sea)
            density = rate0 * ratio
            mean = float((share * density).sum())
            cases[name][f'{sea:g}'] = dict(
                moon_mean=round(mean, 5), land_mean=round(float((share * density * land).sum() / (share * land).sum()), 5),
                sea_mean=round(float((share * density * ~land).sum() / (share * ~land).sum()), 5),
                flashes_per_s=round(mean * moon_km2 / SECONDS_PER_YEAR, 3),
                ground_per_s=round(mean * moon_km2 / SECONDS_PER_YEAR * ground / flashes, 4),
                ratio_to_box_moon_mean=round(mean / rate0, 4),
                sites={s: dict(flash_density=round(bilinear(p[0], p[1], density, lat, lon), 5),
                               ratio_to_box=round(bilinear(p[0], p[1], ratio, lat, lon), 4))
                       for s, p in SITES.items()})
    central_ratio = np.asarray(conv['prc_mean']) / np.mean([bilinear(c[0], c[1], conv['prc_mean'], lat, lon) for c in cells])
    return dict(box=dict(flashes=flashes, ground_strikes=ground, days=round(days, 2), area_km2=round(box_km2, 0),
                         flash_density=round(rate0, 5), gcm_cells=cells),
                lat_deg=conv['lat_deg'], lon_deg=conv['lon_deg'],
                central_flash_density=np.round(rate0 * central_ratio, 5).tolist(),
                central_ratio_to_box=np.round(central_ratio, 4).tolist(), cases=cases)


def main(argv=None) -> int:
    result = occurrence()
    product = dict(schema=SCHEMA, producer=dict(domain='atmosphere', files={'electricity/occurrence.py': digest(__file__)},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (CONVECTION, BOX, BOX_WEATHER)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, **result)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    b = result['box']
    print(f"box: {b['flashes']} flashes in {b['days']} days over {b['area_km2']:.0f} km2 = {b['flash_density']} per km2 a year")
    for name, by_sea in result['cases'].items():
        for sea, c in by_sea.items():
            print(f"{name:10s} sea x{sea:4s}: Moon {c['moon_mean']:.4f} (land {c['land_mean']:.4f}, sea {c['sea_mean']:.4f}) "
                  f"per km2 a year; {c['flashes_per_s']:.2f} flashes/s; summit {c['sites']['summit port (5.4 N, 158.6 W)']['flash_density']:.4f}")
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
