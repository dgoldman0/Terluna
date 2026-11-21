"""The atlas round's first map of candidates for satoyama-like landscapes, on the drainage of the redrawn floor and the
present design climate.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.satoyama_candidates     # about 5 minutes

satoyama.md sets the criteria from the Earth evidence; this map tests them at the drainage's 1.9 km cells:

- the units are side valleys: the land that drains through a cell draining 10-100 km2 into a river that drains more
  than 100 km2, a lake or a sea, the scale of a tank cascade or a village's valley. The large rivers' own corridors
  and slopes of under 10 km2 that run straight to a river or the shore are left out;
- a unit has the order down the slope where gentle valley floor, under 1 in 20, makes up at least 3% of it in the
  lower half of its relief, the irrigated share of the West Javan hill village, and land steeper than 15 degrees or in
  the upper half of its relief makes up at least 30% of it for the woodland above;
- its rain lies between 800 and 3,000 mm a year on the design climate;
- the flooded fields it could carry are the least of its valley floor, the land its own discharge could water for one
  rice season a year beyond the rain that falls in the season, and 28% of it, Hani's terraced share;
- their methane is that area's at the IPCC's rate for fields drained more than once a season and cropped once a year,
  set against the room under the held 1.87 ppm; the same fields with upland crops before the rice give the low end.

The salt south sea's shallows are measured for satoumi beside the candidates. The hybrid land is a fifth of each zone's
land, unmapped, so the candidates' land is set against a fifth of the land in each band of latitude on the turned Moon.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from geography import drainage as dr
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'satoyama_candidates.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/satoyama_candidates')
SCHEMA = 'terluna.research.atlas-round-satoyama-candidates/1'
CODE = ('research/studies/atlas_round/satoyama_candidates.py', 'research/studies/atlas_round/redrawn_drainage.py',
        'geography/drainage.py')
AREA_KM2 = (10.0, 100.0)                # a unit drains this much at its outlet
GENTLE = 0.05                           # paddy blocks gentle to 1 in 20
STEEP_DEG = 15.0                        # rain-fed fields gentle to 15 degrees, steep beyond
VALLEY_SHARE_MIN = 0.03                 # the West Javan village's irrigated share, the least of the three analogues
WOODLAND_SHARE_MIN = 0.30               # land above the fields for the woodland, against Hani's half and Java's 41%
RAIN_MM = (800.0, 3000.0)
SEASON_DAYS = 113.0                     # the IPCC's default season, as satoyama.md takes it
NEED_MM = (1240.0, 1880.0)              # a flooded paddy's water for a season
IRRIGATED_SHARE_MAX = 0.28              # Hani's terraces
SRI_LANKA_SHARE = 0.06                  # the heritage site's paddy
BASE_KG_HA_DAY = 1.19                   # IPCC 2019 Table 5.11, continuously flooded, no organic inputs
DRAINED_MORE_THAN_ONCE = 0.55
DRY_BEFORE = 0.59                       # rice after upland crops or dry fallow
POND_T_KM2_YR = 18.3                    # 183 kg/ha/yr
CHANNEL_T_KM2_YR = 41.6                 # 416 kg/ha/yr
ROOM_TG_YR = (0.35, 2.5)                # under 1.87 ppm, the removal ledger's high- and low-cost ends (lake methane study)
BANDS = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
SHALLOWS_M = (2.0, 5.0, 10.0, 20.0)
TOP = 25


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def rate_t_km2(drained=DRAINED_MORE_THAN_ONCE, dry_before=1.0, days=SEASON_DAYS):
    """Methane from a square kilometre of paddy cropped once a year (t a year)."""
    return BASE_KG_HA_DAY * drained * dry_before * days * 100.0 / 1000.0


def slopes(h, lat, cell_km):
    """The ground's slope (m per m) at each 16 px/deg cell, by central differences."""
    dy = cell_km * 1000.0
    dx = dy * np.maximum(np.cos(np.radians(lat)), 1e-3)[:, None]
    gx = (np.roll(h, -1, axis=1) - np.roll(h, 1, axis=1)) / (2.0 * dx)
    gy = np.gradient(h, axis=0) / dy
    return np.hypot(gx, gy).astype(np.float32)


def units(receiver, drained_km2, water, lo=AREA_KM2[0], hi=AREA_KM2[1]):
    """Each land cell's side valley: the cell its water leaves the valley by, a cell draining at most `hi` km2 whose water
    goes next to a river draining more, to water, or nowhere further (a closed hollow), or -1 for cells on the large
    rivers and on water. Valleys draining less than `lo` at their outlet are dropped (-1)."""
    n = receiver.size
    cells = np.arange(n)
    small = ~water & (drained_km2 <= hi)
    nxt = receiver.astype(np.int64)
    after = drained_km2[nxt]
    outlet = small & (water[nxt] | (after > hi) | (after < drained_km2) | (nxt == cells))
    step = np.where(small & ~outlet, nxt, cells)
    root = dr.terminal(step)
    unit = np.where(small, root, -1)
    unit[small & ~outlet[root]] = -1                     # a path that left the small land without an outlet
    keep = np.zeros(n, bool)
    keep[cells[outlet & (drained_km2 >= lo)]] = True
    unit[(unit >= 0) & ~keep[np.maximum(unit, 0)]] = -1
    return unit


def per_unit(unit, values, weights=None, n=None):
    """Sums of values (times weights) over each unit, indexed by the unit's position in `n` compact ids."""
    w = values if weights is None else values * weights
    return np.bincount(unit, weights=w, minlength=n)


def candidates(h, lat, cell_km, routing, climate, rotation, log=print):
    """The units, their measures and which meet the criteria."""
    started = time.time()
    rows, cols = h.shape
    hf = h.ravel()
    area = np.repeat((cell_km * 1e3) ** 2 * np.cos(np.radians(lat)), cols) / 1e6        # km2
    receiver, drained, q = routing['receiver'], routing['drained_km2'], routing['discharge_m3s']
    lake = np.zeros(hf.size, bool)
    lake[routing['lake_cell']] = True
    water = (routing['sea'] >= 0) | lake
    unit_cell = units(receiver, drained, water)
    ids, compact = np.unique(unit_cell[unit_cell >= 0], return_inverse=True)
    n = ids.size
    member = np.flatnonzero(unit_cell >= 0)
    u = compact
    log(f'{n:,} side valleys of {AREA_KM2[0]:g}-{AREA_KM2[1]:g} km2 over {area[member].sum() / 1e6:.2f} million km2, '
        f'{time.time() - started:.0f} s')
    a = area[member]
    s = slopes(h, lat, cell_km).ravel()[member]
    z = hf[member].astype(float)
    order = np.argsort(u, kind='stable')
    starts = np.r_[0, np.flatnonzero(np.diff(u[order])) + 1]
    zmin = np.minimum.reduceat(z[order], starts)
    zmax = np.maximum.reduceat(z[order], starts)
    relief = zmax - zmin
    frac = (z - zmin[u]) / np.maximum(relief[u], 1.0)
    unit_area = per_unit(u, a, n=n)
    valley = per_unit(u, a * ((s < GENTLE) & (frac < 0.5)), n=n)
    wood = per_unit(u, a * ((s > np.tan(np.radians(STEEP_DEG))) | (frac >= 0.5)), n=n)
    rain = per_unit(u, a * climate['pr'].ravel()[member], n=n) / unit_area * 365.25          # mm a year
    temp = per_unit(u, a * climate['ts'].ravel()[member], n=n) / unit_area - 273.15
    row = per_unit(u, a * (member // cols), n=n) / unit_area
    col_angle = np.angle(per_unit(u, a * np.exp(2j * np.pi * (member % cols) / cols).real, n=n)
                         + 1j * per_unit(u, a * np.exp(2j * np.pi * (member % cols) / cols).imag, n=n))
    col = (col_angle % (2 * np.pi)) / (2 * np.pi) * cols
    lat_c = 90.0 - (row + 0.5) / dr.PPD
    lon_c = ((col + 0.5) / dr.PPD + 180.0) % 360.0 - 180.0
    from climate.gcm.turned import to_turned
    t_lat, t_lon = to_turned(lat_c, lon_c, rotation)
    outflow = q[ids]                                     # m3/s leaving each valley
    season_rain = rain * SEASON_DAYS / 365.25
    net = {k: np.maximum(need - season_rain, 0.0) for k, need in (('low', NEED_MM[0]), ('high', NEED_MM[1]))}
    seconds = 365.25 * 86400.0
    watered = {k: np.where(v > 0, outflow * seconds * 1000.0 / np.maximum(v, 1e-9) / 1e6, np.inf)     # km2 of paddy
               for k, v in net.items()}
    possible = np.minimum(np.minimum(valley, watered['high']), IRRIGATED_SHARE_MAX * unit_area)
    order_ok = (valley >= VALLEY_SHARE_MIN * unit_area) & (wood >= WOODLAND_SHARE_MIN * unit_area)
    rain_ok = (rain >= RAIN_MM[0]) & (rain <= RAIN_MM[1])
    ok = order_ok & rain_ok
    log(f'{int(ok.sum()):,} meet the criteria, {time.time() - started:.0f} s')
    return dict(ids=ids, unit_cell=unit_cell, area=unit_area, valley=valley, wood=wood, relief=relief, rain=rain, temp=temp,
                lat=lat_c, lon=lon_c, t_lat=t_lat, t_lon=t_lon, outflow=outflow, watered=watered, possible=possible,
                order_ok=order_ok, rain_ok=rain_ok, ok=ok, cell_area=area, water=water)


def satoumi(h, level0, routing, area, lat, cols, rotation):
    """The salt south sea's shallows below the common level, by depth, with the length of its shore."""
    seas = routing['meta']['seas']
    south = routing['sea'] == seas.index('south_pole_aitken')
    depth = level0 - h.ravel()
    out = dict(sea_km2=round(float(area[south].sum())))
    for d in SHALLOWS_M:
        out[f'shallower_than_{d:g}_m_km2'] = round(float(area[south & (depth < d)].sum()))
    land = (routing['sea'] < 0).reshape(-1, cols)
    sea2 = south.reshape(-1, cols)
    edge = np.zeros_like(sea2)
    for dr_, dc in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        edge |= sea2 & dr.shifted(land, dr_, dc, False)
    cell_side = np.sqrt(area.reshape(-1, cols))
    out['shore_km'] = round(float(cell_side[edge].sum()))
    from climate.gcm.turned import to_turned
    shallow = np.flatnonzero(south & (depth < 10.0))
    t_lat, _ = to_turned(90.0 - (shallow // cols + 0.5) / dr.PPD, ((shallow % cols) + 0.5) / dr.PPD, rotation)
    bands = {}
    for lo_, hi_ in zip(BANDS[:-1], BANDS[1:]):
        m = (np.abs(t_lat) >= lo_) & (np.abs(t_lat) < hi_)
        if m.any():
            bands[f'{lo_:g}-{hi_:g}'] = round(float(area[shallow[m]].sum()))
    out['shallower_than_10_m_km2_by_turned_latitude'] = bands
    return out


def record(c, routing, level0, area, sat, rotation):
    """Totals, the bands, the methane against the room and the leading candidates."""
    ok = c['ok']
    central, low_end = rate_t_km2(), rate_t_km2(dry_before=DRY_BEFORE)
    land = ~c['water']
    possible = c['possible'][ok]
    t_lat = np.abs(c['t_lat'])
    total_land = float(area[land].sum())
    bands = []
    rows = area.size // routing['meta']['grid']['cols']
    cols = routing['meta']['grid']['cols']
    cell_rows = np.arange(area.size) // cols
    cell_lat = 90.0 - (cell_rows + 0.5) / dr.PPD
    from climate.gcm.turned import to_turned
    cell_t_lat = np.empty(area.size, np.float32)
    for start in range(0, rows, 360):                    # the turned latitude of every cell, in blocks of rows
        sl = slice(start * cols, min(rows, start + 360) * cols)
        lon = ((np.arange(sl.start, sl.stop) % cols) + 0.5) / dr.PPD
        cell_t_lat[sl] = to_turned(cell_lat[sl], lon, rotation)[0]
    for lo_, hi_ in zip(BANDS[:-1], BANDS[1:]):
        in_band = (t_lat >= lo_) & (t_lat < hi_)
        cells = land & (np.abs(cell_t_lat) >= lo_) & (np.abs(cell_t_lat) < hi_)
        band_land = float(area[cells].sum())
        cand = float(c['area'][ok & in_band].sum())
        bands.append(dict(turned_latitude=f'{lo_:g}-{hi_:g}', land_km2=round(band_land),
                          valleys_km2=round(float(c['area'][in_band].sum())), candidates=int((ok & in_band).sum()),
                          candidates_km2=round(cand), share_of_the_hybrid_fifth=round(cand / (0.2 * band_land), 3) if band_land else None,
                          flooded_possible_km2=round(float(c['possible'][ok & in_band].sum())),
                          rain_mm_median=round(float(np.median(c['rain'][in_band]))) if in_band.any() else None))
    room_km2 = {f'{r:g}_tg': round(r * 1e6 / central) for r in ROOM_TG_YR}
    top = np.flatnonzero(ok)[np.argsort(c['possible'][ok])[::-1][:TOP]]
    cand_area = float(c['area'][ok].sum())
    return dict(
        valleys=dict(count=int(c['ids'].size), km2=round(float(c['area'].sum())),
                     share_of_the_land=round(float(c['area'].sum()) / total_land, 3)),
        land_km2=round(total_land),
        failing=dict(order_down_the_slope=int((~c['order_ok']).sum()), rain=int((~c['rain_ok']).sum())),
        candidates=dict(count=int(ok.sum()), km2=round(cand_area), share_of_the_land=round(cand_area / total_land, 3),
                        share_of_a_fifth_of_the_land=round(cand_area / (0.2 * total_land), 3),
                        valley_floor_km2=round(float(c['valley'][ok].sum())),
                        woodland_km2=round(float(c['wood'][ok].sum())),
                        flooded_possible_km2=round(float(possible.sum())),
                        flooded_at_sri_lankas_share_km2=round(SRI_LANKA_SHARE * cand_area),
                        water_limited=int((c['watered']['high'][ok] < np.minimum(c['valley'][ok], IRRIGATED_SHARE_MAX * c['area'][ok])).sum()),
                        rain_mm_p10_p50_p90=[round(float(x)) for x in np.percentile(c['rain'][ok], [10, 50, 90])] if ok.any() else None,
                        temperature_c_p10_p50_p90=[round(float(x), 1) for x in np.percentile(c['temp'][ok], [10, 50, 90])] if ok.any() else None),
        methane=dict(
            paddy_t_km2_yr=dict(drained_more_than_once_cropped_once=round(central, 2), with_upland_crops_before=round(low_end, 2),
                                continuously_flooded_cropped_twice=round(rate_t_km2(drained=1.0) * 2, 1)),
            ponds_t_km2_yr=POND_T_KM2_YR, channels_t_km2_yr=CHANNEL_T_KM2_YR,
            room_tg_yr=list(ROOM_TG_YR), room_holds_paddy_km2=room_km2,
            at_the_possible_area_tg_yr=[round(float(possible.sum()) * x / 1e6, 2) for x in (low_end, central)],
            at_sri_lankas_share_tg_yr=[round(SRI_LANKA_SHARE * cand_area * x / 1e6, 2) for x in (low_end, central)],
            share_of_the_possible_area_the_room_holds=[round(min(1.0, r * 1e6 / central / max(float(possible.sum()), 1.0)), 3)
                                                       for r in ROOM_TG_YR]),
        bands=bands, satoumi=sat,
        leading=[dict(today=[round(float(c['lat'][i]), 2), round(float(c['lon'][i]), 2)],
                      turned=[round(float(c['t_lat'][i]), 2), round(float(c['t_lon'][i]), 2)],
                      km2=round(float(c['area'][i]), 1), relief_m=round(float(c['relief'][i])),
                      valley_floor_km2=round(float(c['valley'][i]), 1), woodland_share=round(float(c['wood'][i] / c['area'][i]), 2),
                      flooded_possible_km2=round(float(c['possible'][i]), 1), rain_mm=round(float(c['rain'][i])),
                      temperature_c=round(float(c['temp'][i]), 1), outflow_m3s=round(float(c['outflow'][i]), 2))
                 for i in top])


def picture(c, routing, h, path, factor=4):
    """A quick-look map at 4 px/deg: the seas and lakes, the valleys that fail and the candidates by their woodland share."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows, cols = h.shape
    cls = np.zeros(rows * cols, np.int8)                 # 0 land outside the valleys
    unit_cell = c['unit_cell']
    pos = np.searchsorted(c['ids'], unit_cell[unit_cell >= 0])
    good = c['ok'][pos]
    member = np.flatnonzero(unit_cell >= 0)
    cls[member] = np.where(good, 3, 2)
    cls[c['water']] = 1
    cls = cls.reshape(rows, cols)[::factor, ::factor]
    cls = np.roll(cls, cls.shape[1] // 2, axis=1)         # columns start at 0 E; the picture centres it
    colours = np.array([[0.86, 0.83, 0.76], [0.36, 0.55, 0.75], [0.74, 0.70, 0.60], [0.25, 0.55, 0.28]])
    img = colours[cls]
    fig = plt.figure(figsize=(14, 7.4), dpi=110)
    ax = fig.add_axes([0.02, 0.06, 0.96, 0.9])
    ax.imshow(img, extent=(-180, 180, -90, 90), interpolation='nearest')
    ax.set_title('Satoyama-like candidates on the redrawn floor (green), other side valleys (tan), water (blue); '
                 'today\'s Moon, 0 E at the centre', fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def build(log=print):
    started = time.time()
    surf = RD.load_surface(log)
    h, level0, lat, cell_km = surf['h'], surf['level0'], surf['lat'], surf['cell_km']
    rows, cols = h.shape
    with np.load(RD.DESIGN_ROUTING, allow_pickle=False) as d:
        routing = {k: d[k] for k in d.files}
    routing['meta'] = json.loads(str(routing.pop('metadata')))
    product = json.loads(RD.OUT.read_text())
    if sha256(RD.DESIGN_ROUTING) != product['design_routing']['sha256']:
        raise SystemExit(f'{RD.DESIGN_ROUTING} does not match {RD.OUT.name}')
    clim = dict(np.load(ROOT / NW.CLIMATOLOGY))
    clim.pop('metadata')
    climate = dict(pr=dr.climate_on_grid(clim, 'pr_mm_day', rows, cols), ts=dr.climate_on_grid(clim, 'ts_k', rows, cols))
    rotation = np.array(json.loads((ROOT / NW.FRAME).read_text())['frame']['rotation_today_to_turned'])
    c = candidates(h, lat, cell_km, routing, climate, rotation, log)
    del climate
    area = c['cell_area']
    sat = satoumi(h, level0, routing, area, lat, cols, rotation)
    out = record(c, routing, level0, area, sat, rotation)
    BULK.mkdir(parents=True, exist_ok=True)
    picture(c, routing, h, BULK / 'candidates_map.png')
    pos = np.full(rows * cols, -1, np.int32)
    members = c['unit_cell'] >= 0
    pos[members] = np.searchsorted(c['ids'], c['unit_cell'][members]).astype(np.int32)
    np.savez_compressed(BULK / 'satoyama_units.npz', unit=pos, candidate=c['ok'], outlet_cell=c['ids'].astype(np.int32),
                        flooded_possible_km2=c['possible'].astype(np.float32), woodland_km2=c['wood'].astype(np.float32),
                        valley_floor_km2=c['valley'].astype(np.float32), area_km2=c['area'].astype(np.float32))
    inputs = [RD.OUT, ROOT / NW.FRAME]
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in inputs},
        ignored_inputs={str(RD.DESIGN_ROUTING): sha256(RD.DESIGN_ROUTING), rel(ROOT / NW.CLIMATOLOGY): sha256(ROOT / NW.CLIMATOLOGY)},
        evidence=('A screen on products: side valleys of the design routing on the redrawn floor at 1.9 km cells, tested '
                  'against criteria set from Earth\'s satoyama, tank cascades and terraces (satoyama.md), with the present '
                  'design climate\'s rain and surface temperature at T21 and the IPCC 2019 Refinement\'s methane factors for '
                  'rice. Slopes at 1.9 km are gentler than the ground within a cell, so the valley floors and the woodland '
                  'are first estimates. Not validated against any observation.'),
        reading_rule=('Areas in km2, rain in mm a year, temperatures the climate model\'s surface temperature in C, methane in '
                      'Tg a year. A candidate is a side valley that meets the order down the slope and the rain; '
                      'flooded_possible is the least of its valley floor, the paddy its outflow waters for one season a year '
                      'at the season\'s higher need, and 28% of it. share_of_the_hybrid_fifth sets the candidates\' land '
                      'against a fifth of the land in each band of latitude on the turned Moon. Places are given on today\'s '
                      'Moon and the turned Moon at 8 degrees.'),
        criteria=dict(valley_km2=list(AREA_KM2), gentle_slope=GENTLE, steep_deg=STEEP_DEG, valley_floor_share_min=VALLEY_SHARE_MIN,
                      woodland_share_min=WOODLAND_SHARE_MIN, rain_mm=list(RAIN_MM), season_days=SEASON_DAYS, need_mm=list(NEED_MM),
                      irrigated_share_max=IRRIGATED_SHARE_MAX),
        **out, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    result = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    print(f'wrote {args.out} in {result["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
