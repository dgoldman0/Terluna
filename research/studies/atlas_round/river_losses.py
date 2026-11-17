"""What the rivers of the atlas round's governed routing would lose in the drylands, and whether any wetland they could
keep is large enough to belong in the GCM's surface beside the lakes.

    ulimit -v 12000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.river_losses   # a few minutes

PlaSim carries its rivers' water to the sea and loses none on the way, while a river crossing dry land loses water
through its bed, over its floodplain and to the air, and keeps a green corridor or a wetland in return. The governed
routing (governed_drainage.py, the routing's own budget, the design climate's rain) holds the lakes the rivers fill,
which the turned surface carries as water; its channels and floodplains are not in the surface. This estimates, on the
redrawn floor's 16 px/deg cells, in the drylands, where open water loses more to the air than the rain brings it (the
routing's own balance, rain less the open water's evaporation at each latitude, lake_net < 0):

- the channels: every cell a river of at least Q_RIVER crosses, its water surface as wide as Q^0.5 times Earth's range
  of coefficients (WIDTH), widened for lunar gravity by (g_Earth / g_Moon)^0.2 (the exponent of gravity in Parker and
  others' dimensionless widths), over the cell's flow length, losing the open water's net evaporation;
- the riparian belt: the plants along the channel drinking the water the bed lets into the ground, as a belt of
  RIPARIAN_KM in all evaporating at the open water's net rate, an upper reading of the bed's seepage;
- the floodplains: dry land within DIST_KM of a river whose height stands 0 to HAND_M metres above the nearest river
  cell's, which a river could wet; its loss if all of it were wet, and capped by what each river can give, since a
  wetland's area settles where its loss meets its inflow (an inland delta's does);
- the GCM test: each floodplain, capped, carried into the settled frame and onto the T21 grid the run on the final
  surface stands on, against the land fraction below which PlaSim's mask makes a cell water.

Earth's alluvial rivers span the width coefficients; the lunar widening, the belt and the floodplain heights are ranges,
not measurements. Heights come from LOLA at 1.9 km cells, so a floodplain narrower than a few cells goes unseen.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from climate.gcm import boundary as bd
from climate.gcm import turned
from geography import drainage as dr
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import turned_surface as TS
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'river_losses.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/river_losses')
GREEN = BULK / 'river_green.npz'      # the rivers' cells and floodplains, for the drylands' reading (drylands.py)
SCHEMA = 'terluna.research.atlas-round-river-losses/1'
CODE = ('research/studies/atlas_round/river_losses.py', 'research/studies/atlas_round/redrawn_drainage.py',
        'research/studies/atlas_round/governed_drainage.py', 'research/studies/atlas_round/turned_surface.py')
Q_RIVER = 100.0                      # m3/s: the rivers counted
Q_CLASSES = (100.0, 1000.0, 10000.0)
WIDTH = (3.0, 15.0)                  # m per (m3/s)^0.5: Earth's alluvial rivers, single channels to braided
G_EARTH, G_MOON = 9.80665, 1.62
LUNAR_WIDENING = (G_EARTH / G_MOON) ** 0.2
RIPARIAN_KM = (0.5, 2.0)             # the belt of plants drinking the bed's seepage, both banks together
DIST_KM = 10.0                       # how far from its river a floodplain cell may lie
HAND_M = (2.0, 5.0, 10.0)            # how far above its river a floodplain cell may stand
HAND_MAIN = 5.0
PAD = 64                             # columns carried across the 0/360 seam for the nearest river
TOP = 12                             # the largest floodplains listed
SECONDS_PER_YEAR = 365.25 * 86400.0


def flow_length_km(receiver, rows, cols, lat, cell_km):
    """Each cell's flow length to its receiver (km); nil at a pit."""
    idx = np.arange(receiver.size)
    r0, c0 = idx // cols, idx % cols
    r1, c1 = receiver // cols, receiver % cols
    dc = (c1 - c0 + cols // 2) % cols - cols // 2
    cos = np.cos(np.radians(lat))[r0]
    return np.where(receiver == idx, 0.0, cell_km * np.sqrt((r1 - r0) ** 2 + (dc * cos) ** 2))


def channel_width_m(q, coefficient):
    """A river's water surface width (m) at discharge q (m3/s): coefficient q^0.5, widened for lunar gravity."""
    return coefficient * LUNAR_WIDENING * np.sqrt(np.maximum(q, 0.0))


def nearest_river(river, pad=PAD):
    """For every cell of the grid, the flat index of the nearest river cell by grid distance, the columns wrapped
    across the 0/360 seam."""
    rows, cols = river.shape
    wide = np.concatenate([river[:, -pad:], river, river[:, :pad]], axis=1)
    _, (ri, ci) = ndimage.distance_transform_edt(~wide, return_indices=True)
    ri, ci = ri[:, pad:pad + cols], (ci[:, pad:pad + cols] - pad) % cols
    return (ri * cols + ci).ravel()


def distance_km(i, j, cols, lat, cell_km):
    """The distance (km) between cells i and j on the sphere, small-angle."""
    ri, ci, rj, cj = i // cols, i % cols, j // cols, j % cols
    dc = (cj - ci + cols // 2) % cols - cols // 2
    cos = np.cos(np.radians(0.5 * (lat[ri] + lat[rj])))
    return cell_km * np.sqrt((rj - ri) ** 2 + (dc * cos) ** 2)


def floodplain(h, near, dry, land, lat, cell_km, hand_m, dist_km=DIST_KM):
    """The floodplain cells (flat bool): dry land (flat bools dry and land) within dist_km of its nearest river cell
    (near, nearest_river's), standing 0 to hand_m above it."""
    rows, cols = h.shape
    cells = np.flatnonzero(dry & land)
    hand = h.ravel()[cells] - h.ravel()[near[cells]]
    ok = (hand >= 0.0) & (hand <= hand_m) & (distance_km(cells, near[cells], cols, lat, cell_km) <= dist_km)
    out = np.zeros(h.size, bool)
    out[cells[ok]] = True
    return out


def capped_area_km2(area_km2, loss_m_s, inflow_m3s):
    """A wetland's area where its loss meets what feeds it: the whole if the inflow covers its loss, else the share
    the inflow can keep wet."""
    need = area_km2 * 1e6 * loss_m_s
    return np.where(need <= inflow_m3s, area_km2, inflow_m3s / np.maximum(loss_m_s, 1e-30) / 1e6)


def t21_cells(lat_deg, lon_deg, rotation):
    """The T21 cell (row, column) of the turned grid each point falls in."""
    lat_t, lon_t = turned.to_turned(lat_deg, lon_deg, rotation)
    _, edges = bd.gaussian_latitudes(32)
    rows = np.clip(np.searchsorted(-edges, -np.sin(np.radians(lat_t)), side='right') - 1, 0, 31)
    cols = np.floor(np.mod(lon_t + 180.0 / 64, 360.0) / (360.0 / 64)).astype(int) % 64
    return rows, cols


def build(log=print):
    started = time.time()
    surf = RD.load_surface(log)
    h = surf['h']
    rows, cols = h.shape
    lat, cell_km = surf['lat'], surf['cell_km']
    with np.load(GD.GOVERNED_ROUTING, allow_pickle=False) as d:
        receiver, q = d['receiver'].astype(np.int64), d['discharge_m3s'].astype(float)
        sea, lake_cell = d['sea'], d['lake_cell']
    governed = json.loads(GD.OUT.read_text())
    factor = governed['cases']['closed_budget']['evaporation']
    area_km2 = np.repeat(cell_km ** 2 * np.cos(np.radians(lat)), cols)
    land = sea < 0
    land[lake_cell] = False
    length = flow_length_km(receiver, rows, cols, lat, cell_km)
    river = land & (q >= Q_RIVER)
    near = nearest_river(river.reshape(h.shape))
    log(f'rivers of at least {Q_RIVER:.0f} m3/s: {int(river.sum())} cells; nearest river found')
    record = dict(budgets={})
    for label, evaporation in (('routing_budget', 1.0), ('closed_budget', factor)):
        _, _, _, _, _, _, lake_net, _, _ = NW.load(1.0, evaporation)
        loss = np.maximum(-lake_net.ravel().astype(float), 0.0) * dr.MM_DAY          # m/s of water the air takes
        del lake_net
        dry = loss > 0.0
        dry_river = river & dry
        out = dict(evaporation_factor=evaporation, rivers={}, floodplains={})
        for qmin in Q_CLASSES:
            m = dry_river & (q >= qmin)
            out['rivers'][f'{qmin:.0f}'] = dict(length_km=round(float(length[m].sum()), 0),
                                                largest_flow_m3s=round(float(q[m].max()), 0) if m.any() else None)
        lengths, losses = length[dry_river], loss[dry_river]
        channel = {}
        for coef in WIDTH:
            width_km = channel_width_m(q[dry_river], coef) / 1000.0
            channel[f'{coef:g}'] = dict(area_km2=round(float((width_km * lengths).sum()), 0),
                                        evaporation_m3s=round(float((width_km * lengths * 1e6 * losses).sum()), 1))
        out['channel'] = channel
        out['riparian_belt'] = {f'{b:g}': dict(area_km2=round(float((b * lengths).sum()), 0),
                                               evaporation_m3s=round(float((b * lengths * 1e6 * losses).sum()), 1))
                                for b in RIPARIAN_KM}
        out['dry_land_share_of_land'] = round(float(area_km2[dry & land].sum() / area_km2[land].sum()), 4)
        for hand in HAND_M:
            plain = floodplain(h, near, dry, land, lat, cell_km, hand)
            whole = float((area_km2[plain] * 1e6 * loss[plain]).sum())
            out['floodplains'][f'{hand:g}'] = dict(area_km2=round(float(area_km2[plain].sum()), 0),
                                                   evaporation_if_all_wet_m3s=round(whole, 1))
            if label == 'routing_budget' and hand == HAND_MAIN:
                main = (plain, loss)
        record['budgets'][label] = out
        log(f'{label}: dry-land rivers {out["rivers"]["100"]["length_km"]:.0f} km, channels '
            f'{channel[f"{WIDTH[0]:g}"]["evaporation_m3s"]}-{channel[f"{WIDTH[1]:g}"]["evaporation_m3s"]} m3/s, '
            f'floodplains (5 m) {out["floodplains"]["5"]["area_km2"]:.0f} km2')
    plain, loss = main
    labels, count = ndimage.label(plain.reshape(h.shape), structure=np.ones((3, 3), bool))
    labels = labels.ravel()
    cells = np.flatnonzero(labels)
    lab = labels[cells]
    area = np.bincount(lab, weights=area_km2[cells], minlength=count + 1)
    need = np.bincount(lab, weights=area_km2[cells] * 1e6 * loss[cells], minlength=count + 1)
    inflow = np.zeros(count + 1)
    np.maximum.at(inflow, lab, q[near[cells]])
    mean_loss = need / np.maximum(area * 1e6, 1e-30)
    capped = capped_area_km2(area, mean_loss, inflow)
    frame = json.loads(TS.record_path(RD.sha256(TS.ST.OUT)).read_text())
    rotation = np.asarray(frame['frame']['rotation_today_to_turned'])
    t21 = bd.read_t21(TS.PRODUCTS / frame['surface'])
    radius_km = 1737.4
    cell_area = t21['area'] / t21['area'].sum() * 4.0 * np.pi * radius_km ** 2
    cut = float(t21['land'][t21['mask'] == 0].max())
    row_c, col_c = cells // cols, cells % cols
    lon = (np.arange(cols) + 0.5) / dr.PPD
    tr, tc = t21_cells(lat[row_c], lon[col_c], rotation)
    keep = capped[lab] / np.maximum(area[lab], 1e-30)                      # each floodplain's wet share, capped
    added = np.zeros(t21['land'].shape)
    np.add.at(added, (tr, tc), area_km2[cells] * keep)
    fraction = added / cell_area
    land_t21 = t21['mask'] > 0
    flips = land_t21 & (t21['land'] - fraction < cut)
    BULK.mkdir(parents=True, exist_ok=True)
    river_all = np.flatnonzero(land & (q >= Q_RIVER))
    np.savez_compressed(GREEN, metadata=json.dumps(dict(
                            schema=SCHEMA + '#green', grid=dict(pixels_per_degree=dr.PPD, rows=rows, cols=cols),
                            reading_rule=('river cells (land, at least Q_RIVER) with their flow length (km) and flow '
                                          '(m3/s); the floodplain cells at 5 m with their area (km2) and the share '
                                          'their river can keep wet; cells index rows from 90 N and columns from 0 E '
                                          'on today\'s Moon'))),
                        river_cell=river_all.astype(np.int32), river_length_km=length[river_all].astype(np.float32),
                        river_flow_m3s=q[river_all].astype(np.float32), plain_cell=cells.astype(np.int32),
                        plain_area_km2=area_km2[cells].astype(np.float32),
                        plain_kept_share=(capped[lab] / np.maximum(area[lab], 1e-30)).astype(np.float32))
    biggest = np.argsort(area[1:])[::-1][:TOP] + 1
    listed = []
    for k in biggest:
        members = cells[lab == k]
        r_, c_ = members // cols, members % cols
        la, lo = float(lat[r_].mean()), float(np.degrees(np.angle(np.exp(1j * np.radians(lon[c_])).mean())))
        lt, ln = turned.to_turned(la, lo, rotation)
        listed.append(dict(area_km2=round(float(area[k]), 0), today=[round(la, 2), round(lo, 2)],
                           turned=[round(float(lt), 2), round(float(ln), 2)], river_flow_m3s=round(float(inflow[k]), 0),
                           net_evaporation_mm_year=round(float(mean_loss[k] * 1000.0 * SECONDS_PER_YEAR), 0),
                           kept_wet_km2=round(float(capped[k]), 0),
                           evaporation_m3s=round(float(min(need[k], inflow[k])), 1)))
    record.update(
        floodplains_at=dict(hand_m=HAND_MAIN, count=int(count), area_km2=round(float(area[1:].sum()), 0),
                            kept_wet_km2=round(float(capped[1:].sum()), 0),
                            evaporation_m3s=round(float(np.minimum(need[1:], inflow[1:]).sum()), 1),
                            largest=listed),
        gcm=dict(surface=frame['surface'], water_cells=int((t21['mask'] == 0).sum()),
                 mask_land_fraction_at_most=round(cut, 4),
                 largest_added_water_fraction=round(float(fraction.max()), 4),
                 cells_with_added_water_over_1_percent=int((fraction > 0.01).sum()),
                 land_cells_that_would_turn_to_water=int(flips.sum()),
                 those_cells=[dict(row=int(i), col=int(j), lat=round(float(t21['lat'][i]), 2),
                                   lon=round(float(t21['lon'][j]), 2), land_fraction=round(float(t21['land'][i, j]), 4),
                                   added=round(float(fraction[i, j]), 4)) for i, j in np.argwhere(flips)]))
    record.update(
        schema=SCHEMA, producer=dict(files={f: RD.sha256(ROOT / f) for f in CODE}),
        inputs={str(p.relative_to(ROOT)): RD.sha256(p) for p in (GD.OUT, RD.OUT, TS.ST.OUT)},
        ignored_inputs={str(GD.GOVERNED_ROUTING): RD.sha256(GD.GOVERNED_ROUTING),
                        str(TS.PRODUCTS / frame['surface']): RD.sha256(TS.PRODUCTS / frame['surface'])},
        green=dict(path=str(GREEN), sha256=RD.sha256(GREEN)),
        settings=dict(q_river_m3s=Q_RIVER, width_coefficients=list(WIDTH), lunar_widening=round(LUNAR_WIDENING, 3),
                      riparian_belt_km=list(RIPARIAN_KM), floodplain_distance_km=DIST_KM, floodplain_heights_m=list(HAND_M)),
        evidence=('An estimate on products: the governed routing\'s discharges on the redrawn floor, the routing\'s own '
                  'open-water balance from the design climate, Earth\'s range of channel widths widened for lunar '
                  'gravity, and floodplains by height above the nearest river at 1.9 km cells. Not validated; the '
                  'belt and the floodplain heights are ranges chosen for scale.'),
        reading_rule=('Drylands are the cells where open water loses more to the air than the rain brings it; flows in '
                      'm3/s, areas in km2. channel gives the water surfaces of the dryland rivers and their loss at each '
                      'width coefficient; riparian_belt the plants drinking the bed\'s seepage; floodplains the dry land '
                      'within reach and height of a river at each height, as if all wet; floodplains_at lists the '
                      'connected floodplains at 5 m with the area each river can keep wet; gcm carries those, capped, '
                      'onto the run\'s T21 grid against the land fraction at which PlaSim\'s mask turns a cell to water.'),
        seconds=round(time.time() - started))
    return record


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    g = record['gcm']
    print(f'wrote {args.out} in {record["seconds"]} s; floodplains at 5 m {record["floodplains_at"]["area_km2"]:.0f} km2, '
          f'{record["floodplains_at"]["kept_wet_km2"]:.0f} km2 kept wet; largest added T21 water fraction '
          f'{g["largest_added_water_fraction"]}, cells that would turn to water {g["land_cells_that_would_turn_to_water"]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
