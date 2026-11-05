"""A PlaSim restart's climate carried onto a new surface: the turned Moon's first state, and the step a simulation of the
build years would take at each change of coastline or orientation (5 November 2026).

    python -m climate.gcm.carry STATE BASE OUT [--frame RECORD]
        # STATE the restart whose climate is carried, BASE a restart PlaSim wrote on the new surface (a fresh start's
        # first lunar day), OUT the carried restart; RECORD a turned surface's record (turned.py) when the ground turns

A PlaSim restart carries the land-sea mask and the orography with the climate (plasim.f90 reads dls from it, landmod and
the ocean and ice modules their own masks), so a run on a new surface either starts afresh or starts from a restart
whose surface is the new one. This builds the second kind from two restarts. Every record PlaSim derives from the
surface at a fresh start (the three masks, the grid and spectral orography and the ground's geopotential, the runoff's
velocities, and the surface climatologies, uniform on this Moon) comes from BASE, which PlaSim wrote on the new surface.
The climate's state comes from STATE: the atmosphere's spectral fields and its cloud as they stand, and every surface
field moved with the ground (CARRIED), each cell of the new surface taking the value of STATE's nearest cell of its own
kind, land or water, to the place the cell's centre stood before the turn (R^T n, R the record's rotation; the identity
without a frame). The atmosphere keeps its own longitudes, which stay tied to the Sun, and finds the new orography under
its surface pressure, as PlaSim's own fresh start, with its uniform surface pressure, does at its first step. The
clock, the accumulation counters and the random seed are STATE's, so the Sun stands where STATE's atmosphere has it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from climate.gcm import boundary as bd
from climate.gcm import levels as lv

SCHEMA = 'terluna.climate.gcm-carried-restart/1'
NLAT, NLON = 32, 64                # the T21 grid the runs use, north to south from 0 E, as the restart stores it
SHAPE = ('nlat', 'nlon', 'nlev', 'nrsp', 'nlsoil', 'nlev_oce')        # must agree between the two restarts
SURFACE = ('so', 'dls', 'xls', 'yls', 'doro', 'groundsg', 'dglacsg', 'dglac', 'darea', 'duroff', 'dvroff', 'dwmax',
           'dforest', 'dz0clim', 'dz0climo', 'dalbcl', 'dalbcl1', 'dalbcl2', 'dtcl', 'dwcl', 'xclicec', 'xcliced',
           'xclsst', 'yclsst')      # from BASE: what PlaSim derives from the surface at a fresh start
KEPT = ('nstep', 'naccuout', 'seed', 'naccuice', 'naccuo', 'nicec2d', 'naccuoce', 'naccua', 'fixedlon', 'zsolars',
        'sz', 'sd', 'st', 'sq', 'sr', 'sp', 'szm', 'sdm', 'stm', 'sqm', 'spm', 'dcc', 'dql', 'dqsat',
        'aast', 'aasqout', 'aasd', 'aasz', 'aadqo3', 'aadq', 'aammr', 'aanrho', 'aadt', 'aadql', 'aadcc')
                                    # from STATE as they stand: the clock, the atmosphere and its accumulators
CARRIED = ('drhs', 'dalb', 'dsalb1', 'dsalb2', 'dz0', 'dicec', 'diced', 'dwatc', 'drunoff', 'dust3', 'dt', 'dq',
           'dtsl', 'dtsm', 'dqs', 'driver', 'dsnowt', 'dsnowz', 'dsoilt', 'dts',
           'xts', 'xicec', 'xiced', 'xsnow', 'ysst', 'yiflux', 'yfldo',
           'cheat', 'cpme', 'croff', 'ctaux', 'ctauy', 'cust3', 'csnow',
           'xflxicea', 'xheata', 'xofluxa', 'xqmelta', 'xcfluxa', 'xcfluxra', 'xcfluxna', 'xsmelta', 'ximelta',
           'xtsfluxa', 'xfluxca', 'xscflxa', 'xcpmea', 'xcroffa', 'xstoia', 'xicecc',
           'yheata', 'yfssta', 'yifluxa', 'ydssta', 'yfldoa', 'yqhda',
           'cheata', 'cpmea', 'cprsa', 'croffa', 'ctauxa', 'ctauya', 'cust3a', 'cshfla', 'cshdta', 'clhfla', 'clhdta',
           'cswfla', 'clwfla',
           'aprl', 'aprc', 'aprs', 'aevap', 'ashfl', 'alhfl', 'aroff', 'asmelt', 'asndch', 'acc', 'assol', 'asthr',
           'atsol', 'atthr', 'ataux', 'atauy', 'atsolu', 'assolu', 'asthru', 'aqvi', 'atsa', 'ats0', 'atsama',
           'atsami', 'azmuz', 'asigrain', 'tempmax', 'tempmin', 'aaso', 'aasp', 'aadmld', 'aadwatc', 'aadsnow',
           'aadust3', 'aadtd5', 'aadls', 'aadz0', 'aadalb', 'aadsalb1', 'aadsalb2', 'aadtsoil', 'aadtd2', 'aadtd3',
           'aadtd4', 'aadicec', 'aadiced', 'aadforest', 'aadwmax', 'aadglac', 'aagroundoro', 'aaglacieroro')
                                    # from STATE, moved with the ground: the surface's state and its accumulators


def grid_centres(nlat=NLAT, nlon=NLON):
    """The T21 cells' centres, north to south from 0 E, as unit vectors (nlat * nlon, 3) and degrees."""
    lat, _ = bd.gaussian_latitudes(nlat)
    lon = 360.0 / nlon * np.arange(nlon)
    la, lo = np.meshgrid(lat, lon, indexing='ij')
    from climate.gcm.turned import unit
    return unit(la.ravel(), lo.ravel()), la.ravel(), lo.ravel()


def source_cells(old_land, new_land, rotation=None, nlat=NLAT, nlon=NLON):
    """For each cell of the new surface, the index of the old surface's cell its value comes from: the nearest old cell
    of its own kind (land or water) to where the cell's centre stood before the turn."""
    vec, _, _ = grid_centres(nlat, nlon)
    before = vec if rotation is None else vec @ np.asarray(rotation, float)          # R^T n, row by row
    old_land, new_land = np.asarray(old_land, bool), np.asarray(new_land, bool)
    index = np.empty(vec.shape[0], int)
    for kind in (True, False):
        pool = np.flatnonzero(old_land == kind)
        if not pool.size:
            raise ValueError(f"the old surface holds no {'land' if kind else 'water'} for the new one's")
        want = np.flatnonzero(new_land == kind)
        index[want] = pool[np.argmax(before[want] @ vec[pool].T, axis=1)]
    return index


def carry(state, base, rotation=None):
    """The carried restart's records from STATE's and BASE's (levels.read_records), with what moved."""
    s = {name.strip(): data for name, data in state}
    b = {name.strip(): data for name, data in base}
    names = [name.strip() for name, _ in state]
    if sorted(names) != sorted(b):
        raise ValueError(f'the restarts hold different records: {sorted(set(names) ^ set(b))}')
    unknown = sorted(set(names) - set(SHAPE) - set(SURFACE) - set(KEPT) - set(CARRIED))
    if unknown:
        raise ValueError(f'records with no rule: {unknown}')
    for name in SHAPE:
        if s[name] != b[name]:
            raise ValueError(f'{name} differs between the restarts')
    nhor = NLAT * NLON
    old_land = np.frombuffer(s['dls'], '<f8') > 0.5
    new_land = np.frombuffer(b['dls'], '<f8') > 0.5
    for mask in ('xls', 'yls'):                  # the three modules' masks agree in each restart
        for r, land in ((s, old_land), (b, new_land)):
            if not np.array_equal(np.frombuffer(r[mask], '<f8') > 0.5, land):
                raise ValueError(f'{mask} disagrees with dls')
    index = source_cells(old_land, new_land, rotation)
    out = []
    for name, _ in state:
        key = name.strip()
        if key in SURFACE or key in SHAPE:
            out.append((name, b[key]))
        elif key in KEPT:
            out.append((name, s[key]))
        else:
            field = np.frombuffer(s[key], '<f8')
            if field.size % nhor:
                raise ValueError(f'{key} is not a grid of {nhor} cells')
            levels = field.reshape(-1, nhor)
            out.append((name, np.ascontiguousarray(levels[:, index]).astype('<f8').tobytes()))
    vec, lat, lon = grid_centres()
    flipped = old_land != new_land
    same = index == np.arange(nhor)
    moved = dict(cells=nhor, land_to_water=int((old_land & ~new_land).sum()),
                 water_to_land=int((~old_land & new_land).sum()), from_itself=int(same.sum()),
                 mean_move_deg=float(np.degrees(np.arccos(np.clip((vec * vec[index]).sum(axis=1), -1, 1))).mean()),
                 largest_move_deg=float(np.degrees(np.arccos(np.clip((vec * vec[index]).sum(axis=1), -1, 1))).max()))
    return out, moved, index


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('state', type=Path, help='the restart whose climate is carried')
    parser.add_argument('base', type=Path, help='a restart PlaSim wrote on the new surface')
    parser.add_argument('out', type=Path, help='the carried restart')
    parser.add_argument('--frame', type=Path, help="a turned surface's record (turned.py) when the ground turns")
    parser.add_argument('--record', type=Path, help="where to write the carry's record (default: beside OUT)")
    args = parser.parse_args(argv)
    frame = json.loads(args.frame.read_text()) if args.frame else None
    rotation = np.asarray(frame['frame']['rotation_today_to_turned'], float) if frame else None
    records, moved, index = carry(lv.read_records(args.state), lv.read_records(args.base), rotation)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    lv.write_records(args.out, records)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    record = dict(
        schema=SCHEMA, producer=dict(domain='climate', files={
            f'gcm/{name}': hashlib.sha256((bd.HERE / name).read_bytes()).hexdigest()[:16]
            for name in ('carry.py', 'levels.py')}),
        inputs={str(args.state): sha(args.state), str(args.base): sha(args.base),
                **({str(args.frame): sha(args.frame)} if args.frame else {})},
        out={str(args.out): sha(args.out)}, moved=moved,
        rules=dict(surface_from_base=list(SURFACE), kept_from_state=list(KEPT), carried_with_the_ground=list(CARRIED)),
        evidence=('A construction: the climate of one PlaSim restart carried onto the surface of another, the surface '
                  'fields moved with the ground cell by cell to the nearest cell of their own kind; an initial state '
                  'for a run whose first years are set aside, not a result.'),
        reading_rule=('moved counts the new surface\'s cells that changed kind and those that take their own old '
                      'value, and the distance each cell\'s value travels, in degrees of arc, from the cell it comes '
                      'from to the cell itself on the grid'))
    path = args.record or args.out.with_name(args.out.name + '.carry.json')
    path.write_text(json.dumps(record, indent=1) + '\n')
    print(f"{args.out}: {moved['land_to_water']} cells land to water, {moved['water_to_land']} water to land; values "
          f"travel {moved['mean_move_deg']:.2f} deg on average; record {path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
