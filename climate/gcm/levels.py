"""PlaSim's vertical levels: its own ten, a set of fourteen with more levels near the ground, and a set of fifteen with a
thin bottom layer (2026-10-26).

    climate/gcm/.venv/bin/python -m climate.gcm.levels show l15                  # the set, its heights and layers
    climate/gcm/.venv/bin/python -m climate.gcm.levels restart SRC DST --to l15  # a restart carried onto the set

PlaSim puts each full level midway in sigma between the half levels that bound its layer (plasim.f90, initpm), from 0
at the top to 1 at the ground, so the half levels fix the full levels. Its own grid, ExoPlaSim's vtype 0, sets ten
half levels by a formula (plasim_half); on the Moon, at the 51 km scale height of the design air, its full levels stand
near 0.9, 4.0, 9.9, 18.2, 29, 42, 59, 79, 108 and 166 km. CM1's settled sea box carries a boundary layer 6 km deep and
night cloud from its level of free convection at 6.2 km to tops at 7-8 km (climate/crm/README.md), all inside the
single layer between 6.4 and 13.7 km, so the set 'l14' adds levels there.

How the set was chosen. The six full levels from 29 km up keep their half levels, so the upper air, the top two levels
the exobase check reads, and the three layers the cold cloud's factors describe stay exactly as they were; the lowest
layer, the ground to 1.73 km, stays too, with its level at 0.86 km where the surface fluxes and the 10 m wind are taken.
Below the half level at 23.3 km six new half levels stand at 15.5, 11.5, 9.0, 7.0, 5.0 and 3.1 km (sigma rounded to four
places), which puts full levels at 19.2, 13.5, 10.2, 8.0, 6.0, 4.0 and 2.4 km. The layers thin smoothly from 7.8 km
to about 2 km below 9 km, never by more than half from one layer to the next in sigma. The present level at 18.2 km
cannot be kept with one at 13.5 km: with the half level at 23.3 km fixed, a full level at 18.2 km puts the next half
level at 13.7 km, and the level below it then falls near 12 km under a layer three times thicker than its own. So the
18.2 km level moves up to 19.2 km and the 9.9 km level to 10.2 km, and the new levels land within 0.04 km of 2.4, 6, 8
and 13.5 km. Heights here are -51 km times ln sigma; the run's own temperatures set its levels' heights.

The set 'l15' keeps every half level of 'l14' and adds one at sigma 0.9883, 0.60 km, splitting the lowest layer, the
ground to 1.73 km, in two. Its lowest full level stands midway in sigma between that half level and the ground, at
0.99415 and 0.30 km, and the level above it midway in the layer from 0.60 to 1.73 km, at 0.97745 and 1.16 km. The
offline mass-flux scheme (massflux.md) found that the 1.73 km layer averages away the near-ground air CM1's updraughts
carry: it holds 1.4-2.3 kJ/kg less moist static energy than CM1's air at 50 m. The new bottom layer holds the lowest
0.6 km, close to the lowest 0.5 km over which CM1's vapour is compared. PlaSim takes its surface fluxes and 10 m wind at
the lowest level (fluxmod.f90, surflx), at the height its hydrostatic relation gives from that level's own temperature,
about 0.31 km in the design climate's 300 K air, so they move down with the level. With one half level added, 0.3 km is
near the lowest the bottom level can stand while the bottom layer keeps at least half the mass of the layer above it,
the limit the 14-level set kept between neighbours: the limit is a half level at sigma 0.9889 with a lowest level at
0.28 km, and the set's ratio is 0.54. The bottom layer's 0.0117 of sigma holds about 880 kg/m2 of air at the design
surface pressure, a third of the 14-level bottom layer's mass and two and a half times that of PlaSim's lowest layer on
Earth.

PlaSim presets Rayleigh friction on the top two levels only for its ten-level grid (20 and 100 days, plasim.f90,
readnl), so a run on another set writes the same friction for its top two levels, which are the same levels
(namelist_entries). Its reference temperature for the semi-implicit scheme is 250 K at every level whatever the grid.

Carrying a restart onto another set (remap_matrix, carry_restart). Sigma counts the air's mass, a layer holding
ps dsigma / g, so the regridding is a conservative remap in sigma: each old layer's profile is reconstructed linearly
about its full level with the centred slope between its neighbours' values (one-sided at the top and the ground), and
each new layer takes the reconstruction's mean over its sigma range, the parts in each old layer weighted by their
mass. A layer the two grids share keeps its values bit for bit. The surface pressure, and with it every column's mass,
is left as it is, and the remap keeps every column's integrals over sigma: its water, its sum of temperature times
mass, and its momentum. It is one linear operator in the vertical, the same at every point, so it acts on the
spectral coefficients of vorticity, divergence, temperature and humidity (both time levels) directly, and on the
gridpoint fields held by level (cloud cover, cloud water, saturation humidity and the accumulated means), the cloud
cover clipped to 0-1 and the water to 0 and above.
"""
from __future__ import annotations

import argparse
import hashlib
import struct
import sys
from pathlib import Path

import numpy as np

SCALE_HEIGHT_KM = 51.0          # the design air's, by which the levels' heights here are counted
DEFAULT = '10'                  # PlaSim's own grid (vtype 0), which every run before 2026-10-26 had


def plasim_half(levels: int = 10) -> np.ndarray:
    """PlaSim's own half levels (vtype 0, initpm's last branch), the top layer's bottom first, ending at 1."""
    s = np.arange(1, levels + 1) / levels
    return 0.75 * s + 1.75 * s ** 3 - 1.5 * s ** 4


def full_levels(half) -> np.ndarray:
    """The full levels PlaSim puts midway in sigma between the half levels (0 above the first)."""
    half = np.asarray(half, float)
    return 0.5 * (np.r_[0.0, half[:-1]] + half)


_TEN = plasim_half(10)
# The sets of half levels a run may take, the top layer's bottom first and the ground's 1 last. 'l14' keeps PlaSim's
# own six upper half levels and its lowest one as the formula gives them, so the layers the grids share are the same
# numbers, and puts six new half levels at 15.5, 11.5, 9.0, 7.0, 5.0 and 3.1 km (sigma to four places). 'l15' is 'l14'
# with one more half level, at 0.60 km (sigma 0.9883), splitting its lowest layer.
_L14 = np.r_[_TEN[:6], 0.7379, 0.7981, 0.8382, 0.8717, 0.9066, 0.9410, _TEN[8], 1.0]
LEVEL_SETS = {
    DEFAULT: tuple(float(x) for x in _TEN),
    'l14': tuple(float(x) for x in _L14),
    'l15': tuple(float(x) for x in np.r_[_L14[:-1], 0.9883, 1.0]),
}
# The levels near which each set was asked to put full levels (km at the 51 km scale height), for its tests.
L14_TARGETS_KM = (13.5, 8.0, 6.0, 2.4)
L15_LOWEST_KM = 0.3
# PlaSim's Rayleigh friction on its top two levels, which it presets only on its ten-level grid (days).
TOP_FRICTION_DAYS = (20.0, 100.0)
DAY_S = 86400.0


def half_levels(name: str = DEFAULT) -> np.ndarray:
    if name not in LEVEL_SETS:
        raise ValueError(f'unknown level set {name!r}; the sets are {sorted(LEVEL_SETS)}')
    half = np.array(LEVEL_SETS[name], float)
    if half[-1] != 1.0 or np.any(np.diff(np.r_[0.0, half]) <= 0.0):
        raise ValueError(f'level set {name!r} does not rise from 0 to 1')
    return half


def heights_km(sigma, scale_height_km: float = SCALE_HEIGHT_KM) -> np.ndarray:
    """Heights of sigma levels at a constant scale height, km."""
    return -scale_height_km * np.log(np.asarray(sigma, float))


def namelist_entries(name: str) -> list:
    """The plasim_namelist entries (name, text) a run on a level set other than PlaSim's own writes: its half levels
    (NEQSIG -1 reads them from SIGH) and PlaSim's friction on its top two levels, in seconds."""
    if name == DEFAULT:
        return []
    half = half_levels(name)
    friction = [d * DAY_S for d in TOP_FRICTION_DAYS] + [0.0] * (half.size - len(TOP_FRICTION_DAYS))
    return [('NEQSIG', '-1'), ('SIGH', ','.join(repr(float(x)) for x in half)),
            ('TFRC', ','.join(repr(float(x)) for x in friction))]


def shared_layers(half_from, half_to) -> dict:
    """The layers two grids share, new index to old: the same two half levels (0 above the first)."""
    a, b = np.r_[0.0, np.asarray(half_from, float)], np.r_[0.0, np.asarray(half_to, float)]
    out = {}
    for k in range(b.size - 1):
        for j in range(a.size - 1):
            if a[j] == b[k] and a[j + 1] == b[k + 1]:
                out[k] = j
    return out


def remap_matrix(half_from, half_to) -> np.ndarray:
    """The conservative remap from one grid's layer values to another's (new by old, the top first): a linear
    reconstruction in each old layer about its full level with the centred slope of its neighbours, averaged over each
    new layer by mass (sigma). Layers the grids share are copied; every column's integral over sigma is kept."""
    a = np.r_[0.0, np.asarray(half_from, float)]
    b = np.r_[0.0, np.asarray(half_to, float)]
    if a[-1] != 1.0 or b[-1] != 1.0:
        raise ValueError('both grids must end at the ground, sigma 1')
    n, m = a.size - 1, b.size - 1
    centre = 0.5 * (a[:-1] + a[1:])
    slope = np.zeros((n, n))                      # slope[j] @ x: the old layer j's slope in sigma
    for j in range(n):
        lo, hi = max(j - 1, 0), min(j + 1, n - 1)
        slope[j, hi] += 1.0 / (centre[hi] - centre[lo])
        slope[j, lo] -= 1.0 / (centre[hi] - centre[lo])
    w = np.zeros((m, n))
    shared = shared_layers(half_from, half_to)
    for k in range(m):
        if k in shared:
            w[k, shared[k]] = 1.0
            continue
        for j in range(n):
            lo, hi = max(a[j], b[k]), min(a[j + 1], b[k + 1])
            if hi <= lo:
                continue
            part = (hi - lo) / (b[k + 1] - b[k])
            w[k, j] += part
            w[k] += part * slope[j] * (0.5 * (lo + hi) - centre[j])
    return w


def remap(values, half_from, half_to, axis: int = 0) -> np.ndarray:
    """Layer values (any shape, the levels along axis) carried from one grid to another by remap_matrix; layers the
    grids share are copied exactly."""
    w = remap_matrix(half_from, half_to)
    x = np.moveaxis(np.asarray(values, float), axis, 0)
    out = np.tensordot(w, x, axes=(1, 0))
    for k, j in shared_layers(half_from, half_to).items():
        out[k] = x[j]
    return np.moveaxis(out, 0, axis)


# The restart's fields held by level (plasim.f90, the restart's write and read): spectral at each level (both time
# levels), gridpoint at each level, and gridpoint at each level and the ground (NLEP), whose ground value is kept.
SPECTRAL_BY_LEVEL = ('sz', 'sd', 'st', 'sq', 'sr', 'szm', 'sdm', 'stm', 'sqm')
GRID_BY_LEVEL = ('dcc', 'dql', 'dqsat', 'aast', 'aasqout', 'aasd', 'aasz', 'aadqo3')
GRID_BY_LEVEL_AND_GROUND = ('aadq', 'aammr', 'aanrho', 'aadt', 'aadql', 'aadcc')
CLIPPED = {'dcc': (0.0, 1.0), 'dql': (0.0, None), 'dqsat': (0.0, None)}
# Records held by something other than the air's levels: the 14 months of the surface's climatologies and the soil's
# layers.
NOT_BY_LEVEL = ('dtcl', 'dwcl', 'dalbcl', 'dalbcl1', 'dalbcl2', 'xclicec', 'xcliced', 'xclsst', 'yclsst', 'dsoilt')


def read_records(path) -> list:
    """A PlaSim restart's records (Fortran sequential, unformatted, little-endian 4-byte markers): pairs of a
    16-character name and its data, as (name, bytes)."""
    data = Path(path).read_bytes()
    raw, i = [], 0
    while i < len(data):
        size = struct.unpack('<i', data[i:i + 4])[0]
        raw.append(data[i + 4:i + 4 + size])
        if struct.unpack('<i', data[i + 4 + size:i + 8 + size])[0] != size:
            raise ValueError(f'{path}: a record at byte {i} does not close as it opens')
        i += size + 8
    if len(raw) % 2:
        raise ValueError(f'{path}: a name without its data')
    return [(raw[k].decode('ascii'), raw[k + 1]) for k in range(0, len(raw), 2)]


def write_records(path, records):
    out = bytearray()
    for name, data in records:
        for rec in (name.ljust(16)[:16].encode('ascii'), data):
            out += struct.pack('<i', len(rec)) + rec + struct.pack('<i', len(rec))
    tmp = Path(str(path) + '.part')
    tmp.write_bytes(bytes(out))
    tmp.replace(path)


def carry_restart(src, dst, half_from, half_to) -> dict:
    """Write src's restart onto another grid at dst (remap_matrix), every other record as it was; returns what moved
    and the column integrals before and after."""
    half_from, half_to = np.asarray(half_from, float), np.asarray(half_to, float)
    n, m = half_from.size, half_to.size
    records = read_records(src)
    names = [name.strip() for name, _ in records]
    nlev = struct.unpack('<i', records[names.index('nlev')][1])[0]
    if nlev != n:
        raise ValueError(f'{src} holds {nlev} levels, the grid it is carried from {n}')
    dsig_from, dsig_to = np.diff(np.r_[0.0, half_from]), np.diff(np.r_[0.0, half_to])
    out, moved = [], {}
    for name, data in records:
        key = name.strip()
        values = np.frombuffer(data, '<f8') if len(data) % 8 == 0 else None
        if key == 'nlev':
            out.append((name, struct.pack('<i', m)))
            continue
        if key in SPECTRAL_BY_LEVEL or key in GRID_BY_LEVEL or key in GRID_BY_LEVEL_AND_GROUND:
            levels = n + (key in GRID_BY_LEVEL_AND_GROUND)
            if values is None or values.size % levels:
                raise ValueError(f'{src}: {key} does not hold {levels} levels')
            x = values.reshape(levels, -1)               # Fortran's (point, level) order: one level after another
            new = remap(x[:n], half_from, half_to)
            if key in GRID_BY_LEVEL_AND_GROUND:
                new = np.concatenate([new, x[n:]])
            if key in CLIPPED:
                lo, hi = CLIPPED[key]
                new = np.clip(new, lo, hi)
            before = float(np.abs(dsig_from @ x[:n]).sum())
            after = float(np.abs(dsig_to @ new[:m]).sum())
            moved[key] = dict(points=int(x.shape[1]), levels=[levels, int(new.shape[0])],
                              column_integral_change=float(np.abs(dsig_to @ new[:m] - dsig_from @ x[:n]).max()),
                              column_integral_scale=max(before, after) / x.shape[1])
            out.append((name, np.ascontiguousarray(new, '<f8').tobytes()))
            continue
        if values is not None and values.size and any(values.size == k * p for k in (n, n + 1) for p in (506, 2048)) \
                and key not in NOT_BY_LEVEL:
            raise ValueError(f'{src}: {key} is sized like a field by level, which this carrier does not know')
        out.append((name, data))
    write_records(dst, out)
    return dict(source=str(src), source_sha256=hashlib.sha256(Path(src).read_bytes()).hexdigest()[:16],
                written=str(dst), written_sha256=hashlib.sha256(Path(dst).read_bytes()).hexdigest()[:16],
                levels=[n, m], fields=moved)


def describe(name: str) -> str:
    half = half_levels(name)
    full = full_levels(half)
    rows = [f'level set {name}: {half.size} levels (heights at the {SCALE_HEIGHT_KM:g} km scale height)',
            '  level  full sigma  km      layer from-to sigma      km           dsigma']
    top = np.r_[0.0, half[:-1]]
    for k in range(half.size):
        z_top = heights_km(top[k]) if top[k] > 0 else float('inf')
        rows.append(f'  {k + 1:5d}  {full[k]:.5f}  {heights_km(full[k]):6.2f}  {top[k]:.4f}-{half[k]:.4f}  '
                    f'{z_top:6.2f}-{heights_km(half[k]):5.2f}  {half[k] - top[k]:.4f}')
    return '\n'.join(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='action', required=True)
    show = sub.add_parser('show')
    show.add_argument('name', nargs='?', default='l14')
    carry = sub.add_parser('restart')
    carry.add_argument('src')
    carry.add_argument('dst')
    carry.add_argument('--to', default='l14')
    carry.add_argument('--from', dest='source_set', default=DEFAULT)
    args = parser.parse_args(argv)
    if args.action == 'show':
        print(describe(args.name))
        return 0
    if Path(args.dst).exists():
        print(f'{args.dst} exists; not overwriting it', file=sys.stderr)
        return 1
    report = carry_restart(args.src, args.dst, half_levels(args.source_set), half_levels(args.to))
    for key, r in report['fields'].items():
        print(f"  {key:8s} {r['levels'][0]} -> {r['levels'][1]} levels, {r['points']} values each; largest change "
              f"of a column integral {r['column_integral_change']:.3e} (typical size {r['column_integral_scale']:.3e})")
    print(f"wrote {report['written']} ({report['written_sha256']}) from {report['source']} ({report['source_sha256']})")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
