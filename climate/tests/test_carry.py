"""Carrying a PlaSim restart's climate onto a new surface (carry.py): every record has a rule, the surface comes from the
base restart, the atmosphere stands as it was, and the surface's state moves with the ground to cells of its own kind;
on synthetic restarts with the real records' names and sizes, and on the design run's restart when present."""
from pathlib import Path

import numpy as np
import pytest

from climate.gcm import carry as cr
from climate.gcm import levels as lv
from climate.gcm.turned import unit

NHOR = cr.NLAT * cr.NLON
DESIGN = Path('/media/projectspace/terluna-research/gcm-runs/A28R_301k_design_l15_snow_s934/model/MOST_REST.00019')


def sizes():
    """The design run's records' names and lengths in values (grid levels, spectral sizes, scalars)."""
    spectral = {'sz': 506 * 15, 'sd': 506 * 15, 'st': 506 * 15, 'sq': 506 * 15, 'sr': 506 * 15, 'sp': 506, 'so': 506,
                'szm': 506 * 15, 'sdm': 506 * 15, 'stm': 506 * 15, 'sqm': 506 * 15, 'spm': 506}
    levels = {'dcc': 15, 'dql': 15, 'dqsat': 15, 'aast': 15, 'aasqout': 15, 'aasd': 15, 'aasz': 15, 'aadqo3': 15,
              'aadq': 16, 'aammr': 16, 'aanrho': 16, 'aadt': 16, 'aadql': 16, 'aadcc': 16, 'dsoilt': 5,
              **{k: 14 for k in ('dtcl', 'dwcl', 'dalbcl', 'dalbcl1', 'dalbcl2', 'xclicec', 'xcliced', 'xclsst',
                                 'yclsst')}}
    out = {}
    for name in cr.SHAPE + ('nstep', 'naccuout', 'naccuice', 'naccuo', 'nicec2d', 'naccuoce', 'naccua'):
        out[name] = ('<i4', 1)
    out['seed'] = ('<i4', 8)
    out['fixedlon'] = ('<f8', 1)
    for name in cr.SURFACE + cr.KEPT + cr.CARRIED:
        if name not in out:
            out[name] = ('<f8', spectral.get(name, NHOR * levels.get(name, 1)))
    return out


def restart(land, seed, step=7):
    """A synthetic restart: random values everywhere, the three masks from `land`, scalars equal but the clock."""
    rng = np.random.default_rng(seed)
    records = []
    for name, (kind, n) in sizes().items():
        if kind == '<i4':
            v = np.full(n, step if name == 'nstep' else 32 if name == 'nlat' else 3, '<i4')
        elif name in ('dls', 'xls', 'yls'):
            v = np.asarray(land, '<f8')
        else:
            v = rng.random(n).astype('<f8')
        records.append((name.ljust(16), v.tobytes()))
    return records


def field(records, name):
    return np.frombuffer(dict((k.strip(), v) for k, v in records)[name], '<f8' if name not in sizes() or
                         sizes()[name][0] == '<f8' else '<i4')


def test_every_record_of_the_design_run_has_one_rule():
    groups = [set(cr.SHAPE), set(cr.SURFACE), set(cr.KEPT), set(cr.CARRIED)]
    assert sum(len(g) for g in groups) == len(set().union(*groups))
    if DESIGN.exists():
        names = {name.strip() for name, _ in lv.read_records(DESIGN)}
        assert names == set().union(*groups)


def test_on_the_same_surface_without_a_turn_the_climate_is_carried_unchanged():
    land = np.random.default_rng(1).random(NHOR) < 0.6
    state, base = restart(land, 2, step=99), restart(land, 3)
    out, moved, index = cr.carry(state, base)
    assert np.array_equal(index, np.arange(NHOR)) and moved['from_itself'] == NHOR
    for name in cr.SURFACE:
        assert np.array_equal(field(out, name), field(base, name)), name
    for name in cr.KEPT + cr.CARRIED:
        assert np.array_equal(field(out, name), field(state, name)), name
    assert [k for k, _ in out] == [k for k, _ in state]


def test_the_surface_state_moves_with_the_ground_to_cells_of_its_own_kind():
    lat = np.repeat(np.arange(cr.NLAT), cr.NLON)
    lon = np.tile(np.arange(cr.NLON), cr.NLAT)
    old = (lon < 32)                                    # land on the first half of the longitudes
    a = np.radians(45.0)                                # a turn of eight cells about the axis
    turn = np.array([[np.cos(a), np.sin(a), 0.0], [-np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
    new = np.roll(old.reshape(cr.NLAT, cr.NLON), -8, axis=1).ravel()
    new[(lat == 10) & (lon == 30)] = True               # one cell that the turn alone would leave as water
    state, base = restart(old, 4), restart(new, 5)
    out, moved, index = cr.carry(state, base, turn)
    assert np.array_equal(field(out, 'dls'), new.astype(float))
    assert np.array_equal((field(state, 'dls') > 0.5)[index], new)            # every value from its own kind
    plain = (lat * cr.NLON + (lon + 8) % cr.NLON)
    shifted = new == np.roll(old.reshape(cr.NLAT, cr.NLON), -8, axis=1).ravel()
    assert np.array_equal(index[shifted], plain[shifted])
    soil = field(state, 'dsoilt').reshape(5, NHOR)
    assert np.array_equal(field(out, 'dsoilt').reshape(5, NHOR), soil[:, index])
    assert np.array_equal(field(out, 'st'), field(state, 'st'))               # the atmosphere stays as it was
    assert np.array_equal(field(out, 'so'), field(base, 'so'))                # the dynamics see the new orography
    assert moved['water_to_land'] >= 1


def test_records_without_a_rule_or_differing_shapes_are_refused():
    land = np.ones(NHOR, bool)
    land[:100] = False
    state, base = restart(land, 6), restart(land, 7)
    with pytest.raises(ValueError, match='no rule'):
        cr.carry(state + [('newfield'.ljust(16), np.zeros(NHOR).tobytes())],
                 base + [('newfield'.ljust(16), np.zeros(NHOR).tobytes())])
    other = [(k, (np.array([64], '<i4').tobytes() if k.strip() == 'nlat' else v)) for k, v in base]
    with pytest.raises(ValueError, match='nlat'):
        cr.carry(state, other)
    broken = [(k, (np.zeros(NHOR).tobytes() if k.strip() == 'yls' else v)) for k, v in base]
    with pytest.raises(ValueError, match='yls'):
        cr.carry(state, broken)


def test_the_design_runs_restart_carried_onto_itself_is_unchanged():
    if not DESIGN.exists():
        pytest.skip('the design run is not on the research drive')
    records = lv.read_records(DESIGN)
    out, moved, _ = cr.carry(records, records)
    assert out == records and moved['from_itself'] == NHOR
