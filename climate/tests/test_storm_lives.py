"""Checks of the storm tracker (climate/crm/storm_lives.py) on cores laid out by hand."""
import numpy as np
import pytest

from climate.crm import storm_lives as sl

N = 64


def layout(*cores):
    """Label map and core list from cores given as (row, first column, last column), one row each."""
    labels = np.zeros((N, N), int)
    table = []
    for c, (row, i0, i1) in enumerate(cores, start=1):
        cols = np.arange(i0, i1 + 1) % N
        labels[row, cols] = c
        r, col = sl.periodic_centroid(np.full(cols.size, row), cols, np.ones(cols.size), N, N)
        table.append(dict(row=r, column=col))
    return labels, table


def test_cores_join_across_the_periodic_edges_by_sides_and_corners():
    mask = np.zeros((N, N), bool)
    mask[5, 0] = mask[5, N - 1] = True              # across x
    mask[0, 20] = mask[N - 1, 21] = True            # across y, corner to corner
    mask[0, 0] = mask[N - 1, N - 1] = False
    mask[30, 30] = mask[30, 32] = True              # two apart, not joined
    labels = sl.periodic_labels(mask)
    assert labels.max() == 4
    assert labels[5, 0] == labels[5, N - 1] and labels[0, 20] == labels[N - 1, 21]
    assert labels[30, 30] != labels[30, 32]
    corner = np.zeros((N, N), bool)
    corner[0, 0] = corner[N - 1, N - 1] = True
    assert sl.periodic_labels(corner).max() == 1


def test_open_edges_join_nothing_across_them():
    mask = np.zeros((N, N), bool)
    mask[5, 0] = mask[5, N - 1] = True
    labels = sl.periodic_labels(mask, periodic=False)
    assert labels.max() == 2
    wide = sl.widen(labels, periodic=False)
    assert (wide > 0).sum() == 2 * 6 and wide[5, 1] == labels[5, 0] and wide[5, N - 2] == labels[5, N - 1]
    assert sl.periodic_centroid(np.array([7, 7]), np.array([N - 1, 0]), np.ones(2), N, N, periodic=False) == (7.0, 31.5)
    assert sl.periodic_distance((0.0, 1.0), (0.0, N - 1.0), N, N, periodic=False) == pytest.approx(N - 2.0)


def test_widening_and_the_centre_wrap_around_the_box():
    labels = np.zeros((N, N), int)
    labels[0, 0] = 1
    wide = sl.widen(labels)
    assert (wide == 1).sum() == 9 and wide[N - 1, N - 1] == 1 and wide[1, 1] == 1
    r, c = sl.periodic_centroid(np.array([7, 7]), np.array([N - 1, 0]), np.ones(2), N, N)
    assert r == pytest.approx(7.0) and c == pytest.approx(N - 0.5)
    assert sl.periodic_distance((0.0, 1.0), (0.0, N - 1.0), N, N) == pytest.approx(2.0)


def test_tracks_carry_on_split_merge_and_jump_as_the_rules_say():
    steps = [layout((10, 10, 12), (40, 40, 45)),
             layout((10, 11, 13), (40, 40, 42), (40, 44, 45)),          # A moves; B splits, its larger part going on
             layout((10, 12, 14), (40, 40, 45)),                        # B's parts merge again
             layout((40, 47, 49), (20, 30, 31))]                        # A is gone; B jumps 5.5 columns; a new core C
    tracks = sl.follow([s[0] for s in steps], [s[1] for s in steps], match_columns=6.0)
    by_start = {(t['outputs'][0], t['labels'][0]): t for t in tracks}
    a, b = by_start[(0, 1)], by_start[(0, 2)]
    assert a['outputs'] == [0, 1, 2] and a['merged_into'] is None
    assert b['outputs'] == [0, 1, 2, 3]                                 # the overlap, then the jump within reach
    part = by_start[(1, 3)]
    assert part['split_from'] == tracks.index(b) and part['outputs'] == [1]
    assert part['merged_into'] == tracks.index(b)
    assert by_start[(3, 2)]['split_from'] is None and len(tracks) == 4
    # out of reach, the jump starts a new track
    assert len(sl.follow([s[0] for s in steps], [s[1] for s in steps], match_columns=4.0)) == 5


def test_flashes_go_to_the_core_at_the_nearest_output_or_beside_it():
    first, second = layout((10, 10, 12))[0], layout((10, 20, 22))[0]
    dx = 6000.0
    f = dict(time_s=np.array([100.0, 1000.0, 1000.0, 1000.0]),
             x_m=np.array([10.5, 21.5, 23.5, 40.5]) * dx, y_m=np.array([10.5, 10.5, 11.5, 10.5]) * dx)
    k, core = sl.assign_flashes(f, [0.0, 600.0, 1200.0], [first, first, second], dx, dx)
    assert k.tolist() == [0, 2, 2, 2]
    assert core.tolist() == [1, 1, 1, 0]                                # in it, beside it, nowhere near


def test_a_track_s_record_times_its_life_lead_and_lags():
    keys = ('area_km2', 'x_km', 'y_km', 'gh_kg', 'gh_zone_kg', 'ice_snow_zone_kg', 'w_max_m_s', 'updraft_km3',
            'top_km', 'rain_mm_h', 'positive_c', 'negative_c')
    w, gh, rain = [4.0, 9.0, 7.0, 5.0, 3.0], [1.0, 3.0, 8.0, 6.0, 2.0], [1.0, 2.0, 5.0, 9.0, 4.0]
    core_lists = [[]] + [[{**{key: 1.0 for key in keys}, 'w_max_m_s': w[i], 'gh_zone_kg': gh[i],
                           'rain_mm_h': rain[i], 'positive_c': 6.0 * i, 'negative_c': -float(i)}]
                         for i in range(5)] + [[]]
    hours = np.arange(7) / 6.0
    flashes = dict(count={(3, 1): 2, (4, 1): 1}, ground={(4, 1): 1},
                   times={(3, 1): [0.48, 0.52], (4, 1): [0.70]})
    track = dict(outputs=[1, 2, 3, 4, 5], labels=[1] * 5, split_from=None, merged_into=None)
    rec = sl.track_record(track, core_lists, hours, flashes, n_out=7, output_h=1.0 / 6.0)
    assert rec['whole'] and rec['life_h'] == pytest.approx(5.0 / 6.0)
    assert rec['flashes'] == 3 and rec['ground_strikes'] == 1
    assert rec['lead_h'] == pytest.approx(0.48 - (1.0 / 6.0 - 1.0 / 12.0))
    assert rec['peak_w_max_m_s_h'] == pytest.approx(2.0 / 6.0) and rec['peak_gh_zone_kg_h'] == pytest.approx(0.5)
    assert rec['graupel_after_updraft_h'] == pytest.approx(1.0 / 6.0, abs=1e-6)      # lags kept to 1e-6 h
    assert rec['flashes_after_graupel_h'] == 0.0
    assert rec['rain_after_flashes_h'] == pytest.approx(1.0 / 6.0, abs=1e-6)          # the rain peaks an output later
    assert rec['charge_lead_h'] == pytest.approx(0.48 - (3.0 / 6.0 - 1.0 / 12.0))      # 12 C at the third output
    assert rec['peak_flashes_per_min'] == pytest.approx(0.2) and rec['peak_negative_c'] == -4.0
    assert rec['after_last_flash_h'] == pytest.approx(5.0 / 6.0 + 1.0 / 12.0 - 0.70)
    assert not sl.track_record(dict(track, merged_into=0), core_lists, hours, flashes, 7, 1.0 / 6.0)['whole']
    out = sl.summarize([rec], flashes_total=4, flashes_assigned=3)
    assert out['flashing_tracks'] == 1 and out['flashes_in_cores'] == 3 and out['ground_strikes_in_cores'] == 1
    # ranks of graupel 1 3 5 4 2 against flashes 2 2 5 4 2 (three tied at none)
    assert out['flashes_spearman']['gh_zone_kg'] == pytest.approx(8.0 / np.sqrt(80.0))
    assert 'area_km2' not in out['flashes_spearman'] and out['flashes_spearman']['outputs'] == 5
