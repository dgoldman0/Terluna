import numpy as np
from scipy.spatial.transform import Rotation

from shared import constants as K
from .active_formation import lattice, solar_frame
from .natural_pattern import mutual_visibility
from .parallel_shadow import source_illumination, minimum_clearance
from .test_cycling import static_sample


def test_dynamic_shadow_candidates_match_all_pairs_for_tilted_two_sided_tiles():
    offsets, _, _, _ = lattice(5, 8500., 4000.)
    q = offsets[:, [1, 2, 0]]+np.array([1e6, 2e6, 15e6])
    source = np.array([0., 0., K.AU])
    for angle in [0., .8, 1.565, 1.575, 2.2, np.pi]:
        frame = Rotation.from_rotvec([angle, .12, .4]).as_matrix()
        reduced = source_illumination(q, frame, source)
        full = source_illumination(q, frame, source, every=True)
        np.testing.assert_allclose(reduced, full, atol=1e-11)
        assert np.all((reduced >= 0) & (reduced <= 1))


def test_sun_facing_limit_agrees_with_the_previous_coupled_pattern():
    sample = static_sample(); solar = solar_frame(sample)
    offsets, neighbour, valid, _ = lattice(5, 8500., 4000.)
    state = np.c_[np.array([15e6, 0., 0.])+offsets@solar.T, np.zeros((25, 3))]
    old, sources = mutual_visibility(state, solar, sample, neighbour, valid)
    u, b, c = solar.T; frame = np.column_stack([b, -c, -u])
    new = np.array([source_illumination(state[:, :3], frame, source) for source in sources])
    np.testing.assert_allclose(new, old, atol=1e-8)


def test_physical_clearance_uses_the_full_square():
    state = np.zeros((2, 6)); state[1, :3] = [9950., 0., 80.]
    distance, pair, projected = minimum_clearance(state, np.eye(3))
    assert distance == 80.
    np.testing.assert_array_equal(pair, [0, 1])
