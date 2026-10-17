import numpy as np
import pytest

from climate.crm import gcm_profiles as g
from climate.crm.cm1_run import plasim_half_levels

SIGMA = np.array([0.0383, 0.1191, 0.2108, 0.3168, 0.4368, 0.5668, 0.6994, 0.8233, 0.9241, 0.9833])


def test_local_time_bins_match_the_ring_analysis_composites():
    assert g.bin_index([-180.0, -144.01, -144.0, 0.0, 179.99, 180.0]).tolist() == [0, 0, 1, 5, 9, 9]


def test_an_isothermal_dry_column_stands_at_its_scale_heights():
    gravity, t = 1.6242, 290.0
    ta = np.full((2, 3, SIGMA.size), t)
    z = g.interface_heights(ta, np.zeros_like(ta), SIGMA, gravity)
    edges = plasim_half_levels(SIGMA)[::-1][:g.LAYERS + 1]                   # the ground's interface first
    assert z.shape == (2, 3, g.LAYERS + 1)
    assert z[1, 2] == pytest.approx(g.RD * t / gravity * np.log(edges[0] / edges), rel=1e-12)
    assert 1500.0 < z[0, 0, 1] < 2500.0 and 5500.0 < z[0, 0, 2] < 7500.0      # the lowest tops, about 2 and 7 km


def test_layer_means_weigh_the_air_by_mass_over_each_span():
    zh = np.arange(50.0, 30000.0, 100.0)
    prs = np.repeat((1.2e5 - zh)[:, None], 2, axis=1)                        # equal mass per level
    spans = np.array([0.0, 2000.0, 7000.0, 14000.0, 24000.0])
    uniform = g.layer_means(np.full_like(prs, 280.0), prs, zh, spans)
    assert uniform.shape == (g.LAYERS, 2) and np.allclose(uniform, 280.0)
    linear = g.layer_means(np.repeat((300.0 - 0.001 * zh)[:, None], 2, axis=1), prs, zh, spans)
    assert linear[:, 0] == pytest.approx(300.0 - 0.001 * 0.5 * (spans[:-1] + spans[1:]), abs=0.06)


def test_contrasts_are_differences_in_height_above_the_ground():
    side = dict(layers_c=[[20.0, 15.0, 5.0, -5.0]] * g.BINS, air_c=[24.0] * g.BINS, ground_c=[27.0] * g.BINS)
    c = g.contrasts(side)
    assert c['air_minus_layer2_k'] == [9.0] * g.BINS and c['layer1_minus_layer2_k'] == [5.0] * g.BINS
    assert c['ground_minus_air_k'] == [3.0] * g.BINS
