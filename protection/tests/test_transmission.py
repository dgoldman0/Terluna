"""Checks on the swarm's ultraviolet transmission."""
import math

import numpy as np
import pytest

from protection.transmission import model as tm


def test_grun_flux_follows_the_published_form():
    # SPENVIS states the model as F = 3.15576e7 (F1 + F2 + F3) impacts per m^2 per year, m in grams
    m = 1e-6
    f1 = (2.2e3 * m ** 0.306 + 15.0) ** -4.38
    f2 = 1.3e-9 * (m + 1e11 * m ** 2 + 1e27 * m ** 4) ** -0.36
    f3 = 1.3e-16 * (m + 1e6 * m ** 2) ** -0.85
    assert tm.grun_flux(m) * 3.15576e7 == pytest.approx((f1 + f2 + f3) * 3.15576e7)
    masses = np.logspace(-15, 2, 50)
    assert np.all(np.diff(tm.grun_flux(masses)) < 0)


def test_holes_come_from_particles_of_tens_of_micrometres():
    rate = tm.holing_rate_per_year()
    assert 1e-7 < rate < 2e-7
    assert 30e-6 < tm.dominant_particle_diameter_m() < 100e-6
    assert tm.hole_fraction(20.0, 2.0, rate) == pytest.approx(rate * 4 * 10)


@pytest.mark.parametrize('overlap_m', [0.0, 3.0])
def test_seam_gaps_match_a_monte_carlo(overlap_m):
    sigma = 3.0
    rng = np.random.default_rng(1)
    separation = rng.normal(0.0, sigma, 2_000_000)
    mean_gap = np.maximum(separation - overlap_m, 0.0).mean()
    expected = 2.0 / tm.CELL_SIDE_M * mean_gap
    assert tm.seam_fraction(sigma, overlap_m, 0.0) == pytest.approx(expected, rel=0.02)


def test_sunlight_slant_eats_into_the_overlap():
    slant = 1000.0 * math.tan(tm.SUN_ANGULAR_RADIUS)
    assert tm.seam_fraction(1.0, 5.0, 1000.0) == pytest.approx(tm.seam_fraction(1.0, 5.0 - slant, 0.0))


def test_missing_cells_are_failures_times_time_to_cover():
    assert tm.missing_fraction(0.01, 7.0) == pytest.approx(0.01 * 7 / 365.25)


def test_levels_are_ordered_and_the_film_itself_is_opaque():
    rate = tm.holing_rate_per_year()
    t = {name: tm.level(p, rate)['transmission_by_hole_factor']['4'] for name, p in tm.LEVELS.items()}
    assert t['tight'] < t['standard'] < t['relaxed']
    assert tm.film_transmission_max() < 1e-40


def test_a_missing_cell_lights_a_wide_dim_patch():
    patch = tm.missing_cell_patch()
    assert patch['radius_km'] == pytest.approx(5.0 + 0.00465 * 78000.0)
    assert patch['share_of_full_sunlight'] == pytest.approx((1e8 / 78e6 ** 2) / (math.pi * 0.00465 ** 2))
    assert patch['share_of_full_sunlight'] < 1e-3
