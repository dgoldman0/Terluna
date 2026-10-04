"""Mass conservation and CM1 dry-air reading conventions."""
import numpy as np
import pytest

from .cloud_columns import dry_density, layer_overlap, water_path


def test_dry_density_recovers_independent_ideal_mixture():
    constants = dict(rd=300., rv=450., cp=1200., p00=100000.)
    temperature = np.array([[280., 240.]])
    pressure = np.array([[90000., 40000.]])
    vapor = np.array([[.01, .002]])
    theta = temperature * (constants['p00'] / pressure) ** .25
    rho, t = dry_density(pressure, theta, vapor, constants)
    assert t == pytest.approx(temperature)
    assert rho * (300 + 450 * vapor) * temperature == pytest.approx(pressure)


def test_partial_layers_conserve_species_mass():
    edges = [0, 100, 500, 2000]
    rho = np.array([[1.2, 1., .5], [1., .8, .4]])
    q = np.array([[.001, .002, .003], [0, .001, 0]])
    all_mass = water_path(rho, q, edges)
    assert all_mass == pytest.approx([.12 + .8 + 2.25, .32])
    pieces = sum(water_path(rho, q, edges, lo, hi) for lo, hi in ((0, 250), (250, 1000), (1000, 2000)))
    assert pieces == pytest.approx(all_mass)
    assert layer_overlap(edges, 250, 1000) == pytest.approx([0, 250, 500])


def test_invalid_columns_rejected():
    with pytest.raises(ValueError):
        water_path([[1, 1]], [[0, -1]], [0, 1, 2])
    with pytest.raises(ValueError):
        layer_overlap([0, 2, 1], 0, 1)
