"""The drop-size tools conserve water and number, separate fixed flux from fixed concentration, and match analytic
encounter integrals."""
import numpy as np
import pytest

from atmosphere.rain import distributions as dist


def test_gamma_spectrum_carries_its_water_flux_and_number():
    d, width = dist.diameter_grid(bins=1024)
    speed = 2 + 700 * d
    spectrum = dist.gamma_flux(d, width, speed, cloud_base_rain_mm_h=10)
    flux = spectrum['flux_m2_s']
    whole = dist.surface_summary(d, flux, d, speed, np.ones(d.size, dtype=bool))
    slower = dist.surface_summary(d, flux, d, speed / 2, np.ones(d.size, dtype=bool))
    shrunk = dist.surface_summary(d, flux, d * 0.8, speed, np.ones(d.size, dtype=bool))
    assert whole['surface_rain_mm_h'] == pytest.approx(10, abs=1e-12)
    assert slower['surface_arrivals_m2_s'] == pytest.approx(whole['surface_arrivals_m2_s'], rel=1e-12)
    assert slower['surface_concentration_m3'] / whole['surface_concentration_m3'] == pytest.approx(2, rel=1e-12)
    assert shrunk['water_survival_fraction'] == pytest.approx(0.8**3, rel=1e-12)
    assert spectrum['realized_dm_mm'] == pytest.approx(1.5, rel=0.01)


def test_lost_classes_leave_the_ledger():
    d, width = dist.diameter_grid(0.1, 6.0, 64)
    speed = np.full(d.size, 1.0)
    flux = dist.gamma_flux(d, width, speed, dm_mm=1.5, shape_mu=3, cloud_base_rain_mm_h=5)['flux_m2_s']
    alive = d > 1e-3
    ground = dist.surface_summary(d, flux, np.where(alive, d, np.nan), np.where(alive, speed, np.nan), alive)
    water = np.pi * d**3 / 6 * flux
    assert ground['water_survival_fraction'] == pytest.approx(water[alive].sum() / water.sum(), rel=1e-12)
    assert ground['number_survival_fraction'] == pytest.approx(flux[alive].sum() / flux.sum(), rel=1e-12)


def test_two_classes_meet_as_the_analytic_column_says():
    z = [100, 50, 0]
    sizes = np.tile([0.001, 0.002], (3, 1))
    speeds = np.tile([2.0, 4.0], (3, 1))
    flux = np.array([100.0, 50.0])
    out = dist.encounters(z, sizes, speeds, speeds, flux)
    kernel = np.pi * 0.0015**2 * 2
    expected = [100 / 2 * kernel * (50 / 4), 100 / 4 * kernel * (100 / 2)]
    assert np.allclose(out['tau'], expected, rtol=1e-13)
    assert np.allclose(out['volume_exposure'], np.array(expected) * [8, 0.125], rtol=1e-13)
    tenth = dist.encounters(z, sizes, speeds, speeds, flux, collision_efficiency=0.1)
    assert np.allclose(tenth['tau'], np.array(expected) * 0.1, rtol=1e-13)


def test_arrival_flux_conserves_water_and_separates_flux_from_concentration():
    sizes = np.array([0.0007, 0.0015, 0.003])
    flux = dist.flux_from_impact_shares(rain_mm_h=10, diameter_m=sizes, impact_number_share=[0.6, 0.3, 0.1],
                                        area_m2=3.25)
    assert flux.rain_mm_h == pytest.approx(10, rel=1e-14)
    assert flux.water_m3_per_second * 3600 == pytest.approx(0.0325, rel=1e-14)
    speed = np.array([3.0, 5.0, 7.0])
    concentration = flux.concentration_m3(speed)
    assert np.allclose(flux.concentration_m3(speed / 2), 2 * concentration)
    back = dist.flux_from_concentration(diameter_m=sizes, concentration_m3=concentration, normal_speed_m_s=speed,
                                        area_m2=3.25)
    assert np.allclose(back.number_flux_m2_s, flux.number_flux_m2_s)
    half = dist.flux_from_concentration(diameter_m=sizes, concentration_m3=concentration, normal_speed_m_s=speed / 2,
                                        area_m2=3.25)
    assert half.rain_mm_h == pytest.approx(5, rel=1e-14)
    single = [dist.flux_from_impact_shares(rain_mm_h=10, diameter_m=[mm * 1e-3], impact_number_share=[1],
                                           area_m2=1).arrivals_per_second for mm in (1, 2, 3)]
    assert np.allclose(single, np.array([1, 1 / 8, 1 / 27]) * (10 / 3.6e6) / (np.pi / 6 * 1e-9))


def test_invalid_inputs_are_refused():
    sizes = np.array([0.0007, 0.0015, 0.003])
    flux = dist.flux_from_impact_shares(rain_mm_h=10, diameter_m=sizes, impact_number_share=[0.6, 0.3, 0.1],
                                        area_m2=1)
    with pytest.raises(ValueError):
        flux.concentration_m3([0, 1, 1])
    with pytest.raises(ValueError):
        dist.flux_from_impact_shares(rain_mm_h=1, diameter_m=[0.001], impact_number_share=[0.5], area_m2=1)
    with pytest.raises(ValueError):
        dist.gamma_flux([0.001], [0.0001], [-1.0])
    assert flux.diameter_m.flags.writeable is False
