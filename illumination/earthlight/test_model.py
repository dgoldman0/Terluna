"""Checks of the earthlight model: its source table, phase curves and the visual calibration."""
import numpy as np
import pytest
from scipy.integrate import quad

from illumination.earthlight import model as E
from shared.constants import EARTH_MOON_DISTANCE


def test_the_table_spans_the_visible_and_near_infrared():
    c = E.coefficients()
    assert len(c["wavelength_nm"]) == 733
    assert c["wavelength_nm"][0] == pytest.approx(300.3) and c["wavelength_nm"][-1] == pytest.approx(2487.562)
    assert np.all((c["width"] >= 0) & (c["width"] <= 1))


def test_phase_curves_are_one_at_full_phase_and_fall():
    c = E.coefficients()
    assert E.phase_curve(0.0, c["width"]) == pytest.approx(np.ones_like(c["width"]))
    phases = np.arange(0, 181, 10.0)
    assert np.all(np.diff(E.visual_phase_curve(phases)) < 0)
    assert E.visual_phase_curve(0.0) == pytest.approx(1.0)
    assert E.visual_phase_curve(90.0) == pytest.approx(0.2576, abs=1e-4)


def test_the_visual_fit_gives_robinsons_phase_integral_and_spherical_albedo():
    """Robinson et al. (2025): phase integral 1.22 and spherical albedo 0.294 from their full model; the analytic
    fit, which levels off near full phase, with its Lambert-shaped tail past their data, gives a little more."""
    q = 2 * quad(lambda a: float(E.visual_phase_curve(np.degrees(a))) * np.sin(a), 0, np.pi,
                 points=[np.radians(144.0)])[0]
    assert q == pytest.approx(1.32, abs=0.01)
    assert q * E.EARTH_VISUAL_PHASE_NORMALISATION == pytest.approx(0.294, rel=0.05)


def test_calibration_reproduces_the_visual_brightness_and_keeps_the_spectral_shape():
    bands = np.array([[400.0 + 10 * k, 410.0 + 10 * k] for k in range(40)])
    solar = np.full(len(bands), 18.0)                     # a flat sun of 1.8 W m-2 nm-1
    for phase in (0.0, 45.0, 100.0):
        e = E.calibrated_irradiance(bands, solar, phase, EARTH_MOON_DISTANCE)
        visual = (bands[:, 0] >= 400) & (bands[:, 1] <= 700)
        level = e[visual].sum() / solar[visual].sum() / E.solid_angle(EARTH_MOON_DISTANCE)
        assert level == pytest.approx(E.visual_intensity(phase), rel=1e-9)
        shape = E.band_radiance(bands, phase)
        assert e / e[0] == pytest.approx(shape / shape[0], rel=1e-9)


def test_the_earth_is_bluer_than_the_sun():
    blue, red = E.band_radiance([[440.0, 460.0], [640.0, 660.0]], 0.0)
    assert blue > red


def test_beyond_the_observed_phases_the_earth_fades_to_dark_without_a_step():
    bands = np.array([[500.0, 510.0], [600.0, 610.0]])
    solar = np.array([19.0, 17.0])
    at, past = (E.calibrated_irradiance(bands, solar, a, EARTH_MOON_DISTANCE) for a in (144.0, 144.0001))
    assert past == pytest.approx(at, rel=1e-4)
    assert E.calibrated_irradiance(bands, solar, 180.0, EARTH_MOON_DISTANCE) == pytest.approx([0.0, 0.0], abs=1e-15)
    assert E.visual_phase_curve(143.9999) == pytest.approx(E.visual_phase_curve(144.0001), rel=1e-4)
    assert E.visual_phase_curve(180.0) == pytest.approx(0.0, abs=1e-15)
