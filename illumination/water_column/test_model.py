"""Checks of the seawater optics against their sources and against Earth's ocean-colour relations."""
import numpy as np
import pytest

from illumination.water_column import model as M

OC4 = (0.32814, -3.20725, 3.22969, -1.36769, -0.81739)    # NASA OC4 (R2018) band-ratio chlorophyll


def test_pure_water_is_clearest_in_the_violet_blue():
    w = np.arange(380.0, 701.0)
    a, b = M.seawater(w)
    assert 410 <= w[np.argmin(a)] <= 430 and a.min() == pytest.approx(0.0045, abs=0.001)
    assert b[w == 400.0][0] / b[w == 600.0][0] == pytest.approx((600 / 400) ** 4.32, rel=0.05)


def test_case_one_waters_follow_the_ocean_colour_band_ratios():
    """Without dissolved organic matter, which Earth's productive waters carry with their plankton and the model
    takes as a separate input, the ratios hold up to a few milligrams of chlorophyll per cubic metre."""
    bands = np.array([443.0, 490.0, 510.0, 555.0])
    for chl in (0.1, 0.3, 1.0, 3.0):
        iop = M.water(bands, chl)
        r = M.remote_sensing_reflectance(iop["a"], iop["bb"])
        x = np.log10(r[:3].max() / r[3])
        retrieved = 10 ** np.polyval(OC4[::-1], x)
        assert 0.4 < retrieved / chl < 2.5, chl


def test_clear_water_attenuation_matches_morel():
    for chl in (0.05, 0.1, 0.3):
        iop = M.water(np.array([490.0]), chl)
        kd = M.diffuse_attenuation(iop["a"], iop["bb"], 30.0)[0]
        assert kd == pytest.approx(0.0166 + 0.07242 * chl ** 0.68955, rel=0.25)


def test_the_hapke_inversion_returns_the_measured_reflectance():
    lam, refl = M.soil_reflectance("12001")
    albedo = M.single_scattering_albedo(refl[::40])
    assert M.reflectance_factor(albedo) == pytest.approx(refl[::40], rel=1e-9)


def test_grains_carry_into_water_between_transparent_and_opaque():
    assert M.albedo_in_water(1.0, 1.7, 1.34) == pytest.approx(1.0)
    se_water = ((1.7 / 1.34 - 1) / (1.7 / 1.34 + 1)) ** 2 + 0.05
    assert M.albedo_in_water(0.01, 1.7, 1.34) == pytest.approx(se_water)


def test_highland_fines_scatter_more_and_mare_fines_absorb_more():
    w = np.arange(400.0, 701.0, 10)
    mare, highland = M.fines(w, 1.0, "12001"), M.fines(w, 1.0, "62231")
    assert np.all(highland[0] < mare[0]) and np.all(highland[2] > mare[2])
    for a, b, bb in (mare, highland):
        assert np.all((bb / b > 0.01) & (bb / b < 0.035))
    assert mare[0][0] > mare[0][-1]          # space-weathered fines absorb most in the blue


def test_constituents_add():
    w = np.array([450.0, 550.0])
    total = M.water(w, 1.0, 0.1, 2.0, "62231")
    assert total["a"] == pytest.approx(sum(p[0] for p in total["parts"].values()))
    assert total["bb"] == pytest.approx(sum(p[2] for p in total["parts"].values()))
