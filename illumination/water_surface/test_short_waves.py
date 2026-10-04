"""Checks of the short-wave spectrum against its paper, Cox and Munk's slopes and its gravity scaling."""
import numpy as np
import pytest

from illumination.water_surface import short_waves as S
from shared.constants import DATA, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

WU = DATA["model_closures"]["swan_41_51_wu"]


def wu_friction_velocity(u10):
    drag = (WU["low_speed_drag_coefficient"] if u10 < WU["transition_speed_m_s"]
            else WU["high_speed_intercept"] + WU["high_speed_slope_s_m"] * u10)
    return float(np.sqrt(drag) * u10)


def test_capillary_scales_follow_gravity():
    assert S.capillary_scales(STANDARD_GRAVITY) == pytest.approx((370.0, 0.2302), rel=1e-3)
    k_m, c_m = S.capillary_scales(MOON_SURFACE_GRAVITY)
    assert k_m == pytest.approx(370.0 * np.sqrt(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY))
    assert 2 * np.pi / k_m == pytest.approx(0.0417, abs=1e-4)
    k = np.geomspace(k_m / 10, k_m * 10, 20001)
    c = S.phase_speed(k, MOON_SURFACE_GRAVITY)
    assert k[np.argmin(c)] == pytest.approx(k_m, rel=1e-3) and c.min() == pytest.approx(c_m, rel=1e-6)


def test_the_short_wave_level_follows_the_two_regime_law():
    _, c_m = S.capillary_scales(STANDARD_GRAVITY)
    assert S.short_wave_level(c_m, STANDARD_GRAVITY) == pytest.approx(0.01)
    assert S.short_wave_level(2 * c_m, STANDARD_GRAVITY) == pytest.approx(0.01 * (1 + 3 * np.log(2)))
    assert S.short_wave_level(0.5 * c_m, STANDARD_GRAVITY) == pytest.approx(0.01 * (1 + np.log(0.5)))
    assert S.short_wave_level(0.3 * c_m, STANDARD_GRAVITY) == 0.0


@pytest.mark.parametrize("u10", [5.0, 7.0, 9.0, 11.0, 14.0])
def test_earth_slopes_match_cox_and_munk_for_a_clean_sea(u10):
    """Cox and Munk (1954): mss = 0.003 + 5.12e-3 W (scatter 0.004), W the wind at 12.5 m; the
    paper's figure 7a shows its spectrum within this from moderate to strong winds."""
    u_star = wu_friction_velocity(u10)
    roughness = 0.018 * u_star ** 2 / STANDARD_GRAVITY
    w = u10 * np.log(12.5 / roughness) / np.log(10 / roughness)
    slopes = S.slope_variances(u10, u_star, S.FULLY_DEVELOPED, STANDARD_GRAVITY)
    assert slopes["total"] == pytest.approx(0.003 + 5.12e-3 * w, abs=0.006)
    assert slopes["along"] > slopes["across"]


def test_strong_wind_anisotropy_matches_cox_and_munk():
    u10 = 14.0
    slopes = S.slope_variances(u10, wu_friction_velocity(u10), S.FULLY_DEVELOPED, STANDARD_GRAVITY)
    assert slopes["along"] - slopes["across"] == pytest.approx(3.16e-3 * u10 - 0.003 - 1.92e-3 * u10, abs=0.002)


def test_short_wave_slopes_depend_only_on_u_star_over_c_m():
    """Above a fixed fraction of k_m the short-wave slope variance is the same at any gravity."""
    results = []
    for gravity in (STANDARD_GRAVITY, MOON_SURFACE_GRAVITY):
        k_m, c_m = S.capillary_scales(gravity)
        u10 = 8.0 * np.sqrt(gravity / STANDARD_GRAVITY)
        results.append(S.slope_variances(u10, 1.6 * c_m, S.FULLY_DEVELOPED, gravity, k_from=0.2 * k_m)["short_wave_part"])
    assert results[0] == pytest.approx(results[1], rel=1e-9)


def test_integration_converges():
    a = S.slope_variances(6.0, 0.25, 1.2, MOON_SURFACE_GRAVITY, points=2000)
    b = S.slope_variances(6.0, 0.25, 1.2, MOON_SURFACE_GRAVITY, points=8000)
    assert a["total"] == pytest.approx(b["total"], rel=1e-4)
    assert a["along"] + a["across"] == pytest.approx(a["total"])


def test_inverse_wave_age_outside_the_fit_is_refused():
    with pytest.raises(ValueError):
        S.slope_variances(5.0, 0.2, 0.5, MOON_SURFACE_GRAVITY)


def test_neutral_wind_inverts_charnock_roughness():
    closure = DATA["model_closures"]["exoplasim_3_4_2"]
    u_star = 0.2
    u10 = S.neutral_wind(u_star, MOON_SURFACE_GRAVITY, closure["charnock_parameter"],
                         closure["sea_roughness_floor_m"], closure["von_karman_constant"])
    roughness = closure["charnock_parameter"] * u_star ** 2 / MOON_SURFACE_GRAVITY
    assert u10 == pytest.approx(u_star / closure["von_karman_constant"] * np.log(10 / roughness))


def test_short_wave_slopes_alone_match_the_full_spectrum_far_above_the_peak():
    gravity = MOON_SURFACE_GRAVITY
    k_m, c_m = S.capillary_scales(gravity)
    u_star, u10 = 0.25, 6.0
    full = S.slope_variances(u10, u_star, S.FULLY_DEVELOPED, gravity, k_from=k_m / 5)
    alone = S.short_wave_slopes(u_star, gravity, k_m / 5, peak_speed=float(S.phase_speed(
        gravity * S.FULLY_DEVELOPED ** 2 / u10 ** 2, gravity)))
    assert full["long_wave_part"] < 0.01 * full["total"]
    assert alone["total"] == pytest.approx(full["short_wave_part"], rel=1e-9)
    assert alone["along"] - alone["across"] == pytest.approx(
        (full["along"] - full["across"]) * full["short_wave_part"] / full["total"], rel=0.01)
    assert S.short_wave_slopes(0.3 * c_m, gravity, 6.0)["total"] == 0.0
