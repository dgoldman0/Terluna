"""Checks of an open sea's slopes against the sea appearance study and of the neutral wind's inverse."""
import pytest

from illumination.water_surface import open_sea as O
from illumination.water_surface import short_waves as S
from shared.constants import DATA, MOON_SURFACE_GRAVITY

CLOSURES = DATA["model_closures"]["exoplasim_3_4_2"]
ARGS = (CLOSURES["charnock_parameter"], CLOSURES["sea_roughness_floor_m"], CLOSURES["von_karman_constant"])


def test_neutral_friction_velocity_inverts_the_neutral_wind():
    for u_star in (0.03, 0.094, 0.25):
        for height in (10.0, 100.0):
            wind = float(S.neutral_wind(u_star, MOON_SURFACE_GRAVITY, *ARGS, height=height))
            assert O.neutral_friction_velocity(wind, MOON_SURFACE_GRAVITY, *ARGS, height=height) == pytest.approx(u_star, rel=1e-6)


def test_open_sea_slopes_join_the_regimes_as_the_sea_appearance_study_does():
    # The study's western Oceanus Procellarum: median friction velocity 0.109378 m/s, short waves above its wave run's
    # 6.0008 rad/m with a median mean square slope of 0.00961 (research/studies/sea_appearance/results/sea_slopes.json).
    sea = O.open_sea_slopes(0.109378, MOON_SURFACE_GRAVITY, 6.0008, *ARGS)
    assert sea["short_waves"] == pytest.approx(0.00961, rel=0.01)
    assert sea["total"] == pytest.approx(sea["short_waves"] + sea["long_waves"])
    assert 0.01 < sea["long_waves"] < 0.025                     # the study's wave runs gave 0.012-0.016 at the coasts
    calm = O.open_sea_slopes(0.04, MOON_SURFACE_GRAVITY, 6.0008, *ARGS)
    assert calm["short_waves"] == 0.0 and calm["total"] < sea["total"]   # below the threshold only long waves tilt it
