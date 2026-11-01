"""Check area conservation, GPP/NPP meaning, the per-area split, and reproducible scenario arithmetic."""
import json

import pytest

from research.studies.light_limitation import biomass_comparison as B


def test_aquatic_partition_conserves_area_and_converts_units():
    case = dict(coast_share=0.25, aquatic_npp_over_gpp=0.5,
                oligotrophic_gpp_g_c_m2_day=1 / B.JULIAN_YEAR_DAYS,
                coast_npp_g_c_m2_year=2)
    result = B.aquatic(1e9, 1e9, case)
    assert result["area_km2"]["open_sea"] + result["area_km2"]["coasts"] == 1e9
    assert result["total_gt_c_year"] == pytest.approx(1.375)
    assert result["oligotrophic_npp_g_c_m2_year"] == pytest.approx(0.5)


def test_npp_cannot_exceed_gpp():
    case = dict(coast_share=0.1, aquatic_npp_over_gpp=1.1,
                oligotrophic_gpp_g_c_m2_day=0.2, coast_npp_g_c_m2_year=300)
    with pytest.raises(ValueError):
        B.aquatic(10, 10, case)


def test_per_area_split_recombines_and_matches_the_totals():
    product = B.results()
    area, split = product["area"], product["per_area"]
    land_share = split["moon_land_share"]
    assert land_share == pytest.approx(area["dry_land_km2"] / area["moon_km2"])
    earth = split["earth"]
    earth_recombined = earth["land_share"] * earth["land_g_c_m2_year"] + (1 - earth["land_share"]) * earth["ocean_g_c_m2_year"]
    assert earth_recombined == pytest.approx(earth["mean_g_c_m2_year"])
    for name, case in split["cases"].items():
        recombined = land_share * case["land_g_c_m2_year"] + (1 - land_share) * case["water_g_c_m2_year"]
        assert recombined == pytest.approx(case["mean_g_c_m2_year"])
        total = product["cases"][name]["surface_total_gt_c_year"]
        assert case["mean_g_c_m2_year"] == pytest.approx(total * 1e15 / (area["moon_km2"] * 1e6))
    # Earth's own rates on the Moon's land share already exceed Earth's mean, because the Moon has more land.
    assert split["earth_rates_at_moon_land_share_over_earth"] > 1


def test_product_reproduces_and_sky_is_not_double_counted():
    product = B.results()
    assert json.loads(B.OUT.read_text()) == product
    for case in product["cases"].values():
        assert case["surface_total_gt_c_year"] == pytest.approx(
            case["land_gt_c_year"] + case["waters"]["total_gt_c_year"])
        water = case["waters"]
        assert water["oligotrophic_npp_g_c_m2_year"] <= water["oligotrophic_gpp_g_c_m2_year"]
    light = product["light_only_sensitivity"]
    assert max(light["provisional_revised_301k"]["fraction_npp_loss_if_only_light_changes"]) < 0.01
    # Only files whose values are read are bound by hash; context documents are cited by path.
    assert set(product["producer"]["inputs"]) == set(B.INPUTS)
    assert "research/decisions.md" not in product["producer"]["inputs"]
