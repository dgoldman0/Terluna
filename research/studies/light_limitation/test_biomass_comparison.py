"""Check area conservation, GPP/NPP meaning, and reproducible scenario arithmetic."""
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
