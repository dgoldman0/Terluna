"""Checks for the independent Python reader and saved post-export audit."""
import copy
import json

import numpy as np
import pytest

from illumination.calendar import validate as V


def flat_product():
    return dict(schema="terluna.illumination.calendar-transfer/1",
                constants=dict(earth_radius_m=1.0, earth_reference_distance_m=100.0),
                direct_elevation_deg=[-90.0,0.0,90.0], diffuse_elevation_deg=[-90.0,0.0,90.0],
                phase_deg=[0.0,90.0,180.0],
                solar=dict(direct_lux=[2.0,2.0,2.0], diffuse_lux=[3.0,3.0,3.0]),
                earth=dict(direct_lux=[[7.0]*3]*3, diffuse_lux=[[11.0]*3]*3))


def flat_geometry():
    return dict(sun_elevation_deg=0.0, sun_radius_deg=0.267, sun_distance_factor=1.0,
                earth_elevation_deg=0.0, earth_radius_deg=1.0, earth_phase_deg=90.0,
                earth_bright_limb_rad=0.0, earth_distance_m=100.0,
                earth_sunlight_factor=1.0, source_separation_deg=90.0)


def test_disk_weights_do_not_apply_the_earth_phase_law_twice():
    reader=V.Response(flat_product())
    for phase in (0, 60, 90, 150, 179):
        result=reader.light({**flat_geometry(),"earth_phase_deg":phase})
        assert result["solar"]==pytest.approx(5.0)
        assert result["earth"]==pytest.approx(18.0)
        assert result["total"]==pytest.approx(23.0)


def test_earth_toggle_and_new_phase_preserve_solar_component():
    reader=V.Response(flat_product())
    for actual in (reader.light(flat_geometry(),False),
                   reader.light({**flat_geometry(),"earth_phase_deg":180.0})):
        assert actual["earth"]==0.0
        assert actual["total"]==actual["solar"]==pytest.approx(5.0)


def test_distance_scaling_uses_solid_angle_for_extended_earth():
    reader=V.Response(flat_product())
    reference=reader.light(flat_geometry())
    value=reader.light({**flat_geometry(), "earth_distance_m":50.0,
                       "earth_sunlight_factor":2.0, "sun_distance_factor":0.5})
    q0,q1=(1/100)**2,(1/50)**2
    ratio=(q1/(1+np.sqrt(1-q1)))/(q0/(1+np.sqrt(1-q0)))
    assert value["earth"]==pytest.approx(reference["earth"]*ratio*2)
    assert value["solar"]==pytest.approx(reference["solar"]*.5)


def test_reader_rejects_invalid_export_schema_and_arrays():
    wrong=flat_product();wrong["schema"]="unknown"
    with pytest.raises(ValueError,match="Unsupported"):
        V.Response(wrong)
    wrong=flat_product();wrong["earth"]=copy.deepcopy(wrong["earth"])
    wrong["earth"]["direct_lux"][0][0]=-1
    with pytest.raises(ValueError,match="nonnegative"):
        V.Response(wrong)


def test_saved_reference_reader_has_unique_epochs_and_separate_range_subset():
    records,epochs,native,direction,ranged,mask=V.read_reference()
    assert len(records)==5010
    assert np.all(np.diff(epochs)>0)
    assert int(mask.sum())==2507
    for body in ("sun","earth"):
        np.testing.assert_array_equal(native[f"{body}_distance_m"],direction[f"{body}_distance_m"])
        np.testing.assert_array_equal(ranged[f"{body}_distance_m"][~mask],native[f"{body}_distance_m"][~mask])
        np.testing.assert_allclose(np.linalg.norm(direction[body],axis=1),1.0,atol=1e-14)


def test_historical_regression_uses_only_saved_elevations_and_fixed_source_geometry():
    class Recorder:
        def __init__(self):
            self.calls=[]
        def light(self,g,include_earth=True,order=12):
            self.calls.append((g,include_earth,order))
            # A deliberately different reader still gets a descriptive report;
            # the historical comparison imposes no universal tolerance gate.
            return dict(solar_direct=2.0,solar_diffuse=3.0,solar=5.0)
    reader=Recorder()
    transfer=dict(producer=dict(files={}),transport_settings=dict(increment_tolerance=1e-6,
                  convergence_min_elevation_deg=-90.0,solar_disk_radius_deg=0.0))
    report=V.historical_solar_regression(reader,transfer)
    published=report["published_solved_sky"]
    sea=report["sea_calendar_solar_curve"]
    assert published["samples"]==30
    assert published["supported_elevation_span_deg"]==[-30.0,90.0]
    assert sea["samples"]==361
    assert sea["supported_elevation_span_deg"]==[-90.0,90.0]
    assert len(reader.calls)==391
    assert all(not earth and order==12 for _,earth,order in reader.calls)
    assert all(g["sun_distance_factor"]==1.0 for g,_,_ in reader.calls)
    assert all(g["sun_radius_deg"]==report["normalization"]["solar_disk_radius_deg"] for g,_,_ in reader.calls)
    historical=json.loads(V.HISTORICAL_CURVE.read_text())["sky"]["ground_lux_by_sun_elevation"]
    assert sea["curve"]["elevation_deg"]==historical["elevation_deg"]
    assert sea["curve"]["historical_solar_lux"]==historical["lux"]
    assert report["historical"]["transport_settings"]["increment_tolerance"]==2e-5
    assert report["calendar"]["transport_settings"]["increment_tolerance"]==1e-6


@pytest.mark.skipif(not V.TRANSFER.exists() or not V.OUTPUT.exists(), reason="Full transport export/audit is not yet available")
def test_saved_export_audit_matches_current_inputs_and_python_golden_cases():
    transfer=json.loads(V.TRANSFER.read_text())
    report=json.loads(V.OUTPUT.read_text())
    assert report["schema"]==V.SCHEMA
    for name,expected in report["inputs"].items():
        assert V.digest(V.ROOT/name)==expected,name
    for name,expected in report["producer"]["files"].items():
        assert V.digest(V.ROOT/name)==expected,name
    reader=V.Response(transfer)
    assert len(report["golden_light"])>=40
    for case in report["golden_light"]:
        actual=reader.light(case["geometry"],case["includeEarth"],case["order"])
        for key in V.COMPONENTS:
            assert actual[key]==pytest.approx(case["expected"][key],rel=2e-12,abs=1e-10)
        assert actual["total"]==pytest.approx(actual["solar"]+actual["earth"],rel=1e-14)
