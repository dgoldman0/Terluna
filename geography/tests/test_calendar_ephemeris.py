"""Independent pinned JPL checks for the lighting calendar's date range."""
import json

import numpy as np
import pytest

from geography import calendar_ephemeris_check as C


@pytest.fixture(scope="module")
def report():
    return C.compare_fixture()


def test_reference_spans_calendar_and_preserves_direction_accuracy(report):
    assert report["samples"] == 5021
    assert report["span_calendar_tt"] == ["2000-01-01T00:00:00", "2500-12-31T23:30:00"]
    for body in ("sun", "earth"):
        assert report["errors"][f"{body}_direction_deg"]["max_abs"] < 0.07
        assert report["errors"][f"{body}_inverse_square_flux_relative"]["max_abs"] < 0.0002
    assert report["errors"]["earth_distance_km"]["max_abs"] < 20
    assert report["errors"]["sun_distance_km"]["max_abs"] < 15000


def test_horizon_and_twilight_timing_distinguish_grazing_events(report):
    assert len(report["horizon_events"]) == 14
    assert report["horizon_event_abs_error_max_minutes"] < 10
    assert report["event_interpolation_sensitivity_max_seconds"] < 0.2
    light = report["photometry"]
    assert len(light["threshold_events"]) == 73
    assert light["event_count_differences"] == []
    grazing = [e for e in light["threshold_events"] if e["site"] == "W Procellarum" and e["threshold_lux"] == 0.1]
    assert len(grazing) == 2
    assert all(e["reference_step_minutes"] == 1 for e in grazing)
    # A small angular residual can move a near-tangent threshold by an hour.
    assert all(45 < abs(e["model_minus_reference_minutes"]) < 75 for e in grazing)
    assert all(e["reference_interpolation_change_seconds"] < 0.5 for e in grazing)


def test_direction_residual_propagates_to_small_photometric_error(report):
    light = report["photometry"]
    assert len(light["sites"]) == 11
    for site in light["sites"]:
        for band in site["bands"]:
            assert band["max_relative_error_percent"] < 1.6


def test_fixture_tampering_is_rejected(tmp_path):
    copy = tmp_path / "reference.csv"
    copy.write_bytes(C.CSV_PATH.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        C.compare_fixture(copy, C.INPUT_PATH)


def test_tt_calendar_arithmetic_and_longitude_wrap():
    assert C.calendar_jd("2000-01-01") == 2451544.5
    assert C.calendar_jd("2000-03-01") - C.calendar_jd("2000-02-28") == 2
    assert C.calendar_jd("2100-03-01") - C.calendar_jd("2100-02-28") == 1
    a = C.ephemeris.unit(np.radians(179.99), 0.0)
    b = C.ephemeris.unit(np.radians(-179.99), 0.0)
    assert abs(C.angular_separation_deg(a, b) - 0.02) < 1e-10


def test_fetch_provenance_records_tt_and_mean_earth_orientation():
    metadata = json.loads(C.INPUT_PATH.read_text())
    assert len(metadata["responses"]) == 21
    for response in metadata["responses"]:
        assert response["parameters"]["TIME_TYPE"] == "TT"
        assert len(response["response_sha256"]) == 64
        if response["parameters"]["EPHEM_TYPE"] == "OBSERVER":
            frames = [s for s in response["headers"] if s.startswith("Target pole/equ")]
            assert len(frames) == 1 and "MOON_ME" in frames[0] and "IAU_MOON" not in frames[0]
