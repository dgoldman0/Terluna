"""Wave-run acceptance, provenance and the numerical AGROW sensitivity.

Small synthetic output fixtures exercise the wrapper's failure handling; they
are not wave-model results or physical validation of SWAN on the Moon.
"""
from dataclasses import asdict, replace
from datetime import timedelta
import hashlib
import json
import math
import os
from pathlib import Path
import re

import numpy as np
import pytest

from climate.waves import build, model
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


@pytest.fixture
def case():
    return model.Case("acceptance", MOON_SURFACE_GRAVITY, 4.0,
                      length_m=2000, cells=2, steps=2, step_s=300,
                      output_every_steps=1, frequency_intervals=4)


def profile_text(case):
    rows = []
    for seconds in range(0, case.duration_s + 1,
                         case.output_every_steps * case.step_s):
        for fetch in np.linspace(0, case.length_m, case.cells + 1):
            height = 0.0 if seconds == 0 else 0.2
            period = -9.0 if seconds == 0 else 5.0
            rows.append(f"{seconds} {fetch} {height} {period} {period} {case.depth_m}")
    return "\n".join(rows) + "\n"


def spectrum_text(case, *, density=(-99.0, 1.0, 2.0, 1.0, -99.0), times=None):
    frequency = np.geomspace(*case.frequencies_hz, case.frequency_intervals + 1)
    if times is None:
        times = range(0, case.duration_s + 1,
                      case.output_every_steps * case.step_s)
    lines = ["SWAN 1", "TIME", "1", "LOCATIONS", "1", f"{case.length_m} 0",
             "AFREQ", str(len(frequency)), *[f"{f:.12g}" for f in frequency],
             "QUANT", "3", "VaDens", "m2/Hz", "-99.0", "CDIR", "degr", "-999",
             "DSPRDEGR", "degr", "-9"]
    for seconds in times:
        date = model.ORIGIN + timedelta(seconds=seconds)
        lines += [date.strftime("%Y%m%d.%H%M%S") + " date and time", "LOCATION 1"]
        values = [-99.0] * len(frequency) if seconds == 0 else density
        lines += [f"{value} 0 20" for value in values]
    return "\n".join(lines) + "\n"


def output_files(case):
    breaking = "\n".join(" ".join(row.split()[:2]) + " 0"
                         for row in profile_text(case).splitlines()) + "\n"
    return {"profile.tbl": profile_text(case), "breaking.tbl": breaking,
            "end.spc": spectrum_text(case),
            "PRINT": "SWAN completed\n", "norm_end": ""}


def fake_executable(path, files):
    """Exit zero after writing supplied files, including misleading norm_end."""
    path.write_text("#!/usr/bin/env python3\nfrom pathlib import Path\n"
                    f"for name, data in {files!r}.items():\n"
                    "    Path(name).write_text(data)\n")
    path.chmod(0o755)
    return path


@pytest.mark.parametrize("failure", ["missing_final", "wrong_time", "negative_height", "negative_period"])
def test_zero_exit_and_norm_end_do_not_accept_invalid_profile(tmp_path, case, failure):
    files = output_files(case)
    rows = [row.split() for row in files["profile.tbl"].splitlines()]
    if failure == "missing_final":
        rows = rows[:-(case.cells + 1)]
    elif failure == "wrong_time":
        for row in rows[-(case.cells + 1):]:
            row[0] = "599"
    elif failure == "negative_height":
        rows[-1][2] = "-0.2"
    else:
        rows[-1][3] = "-9"
    files["profile.tbl"] = "\n".join(" ".join(row) for row in rows) + "\n"
    executable = fake_executable(tmp_path / "fake_swan", files)
    run = tmp_path / "run"
    with pytest.raises(RuntimeError):
        model.run_case(case, executable, run, timeout_s=10)
    assert not (run / "run.json").exists()


def test_zero_initial_spectrum_requires_a_wave_generation_source(case):
    with pytest.raises(ValueError, match="AGROW"):
        replace(case, agrow=False, seed_height_m=0)
    assert "INIT ZERO" in case.input_text()
    assert " AGROW" in case.input_text()
    seeded = replace(case, agrow=False, seed_height_m=0.001)
    assert "INIT PAR" in seeded.input_text()
    assert " AGROW" not in seeded.input_text()


def test_spectral_exception_values_are_bounded_censoring(tmp_path, case):
    path = tmp_path / "end.spc"
    path.write_text(spectrum_text(case))
    result = model.read_spectral_edges(path, case)
    frequency = np.geomspace(*case.frequencies_hz, 5)
    # Only the two edge bins are below SWAN's output floor; negative exception
    # markers must not subtract energy or masquerade as ordinary negative data.
    expected_variance = np.trapezoid([0, 1, 2, 1, 0], frequency)
    expected_bound = math.pi * 1e-12 * (
        frequency[1] - frequency[0] + frequency[-1] - frequency[-2])
    assert result["resolved_variance_m2"] == pytest.approx(expected_variance)
    assert result["below_output_floor_bins"] == 2
    assert result["below_output_floor_variance_bound_m2"] == pytest.approx(expected_bound, rel=1e-9)
    assert 0 < result["low_two_interval_fraction"] < 1
    assert 0 < result["high_two_interval_fraction"] < 1
    assert not result["peak_at_band_edge"]


@pytest.mark.parametrize("times", [(0, 300), (0, 300, 599), (0, 600), (0, 300, 300)])
def test_spectrum_requires_complete_requested_time_coverage(tmp_path, case, times):
    path = tmp_path / "end.spc"
    path.write_text(spectrum_text(case, times=times))
    with pytest.raises(RuntimeError, match="time coverage"):
        model.read_spectral_edges(path, case)


@pytest.mark.parametrize("density", [(-99, 1, -0.1, 1, -99), (-99,) * 5, (0,) * 5])
def test_negative_or_empty_final_spectrum_is_rejected(tmp_path, case, density):
    path = tmp_path / "end.spc"
    path.write_text(spectrum_text(case, density=density))
    with pytest.raises(RuntimeError):
        model.read_spectral_edges(path, case)


@pytest.mark.parametrize("failure", ["missing_row", "nonfinite", "negative", "above_one", "time", "location"])
def test_breaking_fraction_requires_physical_values_and_matching_samples(tmp_path, case, failure):
    files = output_files(case)
    profile_path = tmp_path / "profile.tbl"
    profile_path.write_text(files["profile.tbl"])
    profile = model.read_profile(profile_path, case)
    path = tmp_path / "breaking.tbl"
    path.write_text(files["breaking.tbl"])
    assert model.breaking_fraction(path, profile) == 0
    rows = [row.split() for row in files["breaking.tbl"].splitlines()]
    if failure == "missing_row":
        rows.pop()
    elif failure == "nonfinite":
        rows[-1][2] = "nan"
    elif failure == "negative":
        rows[-1][2] = "-0.01"
    elif failure == "above_one":
        rows[-1][2] = "1.01"
    elif failure == "time":
        rows[-1][0] = "599"
    else:
        rows[-1][1] = "1900"
    path.write_text("\n".join(" ".join(row) for row in rows) + "\n")
    with pytest.raises(RuntimeError):
        model.breaking_fraction(path, profile)


def cached_run(tmp_path, case):
    """Create a cache fixture without claiming that a model was executed."""
    executable = tmp_path / "unused_executable"
    executable.write_text("This cache fixture must never execute.\n")
    directory = tmp_path / "cached"
    directory.mkdir()
    files = output_files(case)
    files["INPUT"] = case.input_text()
    files["depth.bot"] = f"{case.depth_m:.12g} {case.depth_m:.12g}\n"
    for name, text in files.items():
        (directory / name).write_text(text)
    record = {"schema": "terluna.climate.swan-run/1",
              "identity": {"case": asdict(case), "executable_sha256": model.sha256(executable),
                           "input_sha256": hashlib.sha256(case.input_text().encode()).hexdigest(),
                           "depth_sha256": hashlib.sha256(files["depth.bot"].encode()).hexdigest()},
              "output_sha256": {name: model.sha256(directory / name)
                                for name in output_files(case)},
              "spectral_edges": {"stale_analysis": True}}
    (directory / "run.json").write_text(json.dumps(record))
    return executable, directory


@pytest.mark.parametrize("changed", ["executable", "case", "profile", "INPUT", "depth.bot"])
def test_cache_rejects_stale_inputs_or_changed_bytes(tmp_path, case, changed):
    executable, directory = cached_run(tmp_path, case)
    requested = case
    if changed == "executable":
        executable.write_text("Different build.\n")
    elif changed == "case":
        requested = replace(case, wind_m_s=case.wind_m_s + 1)
    elif changed == "profile":
        with (directory / "profile.tbl").open("a") as stream:
            stream.write("0 0 0 0 0 0\n")
    elif changed == "INPUT":
        (directory / changed).write_text(case.input_text().replace("WIND 4 0", "WIND 40 0"))
    else:
        (directory / changed).write_text("100 100\n")
    with pytest.raises((RuntimeError, FileExistsError)):
        model.run_case(requested, executable, directory)


def test_cached_spectrum_is_reparsed_instead_of_trusting_old_analysis(tmp_path, case):
    executable, directory = cached_run(tmp_path, case)
    record = model.run_case(case, executable, directory)
    assert "stale_analysis" not in record["spectral_edges"]
    assert record["spectral_edges"]["resolved_variance_m2"] > 0
    # Rehashing truncated output does not make incomplete time coverage valid.
    (directory / "end.spc").write_text(spectrum_text(case, times=(0, 300)))
    record["output_sha256"]["end.spc"] = model.sha256(directory / "end.spc")
    (directory / "run.json").write_text(json.dumps(record))
    with pytest.raises(RuntimeError, match="time coverage"):
        model.run_case(case, executable, directory)


def test_source_hash_mismatch_leaves_original_file_untouched(tmp_path):
    source = tmp_path / "swancom3.ftn"
    original = b"            FREQ1 = SIGMA/PI2                                             40.88\n"
    source.write_bytes(original)
    with pytest.raises(ValueError, match="pinned original"):
        build.patch_agrow(tmp_path)
    assert source.read_bytes() == original


def test_agrow_frequency_rescaling_preserves_earth_and_fixed_f_over_g(tmp_path, monkeypatch):
    # An isolated source fragment exercises patch arithmetic; its local hash
    # override does not admit this fragment as the pinned upstream source.
    source = tmp_path / "swancom3.ftn"
    source.write_text("            FREQ1 = SIGMA/PI2                                             40.88\n")
    before = build.sha256(source)
    monkeypatch.setattr(build, "AGROW_SOURCE_SHA256", before)
    record = build.patch_agrow(tmp_path)
    expression = re.search(r"FREQ1\s*=\s*(.+)", source.read_text()).group(1)

    def knee_argument(frequency, gravity):
        return eval(expression, {"__builtins__": {}},
                    {"SIGMA": 2 * math.pi * frequency, "PI2": 2 * math.pi, "GRAV": gravity})

    for frequency in (0.01, 0.1, 1.0):
        assert knee_argument(frequency, STANDARD_GRAVITY) == pytest.approx(frequency)
        lunar_frequency = frequency * MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
        assert knee_argument(lunar_frequency, MOON_SURFACE_GRAVITY) == pytest.approx(frequency)
    assert record["before_sha256"] == before
    assert record["after_sha256"] == build.sha256(source)
    assert record["after_sha256"] != before


@pytest.mark.skipif(not os.environ.get("TERLUNA_SWAN_EXECUTABLE"),
                    reason="Set TERLUNA_SWAN_EXECUTABLE for a real bounded SWAN run")
def test_real_swan_writes_complete_nonnegative_wave_growth(tmp_path):
    executable = Path(os.environ["TERLUNA_SWAN_EXECUTABLE"])
    case = model.Case("integration", STANDARD_GRAVITY, 8.0,
                      length_m=10_000, cells=10, steps=12, step_s=300,
                      output_every_steps=6, frequency_intervals=24)
    record = model.run_case(case, executable, tmp_path / "integration", timeout_s=60)
    profile = model.read_profile(tmp_path / "integration/profile.tbl", case)
    assert profile[-1, -1, 2] > 0
    assert record["spectral_edges"]["resolved_variance_m2"] > 0
    assert record["maximum_breaking_fraction"] == 0
    assert record["threads"] == 1
