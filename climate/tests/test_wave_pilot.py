"""Independent thermodynamics/propagation checks and refusal of corrupt inputs."""
from pathlib import Path
import hashlib
import json
import os

import numpy as np
import pytest

from climate.crm.wave_forcing import FIELDS, moist_density, read_surface
from climate.waves import pilot
from climate.waves.pilot_checks import spectrum_diagnostics
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


def test_moist_density_dry_limit_and_vapour_displacement():
    dry = moist_density(120000., 300., 0.)
    assert dry == pytest.approx(120000 / (CM1_DRY_AIR_GAS_CONSTANT * 300))
    assert 0 < moist_density(120000., 300., .02) < dry
    assert moist_density(240000., 300., .02) == pytest.approx(2 * moist_density(120000., 300., .02))
    with pytest.raises(ValueError):
        moist_density(120000., 300., -.001)


def test_surface_reader_hashes_only_identified_fields_and_rejects_truncation(tmp_path):
    layout = [("unused", 3), *[(key, 1) for key in reversed(FIELDS)]]
    values = np.arange(sum(n for _, n in layout) * 4, dtype="<f4")
    path = tmp_path / "snapshot.dat"
    values.tofile(path)
    fields, record = read_surface(path, 4, layout)
    assert np.array_equal(fields["q2"], values[12:16])
    assert record["selected_bytes_sha256"] == hashlib.sha256(b"".join(fields[k].tobytes() for k in FIELDS)).hexdigest()
    values[:-1].tofile(path)
    with pytest.raises(ValueError):
        read_surface(path, 4, layout)


def test_basin_rejects_a_foreign_schema_and_nonpositive_depth(tmp_path):
    path = tmp_path / "atlas.npz"
    a = np.zeros((8, 8), int)
    a[2:6, 2:6] = 874
    fields = dict(lat_deg=np.arange(8)[::-1], lon_deg=np.arange(8), water_label=a, height_m=np.full((8, 8), -5.))
    np.savez(path, metadata=json.dumps(dict(schema="foreign", sea_level_m=0)), **fields)
    with pytest.raises(ValueError, match="schema|atlas"):
        pilot.load_basin(path, stride=1)
    np.savez(path, metadata=json.dumps(dict(schema="terluna.geography.atlas-grid/1", sea_level_m=-10)), **fields)
    with pytest.raises(ValueError, match="depth"):
        pilot.load_basin(path, stride=1)


def test_basin_table_rejects_missing_times_wrong_bottom_and_missing_wet_points(tmp_path):
    case = pilot.BasinCase("synthetic", hours=1)
    basin = dict(depth=np.array([[10., 20.], [30., -99999.]]), wet=np.array([[True, True], [True, False]]),
                 points=np.array([[0., 0.], [1., 0.], [0., 1.]]))
    rows = []
    for t in (0, 3600):
        for j in range(2):
            for i in range(2):
                rows.append([t, i, j, 1., 10., 8., basin["depth"][j, i], 0., 0.])
    values = np.asarray(rows)
    path = tmp_path / "basin.tbl"
    np.savetxt(path, values)
    assert pilot.read_basin(path, case, basin).shape == (2, 3, 9)
    for bad in (values[:-1], values.copy(), values.copy()):
        if bad.shape == values.shape:
            bad[-4, 6] = -99 if bad[-4, 6] == 10 else 15
        np.savetxt(path, bad)
        with pytest.raises(ValueError):
            pilot.read_basin(path, case, basin)


def test_offshore_spectrum_requires_final_time_and_valid_energy(tmp_path):
    path = tmp_path / "offshore.spc"
    text = "\n".join(["SWAN 1", "TIME", "1", "LONLAT", "1", "89.125 -1.125",
                       "AFREQ", "5", ".01", ".02", ".04", ".08", ".16",
                       "QUANT", "3", "VaDens", "m2/Hz", "-99", "CDIR", "degr", "-999",
                       "DSPRDEGR", "degr", "-9", "20000101.120000 date and time", "LOCATION 1",
                       "-99 -999 -9", "1 0 20", "2 0 20", ".1 0 20", "-99 -999 -9"]) + "\n"
    path.write_text(text)
    assert spectrum_diagnostics(path, "20000101.120000")["resolved_hs_m"] > 0
    with pytest.raises(ValueError, match="final time"):
        spectrum_diagnostics(path, "20000101.130000")
    path.write_text(text.replace("2 0 20", "-2 0 20"))
    with pytest.raises(ValueError, match="variance"):
        spectrum_diagnostics(path, "20000101.120000")


@pytest.fixture
def executable():
    value = os.environ.get("TERLUNA_SWAN_EXECUTABLE")
    if not value:
        pytest.skip("Set TERLUNA_SWAN_EXECUTABLE for numerical SWAN controls")
    return Path(value).resolve(strict=True)


def test_lunar_spherical_radius_matches_cartesian_length(executable, tmp_path):
    """Same small physical basin expressed in metres and lunar degrees."""
    outputs = []
    for coordinate in ("cartesian", "spherical"):
        size = 2000. if coordinate == "cartesian" else np.degrees(2000. / MOON_RADIUS)
        coordinates = "CARTESIAN" if coordinate == "cartesian" else f"SPHERICAL {MOON_RADIUS} CCM"
        ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
        text = f"""PROJECT 'Radius control' '004'
MODE NONSTAT TWODIMENSIONAL
COORD {coordinates}
SET GRAV {MOON_SURFACE_GRAVITY:.12g}
OUTPUT OPTIONS TABLE 16
CGRID REG 0 0 0 {size:.12g} {size:.12g} 4 4 CIRCLE 36 {0.03*ratio} {3*ratio} 48
INPGRID BOTTOM 0 0 0 1 1 {size:.12g} {size:.12g}
READINP BOTTOM 1 'depth.bot' 3 0 FREE
WIND 4 0
GEN3 KOMEN DRAG WU AGROW
BREAKING CONSTANT 1 .73
PROP BSBT
INIT ZERO
TABLE 'COMPGRID' NOHEAD 'control.tbl' HS TM01 DEP OUTPUT 20000101.010000 1 HR
COMPUTE NONSTAT 20000101.000000 30 SEC 20000101.010000
STOP
"""
        directory = tmp_path / coordinate
        pilot.execute(executable, directory, {"INPUT": text, "depth.bot": "1000 1000\n1000 1000\n"}, ["control.tbl"])
        outputs.append(np.loadtxt(directory / "control.tbl"))
    assert outputs[0].shape == outputs[1].shape == (25, 3)
    mask = outputs[0][:, 0] > .01
    assert mask.any()
    np.testing.assert_allclose(outputs[0][mask, :2], outputs[1][mask, :2], rtol=.01)


def test_slope_conserves_linear_wave_energy_flux_before_breaking(executable, tmp_path):
    """Small single-bin waves: Hs² Cg should be constant without dissipation."""
    files = pilot.slope_files("", .01, cells=800)
    files.pop("boundary.spc")
    files["INPUT"] = files["INPUT"].replace("BOUNDSPEC SIDE WEST CONSTANT FILE 'boundary.spc' 1",
        "BOUND SHAPE BIN PEAK DSPR POWER\nBOUNDSPEC SIDE WEST CONSTANT PAR .02 10 0 128")
    files["INPUT"] = files["INPUT"].replace("CIRCLE 36", "CIRCLE 72")
    pilot.execute(executable, tmp_path / "shoaling", files, ["slope.tbl"])
    values = pilot.read_slope(tmp_path / "shoaling/slope.tbl", files, 800)
    mask = values[:, 4] >= 5
    assert np.max(values[mask, 6]) == 0
    depth = values[mask, 4]
    omega = 2 * np.pi / values[mask, 2]
    low, high = np.zeros_like(depth), np.full_like(depth, 100.)
    for _ in range(80):
        k = .5 * (low + high)
        residual = MOON_SURFACE_GRAVITY * k * np.tanh(k * depth) - omega**2
        high = np.where(residual > 0, k, high)
        low = np.where(residual <= 0, k, low)
    kd = k * depth
    cg = .5 * omega / k * (1 + 2 * kd / np.sinh(2 * kd))
    flux = values[mask, 1]**2 * cg
    np.testing.assert_allclose(flux / flux[0], 1, rtol=.02)
