"""Check geographic wind rotation and matched coastal comparisons."""
import os
from pathlib import Path

import numpy as np
import pytest

from climate.crm.cm1_run import tilted_path
from climate.crm.wave_coverage import geographic_vectors
from climate.gcm.wave_coverage import mean_output_settings
from climate.waves.coastal import compare
from climate.waves.coastal_checks import nesting_spectra
from climate.waves.pilot import execute
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


def test_ring_rotation_matches_cartesian_tangent_geometry():
    count = 720
    track = dict(tilted_path(count, 70, 45), tilt_deg=70, node_deg=45)
    latitude, longitude = np.radians(track["lat"]), np.radians(track["lon"])
    r = np.column_stack((np.cos(latitude)*np.cos(longitude), np.cos(latitude)*np.sin(longitude), np.sin(latitude)))
    tangent = np.roll(r,-1,axis=0)-np.roll(r,1,axis=0)
    tangent /= np.linalg.norm(tangent,axis=1)[:,None]
    east_basis = np.column_stack((-np.sin(longitude),np.cos(longitude),np.zeros(count)))
    north_basis = np.cross(r,east_basis)
    east,north = geographic_vectors(np.ones(count),np.zeros(count),track)
    np.testing.assert_allclose(east,np.sum(tangent*east_basis,axis=1),atol=1e-12)
    np.testing.assert_allclose(north,np.sum(tangent*north_basis,axis=1),atol=1e-12)
    left = np.cross(r,tangent)
    east,north = geographic_vectors(np.zeros(count),np.ones(count),track)
    np.testing.assert_allclose(east,np.sum(left*east_basis,axis=1),atol=1e-12)
    np.testing.assert_allclose(north,np.sum(left*north_basis,axis=1),atol=1e-12)


def test_equatorial_rotation_and_corrupt_positions():
    track = dict(tilted_path(32,0,0),tilt_deg=0,node_deg=0)
    east,north = geographic_vectors(np.full((2,32),3.),np.full((2,32),4.),track)
    np.testing.assert_array_equal(east,3.)
    np.testing.assert_array_equal(north,4.)
    track["lon"][2]+=1
    with pytest.raises(ValueError,match="positions"):
        geographic_vectors(np.ones(32),np.ones(32),track)


def test_gcm_mean_export_refuses_instantaneous_output_settings():
    text = " NLOWIO = 1\n NSNAPSHOT = 0\n NHCADENCE = 0\n"
    assert mean_output_settings(text)["NLOWIO"] == 1
    with pytest.raises(ValueError, match="NLOWIO=1"):
        mean_output_settings(text.replace("NLOWIO = 1", "NLOWIO = 0"))
    with pytest.raises(ValueError, match="Missing"):
        mean_output_settings("NSNAPSHOT = 0\nNHCADENCE = 0\n")


def test_coastal_error_excludes_boundary_and_pairs_only_common_wet_nodes():
    coarse=np.ones((13,9,9,9));fine=np.ones((13,17,17,9))
    coarse[...,3]=2.;fine[...,3]=2.
    coarse[...,5]=8.;fine[...,5]=8.
    coarse[:,:2,:,3]=20.
    wet=np.ones((9,9),bool);refined_wet=np.ones((17,17),bool)
    wet[4,4]=False
    coarse[:,4,4,3]=200.
    result=compare(coarse,fine,wet,refined_wet)
    assert result["matched_interior_nodes"]==24
    assert result["hs_relative_max"]==0
    assert result["within_five_percent"]
    fine[:,6,6,5]=10.
    result=compare(coarse,fine,wet,refined_wet)
    assert result["tm01_relative_max"]==pytest.approx(.2)
    assert not result["within_five_percent"]


def test_nested_spectra_reject_truncated_and_negative_energy(tmp_path):
    path = tmp_path/"boundary.nest"
    header = "\n".join(("SWAN 1", "TIME", "1", "LONLAT", "2", "0 0", "1 0", "RFREQ", "3",
                        ".01", ".02", ".04", "CDIR", "4", "0", "90", "180", "270",
                        "QUANT", "1", "VaDens", "m2/Hz/degr", "-99", "20000101.000000"))
    text = header+"\nFACTOR\n0.01\n1 2 3 4\n4 3 2 1\n0 0 0 0\nNODATA\n"
    path.write_text(text)
    data = nesting_spectra(path)
    assert data["available"].tolist() == [[True, False]]
    assert data["spectra"][0, 0, 0, 3] == pytest.approx(.04)
    path.write_text(text.replace("1 2 3 4", "-1 2 3 4"))
    with pytest.raises(ValueError, match="variance"):
        nesting_spectra(path)
    path.write_text(text.removesuffix("NODATA\n"))
    with pytest.raises(ValueError, match="Truncated"):
        nesting_spectra(path)


def test_full_directional_nesting_preserves_uniform_depth_wind_sea(tmp_path):
    value = os.environ.get("TERLUNA_SWAN_EXECUTABLE")
    if not value:
        pytest.skip("Set TERLUNA_SWAN_EXECUTABLE for the nesting control")
    executable = Path(value).resolve(strict=True)
    size = np.degrees(4000/MOON_RADIUS)
    ratio = MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    results = []
    for name, origin, length, cells in (("parent", 0., size, 8), ("child", size/4, size/2, 4)):
        boundary = "" if name == "parent" else "BOUNDNEST1 NEST 'boundary.nest' CLOSED\n"
        output = (f"NGRID 'CHILD' {size/4} {size/4} 0 {size/2} {size/2} 4 4\n"
                  "NESTOUT 'CHILD' 'boundary.nest' OUTPUT 20000101.000000 5 MIN\n") if name == "parent" else ""
        text = f"""PROJECT 'Nesting control' '006'
MODE NONSTAT TWODIMENSIONAL
COORD SPHERICAL {MOON_RADIUS} CCM
SET GRAV {MOON_SURFACE_GRAVITY:.12g}
OUTPUT OPTIONS TABLE 16
CGRID REG {origin} {origin} 0 {length} {length} {cells} {cells} CIRCLE 36 {0.03*ratio} {3*ratio} 48
INPGRID BOTTOM 0 0 0 1 1 {size} {size}
READINP BOTTOM 1 'depth.bot' 3 0 FREE
{boundary}WIND 4 0
GEN3 KOMEN DRAG WU AGROW
BREAKING CONSTANT 1 .73
PROP BSBT
INIT ZERO
TABLE 'COMPGRID' NOHEAD 'control.tbl' HS TM01 DEP OUTPUT 20000101.020000 1 HR
{output}COMPUTE NONSTAT 20000101.000000 30 SEC 20000101.020000
STOP
"""
        files = {"INPUT": text, "depth.bot": "1000 1000\n1000 1000\n"}
        outputs = ["control.tbl"]
        if name == "child":
            files["boundary.nest"] = (tmp_path/"parent/boundary.nest").read_text()
        else:
            outputs.append("boundary.nest")
        execute(executable, tmp_path/name, files, outputs)
        results.append(np.loadtxt(tmp_path/name/"control.tbl").reshape(cells+1, cells+1, 3))
    parsed = nesting_spectra(tmp_path/"parent/boundary.nest")
    assert parsed["spectra"].shape == (25, 16, 49, 36)
    assert parsed["available"].all()
    assert np.max(results[0][..., 0]) > .01
    np.testing.assert_allclose(results[1][1:-1, 1:-1, :2], results[0][3:6, 3:6, :2], rtol=.02)
