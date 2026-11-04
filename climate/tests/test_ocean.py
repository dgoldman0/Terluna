"""The ocean models' set-up: lunar constants, grids, meshes, the tide's forcing files and the SCHISM patch.

None of these tests runs a model; those that need the pinned sources on the research drive skip without them.
"""
import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

from climate.ocean import builds, checks, forcing, lunar, mesh, mitgcm
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, SIDEREAL_MONTH_DAYS

ROOT = Path(__file__).resolve().parents[2]


def test_lunar_constants_come_from_shared():
    assert lunar.GRAVITY == MOON_SURFACE_GRAVITY
    assert lunar.RADIUS == MOON_RADIUS
    assert lunar.ROTATION_PERIOD_S == pytest.approx(SIDEREAL_MONTH_DAYS * 86400)
    assert lunar.OMEGA == pytest.approx(2.6617e-6, rel=1e-4)
    assert 995.0 < lunar.WATER_DENSITY < 997.0          # fresh water near 30 C


def test_schism_patch_compiles_the_shared_constants():
    patch = (ROOT / "climate/ocean/schism/lunar.patch").read_text()
    grav = float(re.search(r"^\+\s*real\(rkind\),parameter :: grav=([0-9.eEdD+-]+)_rkind", patch, re.M).group(1))
    omega = float(re.search(r"^\+.*omega_e=real\(([0-9.eEdD+-]+),rkind\)", patch, re.M).group(1).replace("d", "e"))
    assert grav == pytest.approx(MOON_SURFACE_GRAVITY, rel=1e-15)
    assert omega == pytest.approx(2 * math.pi / (SIDEREAL_MONTH_DAYS * 86400), rel=1e-15)
    assert "rearth_pole=1737400._rkind" in patch and str(int(MOON_RADIUS)) == "1737400"
    # Earth's constituents and their 0.69 Love factor are refused, and sflux's Earth bulk formulas too
    assert "lunar build, set ntip=0" in patch and "sflux (nws=2) carries Earth constants" in patch


def test_sources_pin_commits_and_hashes():
    sources = json.loads((ROOT / "climate/ocean/sources.json").read_text())
    for model in ("mitgcm", "schism"):
        assert re.fullmatch(r"[0-9a-f]{40}", sources["models"][model]["commit"])
    assert all(re.fullmatch(r"[0-9a-f]{64}", p["sha256"]) for p in sources["libraries"]["packages"])


@pytest.mark.skipif(not (builds.DRIVE / "src" / "schism" / ".git").exists(), reason="SCHISM checkout not on this machine")
def test_lunar_patch_applies_to_the_pinned_source(tmp_path):
    import shutil
    import subprocess
    source = builds.checkout("schism")
    for name in ("Core/schism_glbl.F90", "Core/scribe_io.F90", "Hydro/schism_init.F90", "Hydro/schism_step.F90"):
        (tmp_path / "src" / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(source / "src" / name, tmp_path / "src" / name)
    with open(builds.PATCH) as patch:
        done = subprocess.run(["patch", "-p1", "--dry-run", "-d", str(tmp_path)], stdin=patch, capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr


def test_stretched_edges():
    e = mitgcm.stretched_edges(-10.0, 20.0, (0.0, 5.0), 1 / 16, 1 / 4)
    d = np.diff(e)
    assert e[0] == -10.0 and e[-1] == 20.0 and np.all(d > 0)
    inside = (e[:-1] >= 0.0) & (e[1:] <= 5.0)
    assert np.allclose(d[inside], 1 / 16)
    assert d[1:-1].max() <= 0.25 + 1e-9 and d.max() <= 1.5 * 0.25   # the end cells absorb the remainder
    growth = d[1:] / d[:-1]
    assert growth[(d[1:] < 0.25 - 1e-9) & (d[:-1] < 0.25 - 1e-9)].max() <= 1.1 + 1e-9


def test_aggregate_keeps_cells_with_most_water():
    lat = np.arange(0.0, 2.0, 0.25) + 0.125
    lon = np.arange(0.0, 2.0, 0.25) + 0.125
    height = np.full((8, 8), -50.0)
    height[:, :3] = 10.0                                # three dry columns
    box = dict(lat=lat, lon=lon, height=height, level=0.0)
    depth = mitgcm.aggregate(box, np.ones((8, 8), bool), np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
    assert depth.shape == (2, 2)
    assert np.all(depth[:, 0] == 0.0)                   # 3 of 4 columns dry
    assert np.allclose(depth[:, 1], 50.0)


def test_mitgcm_namelists_carry_the_moon(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    grid = mitgcm.Grid(xg=np.arange(0.0, 2.01, 0.25), yg=np.arange(-1.0, 1.01, 0.25), depth=np.full((8, 8), 100.0))
    exp = mitgcm.Experiment(name="t", grid=grid, deltat=60.0, tide=dict(lon0=-3.0, lat0=-4.0, step=1.0, nlon=9, nlat=9,
                                                                      period_s=3600.0, file="tide.bin"))
    exp.write()
    data = (exp.folder / "data.template").read_text()
    assert f"gravity={lunar.GRAVITY:.12g}" in data and f"rotationPeriod={lunar.ROTATION_PERIOD_S:.12g}" in data
    assert "rSphere=1737400" in data and "momTidalForcing=.TRUE." in data
    exf = (exp.folder / "data.exf").read_text()
    assert f"exf_inscal_tidePot={-lunar.GRAVITY:.12g}" in exf and "tidePotstartdate1=20260101" in exf
    assert "selectImplicitDrag=2" in data
    bathy = np.fromfile(exp.folder / "bathy.bin", dtype=">f8").reshape(8, 8)
    assert np.all(bathy == -100.0)
    # a resumed run must find the setup it started with
    from dataclasses import replace
    assert exp.setup_digest("build") == exp.setup_digest("build")
    assert exp.setup_digest("build") != replace(exp, deltat=30.0).setup_digest("build")
    assert exp.setup_digest("build") != exp.setup_digest("other build")


def test_mitgcm_compiles_implicit_drag_and_zstar_on_request(tmp_path, monkeypatch):
    inc = tmp_path / "model" / "inc"
    inc.mkdir(parents=True)
    (inc / "CPP_OPTIONS.h").write_text("#undef ALLOW_SOLVE4_PS_AND_DRAG\n#define INCLUDE_IMPLVERTADV_CODE\n"
                                       "#undef NONLIN_FRSURF\nc#define DISABLE_RSTAR_CODE\n")
    linear = mitgcm.cpp_options(tmp_path)
    assert "#define ALLOW_SOLVE4_PS_AND_DRAG" in linear and "#undef NONLIN_FRSURF" in linear
    assert "#define NONLIN_FRSURF" in mitgcm.cpp_options(tmp_path, nonlinear_surface=True)
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    grid = mitgcm.Grid(xg=np.arange(0.0, 2.01, 0.25), yg=np.arange(-1.0, 1.01, 0.25), depth=np.full((8, 8), 100.0))
    texts, _ = mitgcm.Experiment(name="z", grid=grid, deltat=60.0, nonlinear=True, nonlinear_surface=True).inputs()
    assert "nonlinFreeSurf=4" in texts["data.template"] and "select_rStar=2" in texts["data.template"]
    texts, _ = mitgcm.Experiment(name="l", grid=grid, deltat=60.0).inputs()
    assert "nonlinFreeSurf=0" in texts["data.template"] and "select_rStar" not in texts["data.template"]


def test_mds_reader(tmp_path):
    values = np.arange(12, dtype=">f4").reshape(3, 4)
    np.concatenate([values, 2 * values]).astype(">f4").tofile(tmp_path / "state.0000000010.data")
    (tmp_path / "state.0000000010.meta").write_text(
        " nDims = [   2 ];\n dimList = [\n     4,    1,    4,\n     3,    1,    3\n ];\n dataprec = [ 'float32' ];\n"
        " nrecords = [     2 ];\n timeStepNumber = [         10 ];\n nFlds = [    2 ];\n fldList = {\n 'ETAN    ' 'UVEL    '\n };\n")
    out = mitgcm.read_mds(tmp_path / "state.0000000010")
    assert np.array_equal(out["ETAN"], values) and np.array_equal(out["UVEL"], 2 * values)


def test_cells_mesh_splits_corner_contacts_and_keeps_area():
    depth = np.array([[10.0, 0.0], [0.0, 10.0]])        # two cells touching at one corner
    m = mesh.from_cells(depth, np.array([0.125, 0.375]), np.array([0.125, 0.375]))
    assert m.elements == 4 and m.nodes == 8             # the shared corner is two nodes
    assert len(m.boundary_loops()) == 2
    cells = np.ones((3, 3)) * 20.0
    m = mesh.from_cells(cells, np.arange(3) * 0.25, np.arange(3) * 0.25)
    assert m.elements == 18 and len(m.boundary_loops()) == 1
    expected = mitgcm.Grid(xg=np.arange(4) * 0.25 - 0.125, yg=np.arange(4) * 0.25 - 0.125, depth=cells).area.sum()
    assert m.element_areas(lunar.RADIUS).sum() == pytest.approx(expected, rel=1e-3)


def test_quadtree_mesh_has_no_hanging_nodes():
    ny = nx = 32                                         # 4 x 4 base cells of 8 pixels
    depth = np.full((ny, nx), 100.0)
    sea = np.ones((ny, nx), bool)
    sea[10:14, 20:25] = False                           # an island
    level = np.zeros((4, 4), int)
    level[1, 1] = 3                                     # one fine base cell; balance raises its neighbours
    lat = (np.arange(ny) + 0.5) / 16
    lon = (np.arange(nx) + 0.5) / 16
    m = mesh.build(depth, sea, level, lat, lon)
    edges = np.sort(np.concatenate([m.triangles[:, [0, 1]], m.triangles[:, [1, 2]], m.triangles[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert counts.max() == 2                             # every edge shared by at most two triangles
    loops = m.boundary_loops()
    assert len(loops) == 2                               # the outer coast and the island
    water = sea.sum() * (lunar.RADIUS * math.radians(1 / 16)) ** 2
    assert m.element_areas(lunar.RADIUS).sum() == pytest.approx(water, rel=0.02)


def test_forcing_matches_geography_and_writes_schism_files(tmp_path):
    from geography import tides
    lat, lon = np.array([-6.0, 36.0]), np.array([54.0, -16.0])
    seconds = np.array([0.0, 3600.0, 7200.0])
    c = forcing.coefficients(seconds)
    assert np.allclose(c, tides.forcing_coefficients(tides.JD_START + seconds / 86400.0))
    assert np.allclose(forcing.monomials(lat, lon) @ c[0],
                       tides.monomial_values(tides.unit_points(np.radians(lat), np.radians(lon))) @ c[0])
    forcing.write_schism(tmp_path, lat, lon, seconds, c, attraction=False)
    head = (tmp_path / "lunar_tide.th").read_text().splitlines()[0].split()
    assert head == ["3", "20", "3600.000000"]
    basis = np.loadtxt(tmp_path / "lunar_tide_basis.in", skiprows=1)
    assert basis.shape == (2, 20) and np.allclose(basis @ c[1], forcing.monomials(lat, lon) @ c[1], atol=1e-8)
    assert forcing.ramp(np.array([0.0]), 5.0)[0] == 0.0 and forcing.ramp(np.array([1e9]), 5.0)[0] == pytest.approx(1.0)


def test_line_fit_and_frequency():
    t = np.arange(0.0, 20 * 86400.0, 900.0)
    x = 0.3 * np.cos(2 * math.pi * 0.6 * t / 86400.0 + 0.4) + 0.05
    mean, amp = checks.fit_lines(t, x, [0.6])
    assert mean == pytest.approx(0.05, abs=1e-9) and abs(amp[0]) == pytest.approx(0.3, rel=1e-9)
    assert np.angle(amp[0]) == pytest.approx(0.4, abs=1e-9)
    assert checks.dominant_frequency(t, x, 0.55) == pytest.approx(0.6, rel=1e-6)


def test_reducer_counts_each_snapshot_once(tmp_path):
    from climate.ocean import fecunditatis
    lat, lon = np.array([0.0, 7.0, 40.0]), np.array([50.0, 58.0, 40.0])
    r = fecunditatis.Reducer(tmp_path / "r.npz", lat, lon, np.ones(3), np.array([1, 1, 2]), 1)
    t0 = fecunditatis.SPINUP_DAYS * 86400.0
    r.add(t0, np.array([1.0, 3.0, 5.0]), np.array([0.1, 2.0, 0.2]))
    r.add(t0, np.array([9.0, 9.0, 9.0]), np.array([9.0, 9.0, 9.0]))     # a repeat is ignored
    r.add(t0 + 10800.0, np.array([0.0, 1.0, 2.0]), np.array([0.3, 0.5, 0.1]))
    r.save()
    s = fecunditatis.Reducer(tmp_path / "r.npz", lat, lon, np.ones(3), np.array([1, 1, 2]), 1).state
    assert np.allclose(s["fec_mean"], [2.0, 0.5]) and s["count"][0] == 2
    assert np.allclose(s["month_max"][0], [1.0, 3.0, 5.0]) and np.allclose(s["month_min"][0], [0.0, 1.0, 2.0])
    assert np.allclose(s["zone_speed_1"], [2.0, 0.5])                   # the point at 7 N 58 E is in the Crisium passage


RESULTS = ROOT / "climate/results/ocean"


@pytest.mark.skipif(not (RESULTS / "checks_mitgcm.json").exists(), reason="no MITgcm checks recorded")
def test_recorded_mitgcm_checks_agree_with_the_solver():
    record = json.loads((RESULTS / "checks_mitgcm.json").read_text())
    if "channel" in record:
        assert abs(record["channel"]["model_over_analytic"]) < 0.005
    if "seiches_smythii" in record:
        assert all(abs(m["relative_difference"]) < 0.001 for m in record["seiches_smythii"]["modes"])
    for case in record.get("forced_fecunditatis", {}).get("results", []):
        f = case["fecunditatis_ratio"]
        assert f["mitgcm"] == pytest.approx(f["solver"], rel=0.01)
        assert case["rms_difference_over_rms_solver"] < 0.02


@pytest.mark.skipif(not (ROOT / "geography/products/atlas_28pct_4ppd.npz").exists(), reason="no atlas grid")
def test_solver_cells_map_one_to_one_onto_the_grids_water():
    from climate.ocean import checks
    grid, sea, (jj, ii), _ = checks.quarter_degree(checks.SMYTHII)
    assert (grid.depth[jj, ii] > 0).all() and len(set(zip(jj.tolist(), ii.tolist()))) == sea.cells
    assert int((grid.depth > 0).sum()) == sea.cells


@pytest.mark.skipif(not (RESULTS / "checks_schism.json").exists(), reason="no SCHISM checks recorded")
def test_recorded_schism_checks():
    record = json.loads((RESULTS / "checks_schism.json").read_text())
    if "channel" in record:
        assert abs(record["channel"]["model_over_analytic"]) < 0.005
    if "seiches_smythii_uniform" in record:
        assert all(abs(m["relative_difference"]) < 0.06 for m in record["seiches_smythii_uniform"]["modes"])
