"""The ocean models' set-up: lunar constants, grids, meshes, the tide's forcing files and the SCHISM patch.

None of these tests runs a model; those that need the pinned sources on the research drive skip without them.
"""
import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

from climate.ocean import builds, checks, forcing, lunar, mesh, mitgcm, stratified
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


def test_pickup_reader_gives_each_field_its_levels(tmp_path):
    # a pickup of two levels on a 4 by 3 grid: two whole-depth fields, then the surface's level
    records = np.arange(5 * 12, dtype=">f8").reshape(5, 3, 4)
    records.tofile(tmp_path / "pickup.0000000020.data")
    meta = (" nDims = [   2 ];\n dimList = [\n     4,    1,    4,\n     3,    1,    3\n ];\n dataprec = [ 'float64' ];\n"
            " nrecords = [     5 ];\n timeStepNumber = [         20 ];\n nFlds = [    3 ];\n fldList = {\n"
            " 'Theta   ' 'Salt    ' 'EtaN    '\n };\n")
    (tmp_path / "pickup.0000000020.meta").write_text(meta)
    out = mitgcm.read_pickup(tmp_path / "pickup.0000000020", 2)
    assert out["Theta"].shape == (2, 3, 4) and np.array_equal(out["Salt"], records[2:4])
    assert np.array_equal(out["EtaN"], records[4])
    with pytest.raises(ValueError):
        mitgcm.read_pickup(tmp_path / "pickup.0000000020", 3)


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
    assert "station_eta" not in s
    r = fecunditatis.Reducer(tmp_path / "s.npz", lat, lon, np.ones(3), np.array([1, 1, 2]), 1, np.array([2, 0]), ["a", "b"])
    r.add(t0, np.array([1.0, 3.0, 5.0]), np.array([0.1, 2.0, 0.2]))
    r.add(t0 + 10800.0, np.array([0.0, 1.0, 2.0]), np.array([0.3, 0.5, 0.1]))
    r.save()
    s = fecunditatis.Reducer(tmp_path / "s.npz", lat, lon, np.ones(3), np.array([1, 1, 2]), 1, np.array([2, 0]), ["a", "b"]).state
    assert np.array_equal(s["station_eta"], [[5.0, 1.0], [2.0, 0.0]]) and list(s["station_names"]) == ["a", "b"]


def test_points_take_the_mare_they_lie_in():
    from climate.ocean import fecunditatis
    from geography import nomenclature
    try:
        nomenclature.features()
    except (FileNotFoundError, OSError, KeyError):
        pytest.skip("no nomenclature archive")
    idx, names = fecunditatis.assign_maria(np.array([34.7, -7.8, 16.2, -60.0]), np.array([-14.9, 53.7, 59.1, 0.0]))
    assert [names[i] for i in idx[:3]] == ["Mare Imbrium", "Mare Fecunditatis", "Mare Crisium"] and idx[3] == -1
    assert fecunditatis.weighted_median(np.array([1.0, 2.0, 10.0]), np.array([1.0, 1.0, 5.0])) == 10.0


def test_numerics_comparison_reads_largest_relative_differences():
    from climate.ocean import strait_numerics
    ref = dict(t_h=[3.0, 6.0, 9.0], fecunditatis_level_m=[-0.01, -0.02, -0.04], western_channels_m_s=[0.1, 0.2, 0.4],
               crisium_passage_m_s=[0.1, 0.2, 0.5], largest_speed_m_s=0.5, p99_largest_speed_m_s=0.05,
               speed_max=np.array([0.1, 0.5]))
    case = dict(ref, t_h=[0.0, 3.0, 6.0, 9.0], fecunditatis_level_m=[0.0, -0.01, -0.02, -0.05],
                western_channels_m_s=[0.0, 0.1, 0.2, 0.44], crisium_passage_m_s=[0.0, 0.1, 0.2, 0.5],
                largest_speed_m_s=0.55, speed_max=np.array([0.1, 0.55]))
    out = strait_numerics.compare(case, ref)
    assert out["fecunditatis_level_m"] == pytest.approx(0.25) and out["western_channels_m_s"] == pytest.approx(0.1)
    assert out["crisium_passage_m_s"] == 0.0 and out["largest_speed_m_s"] == pytest.approx(0.1)
    assert out["rms_largest_speed"] == pytest.approx(np.sqrt(0.05 ** 2 / 2) / np.sqrt((0.1 ** 2 + 0.5 ** 2) / 2))


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


def test_schism_last_shorter_chunk_lines_up_its_files_and_ends_on_a_hotstart():
    from climate.ocean import schism
    assert schism.chunk_steps(0.0, 864000.0, 300.0, 10800.0, 864000.0) == (2880, 2880)          # a full chunk
    assert schism.chunk_steps(17280000.0, 864000.0, 300.0, 10800.0, 864000.0) == (2880, 2880)
    # the strait run's last 2.875 days from day 220: one file per snapshot, a hotstart at its last step only
    assert schism.chunk_steps(19008000.0, 248400.0, 300.0, 10800.0, 864000.0) == (36, 64188)
    with pytest.raises(ValueError):
        schism.chunk_steps(19008000.0, 249000.0, 300.0, 10800.0, 864000.0)                    # ends between two snapshots


def test_a_seas_own_frame_keeps_places_and_directions():
    frame = stratified.Frame(-47.8, -176.9)
    glat, glon = frame.to_grid(np.array([-47.8, -80.0, -14.0]), np.array([-176.9, 30.0, 170.0]))
    assert glat[0] == pytest.approx(0.0, abs=1e-9) and glon[0] == pytest.approx(0.0, abs=1e-9)
    lat, lon = frame.to_true(glat, glon)
    assert lat == pytest.approx([-47.8, -80.0, -14.0]) and np.mod(lon, 360.0) == pytest.approx([183.1, 30.0, 170.0])
    # a step to the true north of the centre is a step to the grid's north
    north, _ = frame.to_grid(np.array([-47.7]), np.array([-176.9]))
    assert north[0] == pytest.approx(0.1, abs=1e-6)
    # a wind keeps its speed everywhere and its direction at the centre
    u, v = frame.vector_to_grid(glat, glon, np.array([1.0, 3.0, 0.0]), np.array([0.0, 4.0, 2.0]))
    assert np.hypot(u, v) == pytest.approx([1.0, 5.0, 2.0])
    assert (u[0], v[0]) == pytest.approx((1.0, 0.0), abs=1e-9)
    # and points where a short true-eastward step lands on the grid
    step = 0.01
    for la, lo in ((-70.0, 120.0), (-20.0, -150.0)):
        g0, g1 = frame.to_grid(np.array([la]), np.array([lo])), frame.to_grid(np.array([la]), np.array([lo + step]))
        de = (g1[1][0] - g0[1][0]) * math.cos(math.radians(g0[0][0]))
        dn = g1[0][0] - g0[0][0]
        ue, vn = frame.vector_to_grid(g0[0], g0[1], np.array([1.0]), np.array([0.0]))
        assert (ue[0], vn[0]) == pytest.approx((de / math.hypot(de, dn), dn / math.hypot(de, dn)), abs=1e-4)


def test_levels_thicken_downwards_to_hold_the_depth():
    dz = stratified.levels(5576.0)
    assert dz[0] == 5.0 and dz.max() == 500.0 and np.all(np.diff(dz) >= 0)
    assert dz.sum() >= 5576.0 > dz[:-1].sum()
    assert dz[1] / dz[0] == pytest.approx(1.2)


def test_a_stratified_sea_carries_the_moon_and_reads_its_fluxes(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    grid = mitgcm.Grid(xg=np.arange(-2.0, 2.01, 1.0), yg=np.arange(-2.0, 2.01, 1.0), depth=np.full((4, 4), 60.0))
    flat = np.zeros((4, 4))
    sea = stratified.Sea(name="s", grid=grid, deltat=1800.0, delr=np.array([5.0, 6.0, 50.0]),
                         coriolis=dict(centre=flat - 3e-6, corner=flat - 3e-6, cosine=flat + 4e-6), theta0=27.0,
                         salt0=35.0, forcing=dict(ustress=np.zeros((3, 4, 4)), hflux=np.ones((3, 4, 4))),
                         runoff=flat + 1e-7, forcing_period_s=23400.0, sst_relax_s=0.0, output_s=86400.0)
    sea.write()
    data = (sea.folder / "data.template").read_text()
    assert f"gravity={lunar.GRAVITY:.12g}" in data and "rSphere=1737400" in data
    assert "selectCoriMap=3" in data and "eosType='JMD95Z'" in data and "selectP_inEOS_Zc=0" in data
    assert "tRef=3*27," in data and "sRef=3*35," in data
    assert "useRealFreshWaterFlux=.TRUE." in data and "nonlinFreeSurf=4" in data and "select_rStar=2" in data
    assert "implicitDiffusion=.TRUE." in data and "implicitViscosity=.TRUE." in data     # the closure needs both
    exf = (sea.folder / "data.exf").read_text()
    assert "repeatPeriod=70200," in exf and "ustressperiod=23400," in exf and "ustressStartTime=11700," in exf
    assert "hfluxfile='hflux.bin'" in exf and "runofffile='runoff.bin'" in exf and "climsst" not in exf
    assert "useGGL90=.TRUE." in (sea.folder / "data.pkg").read_text()
    # the old-water dye: 1 everywhere at the start, drawn to 0 in the top layer alone
    assert "usePTRACERS=.TRUE." in (sea.folder / "data.pkg").read_text()
    assert "PTRACERS_ref(1:3,1)=3*1.," in (sea.folder / "data.ptracers").read_text()
    assert "tauRelaxPTR(1)=86400," in (sea.folder / "data.rbcs").read_text()
    mask = np.fromfile(sea.folder / "rbcs_mask.bin", dtype=">f8").reshape(3, 4, 4)
    assert mask[0].all() and not mask[1:].any()
    assert np.fromfile(sea.folder / "fCoriC.bin", dtype=">f8").tolist() == [-3e-6] * 16
    assert np.fromfile(sea.folder / "hflux.bin", dtype=">f8").size == 3 * 16
    sea.write_diagnostics()
    diag = (sea.folder / "data.diagnostics").read_text()
    assert "fileName(1) = 'surface'" in diag and "levels(1,2) = 1.0" in diag and "'deep'" not in diag
    sea.deep_output_s = 864000.0
    sea.write_diagnostics()
    assert "'TRAC01  '" in (sea.folder / "data.diagnostics").read_text()
    # a run on mixed cell sizes: viscosity by the cell's size, and a start carried over from another run
    cold = np.full((3, 4, 4), 26.0)
    patch = stratified.Sea(name="p", grid=grid, deltat=900.0, delr=np.array([5.0, 6.0, 50.0]), coriolis=sea.coriolis,
                           visc_h=0.0, visc_grid=0.0012, start=dict(theta=cold, salt=cold + 9.0, dye=cold * 0.0))
    patch.write()
    data = (patch.folder / "data.template").read_text()
    assert "viscAhGrid=0.0012," in data and "viscAh=0," in data
    assert "hydrogThetaFile='start_theta.bin'" in data and "hydrogSaltFile='start_salt.bin'" in data
    assert "PTRACERS_initialFile(1)='start_dye.bin'" in (patch.folder / "data.ptracers").read_text()
    assert np.fromfile(patch.folder / "start_salt.bin", dtype=">f8").tolist() == [35.0] * 48
    assert "hydrogThetaFile" not in (sea.folder / "data.template").read_text()
    with pytest.raises(ValueError):
        stratified.Sea(name="b", grid=grid, deltat=1800.0, delr=np.array([5.0]), coriolis=sea.coriolis,
                       forcing=dict(ustress=np.zeros((3, 4, 4)), hflux=np.ones((2, 4, 4))),
                       forcing_period_s=1.0).inputs()


def test_a_stratified_sea_takes_the_tide_as_a_potential(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    grid = mitgcm.Grid(xg=np.arange(-2.0, 2.01, 1.0), yg=np.arange(-2.0, 2.01, 1.0), depth=np.full((4, 4), 60.0))
    flat = np.zeros((4, 4))
    coriolis = dict(centre=flat, corner=flat, cosine=flat)
    height = np.arange(3 * 16, dtype=float).reshape(3, 4, 4) * 0.01
    sea = stratified.Sea(name="t", grid=grid, deltat=1800.0, delr=np.array([5.0, 50.0]), coriolis=coriolis,
                         tide_m=height, tide_period_s=86400.0)
    sea.write()
    assert "momTidalForcing=.TRUE." in (sea.folder / "data.template").read_text()
    exf = (sea.folder / "data.exf").read_text()
    assert "tidePotFile='tidepot.bin'" in exf and "tidePotPeriod=86400," in exf and "tidePotStartTime=0," in exf
    assert "tidePotRepCycle=259200," in exf                       # its own cycle, not the fluxes' year
    assert f"exf_inscal_tidePot={-lunar.GRAVITY:.12g}" in exf       # the potential is -g times the height
    assert np.array_equal(np.fromfile(sea.folder / "tidepot.bin", dtype=">f8").reshape(3, 4, 4), height)
    assert json.loads((sea.folder / "experiment.json").read_text())["tide"] == dict(records=3, period_s=86400.0)
    still = stratified.Sea(name="u", grid=grid, deltat=1800.0, delr=np.array([5.0, 50.0]), coriolis=coriolis)
    still.write()
    assert "momTidalForcing" not in (still.folder / "data.template").read_text()
    assert "tidePot" not in (still.folder / "data.exf").read_text() and not (still.folder / "tidepot.bin").exists()


def test_a_stratified_sea_carries_tracers_with_their_own_sources(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    grid = mitgcm.Grid(xg=np.arange(-2.0, 2.01, 1.0), yg=np.arange(-2.0, 2.01, 1.0), depth=np.full((4, 4), 60.0))
    grid.depth[0, 0] = 0.0
    flat = np.zeros((4, 4))
    coriolis = dict(centre=flat - 3e-6, corner=flat - 3e-6, cosine=flat + 4e-6)
    delr = np.array([5.0, 6.0, 50.0])
    source = np.full((3, 4, 4), 2e-8)
    tracers = (stratified.Tracer("old", start=np.full((3, 4, 4), 0.7)),
               stratified.Tracer("used", start=0.0, top=0.0, top_s=86400.0, source=source))
    sea = stratified.Sea(name="o", grid=grid, deltat=1800.0, delr=delr, coriolis=coriolis, dye=False,
                         tracers=tracers, tracers_from=48)
    sea.write()
    ptr = (sea.folder / "data.ptracers").read_text()
    assert "PTRACERS_numInUse=2," in ptr and "PTRACERS_Iter0=48," in ptr
    assert "PTRACERS_initialFile(2)='start_tr02.bin'" in ptr and "PTRACERS_names(2)='used'" in ptr
    rbcs = (sea.folder / "data.rbcs").read_text()
    assert "relaxMaskFile(3)='rbcs_mask_tr01.bin'" in rbcs and "relaxMaskFile(4)='rbcs_mask_tr02.bin'" in rbcs
    assert "tauRelaxPTR(2)=86400," in rbcs and "relaxPtracerFile(2)='rbcs_target_tr02.bin'" in rbcs
    start = np.fromfile(sea.folder / "start_tr01.bin", dtype=">f8").reshape(3, 4, 4)
    assert start[1, 1, 1] == 0.7 and start[:, 0, 0].tolist() == [0.0] * 3        # land holds none
    mask = np.fromfile(sea.folder / "rbcs_mask_tr02.bin", dtype=">f8").reshape(3, 4, 4)
    target = np.fromfile(sea.folder / "rbcs_target_tr02.bin", dtype=">f8").reshape(3, 4, 4)
    assert mask[0, 1, 1] == 1.0 and target[0, 1, 1] == 0.0                    # the top layer is held at 0
    # below it, the relaxation package's pull -mask (c - target) / tau adds the source to within c / BEYOND
    for c in (0.0, 5.0, -3.0):
        pull = -mask[1:, 1, 1] * (c - target[1:, 1, 1]) / 86400.0
        assert np.allclose(pull, 2e-8, rtol=1e-5)
    assert not mask[:, 0, 0].any()
    assert json.loads((sea.folder / "experiment.json").read_text())["tracers"] == ["old", "used"]
    sea.deep_output_s = 864000.0
    sea.extra_lists = ((("UVEL", "VVEL"), "daily", 86400.0, 1),)
    sea.write_diagnostics()
    diag = (sea.folder / "data.diagnostics").read_text()
    assert "'TRAC01  ','TRAC02  '" in diag and "fileName(4) = 'daily'" in diag and "levels(1,4) = 1.0" in diag
    with pytest.raises(ValueError):
        stratified.Sea(name="d", grid=grid, deltat=1800.0, delr=delr, coriolis=coriolis, tracers=tracers)


def test_a_stratified_sea_takes_internal_wave_power_and_the_floors_heat(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    grid = mitgcm.Grid(xg=np.arange(-2.0, 2.01, 1.0), yg=np.arange(-2.0, 2.01, 1.0), depth=np.full((4, 4), 60.0))
    grid.depth[0, 0] = 0.0
    flat = np.zeros((4, 4))
    coriolis = dict(centre=flat - 3e-6, corner=flat - 3e-6, cosine=flat + 4e-6)
    sea = stratified.Sea(name="i", grid=grid, deltat=1800.0, delr=np.array([5.0, 55.0]), coriolis=coriolis,
                         idemix=dict(tidal=flat + 1e-6, wind=flat + 2e-4), geothermal=flat + 0.011, diff_r=1.4e-7)
    sea.write()
    ggl90 = (sea.folder / "data.ggl90").read_text()
    first, second = ggl90.split("&GGL90_PARM02")                                # IDEMIX reads its own group
    assert "useIDEMIX=.TRUE." in first and "IDEMIX_wind_file='idemix_wind.bin'" in second
    assert "geothermalFile='geothermal.bin'" in (sea.folder / "data.template").read_text()
    assert "diffKrT=1.4e-07" in (sea.folder / "data.template").read_text()
    wind = np.fromfile(sea.folder / "idemix_wind.bin", dtype=">f8").reshape(4, 4)
    assert wind[0, 0] == 0.0 and wind[1, 1] == 2e-4                            # nothing over land
    record = json.loads((sea.folder / "experiment.json").read_text())
    assert record["idemix_mw"]["wind"] == pytest.approx(2e-4 * grid.area[grid.wet].sum() / 1e6)
    plain = stratified.Sea(name="p", grid=grid, deltat=1800.0, delr=np.array([5.0, 55.0]), coriolis=coriolis)
    plain.write()
    assert "IDEMIX" not in (plain.folder / "data.ggl90").read_text()
    assert "idemix_mw" not in json.loads((plain.folder / "experiment.json").read_text())
    assert stratified.options_with("#undef ALLOW_GGL90_IDEMIX\n", "ALLOW_GGL90_IDEMIX", "x") == \
        "#define ALLOW_GGL90_IDEMIX\n"
    with pytest.raises(RuntimeError):
        stratified.options_with("#define ALLOW_GGL90_IDEMIX\n", "ALLOW_GGL90_IDEMIX", "x")


def test_a_restart_carries_the_parents_state_and_starts_the_waves_at_rest(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, "RUNS", tmp_path)
    parent = tmp_path / "parent"
    parent.mkdir()
    tke = np.arange(2 * 16, dtype=">f8").tobytes()
    (parent / "pickup.0000000096.data").write_bytes(b"state")
    (parent / "pickup.0000000096.meta").write_text("meta")
    (parent / "pickup_ggl90.0000000096.data").write_bytes(tke)
    (parent / "pickup_ggl90.0000000096.meta").write_text(" nrecords = [          1 ];\n")
    grid = mitgcm.Grid(xg=np.arange(-2.0, 2.01, 1.0), yg=np.arange(-2.0, 2.01, 1.0), depth=np.full((4, 4), 60.0))
    flat = np.zeros((4, 4))
    coriolis = dict(centre=flat - 3e-6, corner=flat - 3e-6, cosine=flat + 4e-6)
    sea = stratified.Sea(name="r", grid=grid, deltat=1800.0, delr=np.array([5.0, 55.0]), coriolis=coriolis,
                         dye=False, idemix=dict(tidal=flat, wind=flat))
    sea.write(fresh=True)
    sea.restart_from(parent, 96, build="b")
    assert (sea.folder / "pickup.0000000096.data").read_bytes() == b"state"
    data = (sea.folder / "pickup_ggl90.0000000096.data").read_bytes()
    assert data[:len(tke)] == tke and data[len(tke):] == bytes(len(tke))       # the waves' energy starts at 0
    assert "nrecords = [          2 ]" in (sea.folder / "pickup_ggl90.0000000096.meta").read_text()
    assert sea.pickups() == [96]
    assert (sea.folder / "setup.sha256").read_text().split()[0] == sea.setup_digest("b")
    late = stratified.Sea(name="l", grid=grid, deltat=1800.0, delr=np.array([5.0, 55.0]), coriolis=coriolis,
                          dye=False, tracers=(stratified.Tracer("new"),), tracers_from=48)
    late.write(fresh=True)
    with pytest.raises(ValueError):
        late.restart_from(parent, 96, build="b")


def test_a_stratified_sea_compiles_out_the_bulk_formulas(tmp_path):
    for sub in ("pkg/exf", "model/inc"):
        (tmp_path / sub).mkdir(parents=True)
    (tmp_path / "pkg/exf/EXF_OPTIONS.h").write_text(
        "#define ALLOW_ATM_TEMP\n#define ALLOW_ATM_WIND\n#define ALLOW_DOWNWARD_RADIATION\n#ifdef ALLOW_ATM_TEMP\n"
        "# define ALLOW_BULKFORMULAE\n#endif\n#define ALLOW_RUNOFF\n#define ALLOW_CLIMSST_RELAXATION\n"
        "#undef EXF_ALLOW_TIDES\n")
    text = stratified.exf_options(tmp_path)
    assert "#undef ALLOW_ATM_TEMP" in text and "#undef ALLOW_ATM_WIND" in text and "#define ALLOW_RUNOFF" in text
    assert not re.search(r"^#define\s+ALLOW_ATM_TEMP", text, re.M)
    assert "#undef EXF_ALLOW_TIDES" in text
    assert "#define EXF_ALLOW_TIDES" in stratified.exf_options(tmp_path, tides=True)
    (tmp_path / "model/inc/CPP_OPTIONS.h").write_text(
        "#undef SHORTWAVE_HEATING\n#undef ALLOW_SOLVE4_PS_AND_DRAG\n#define INCLUDE_IMPLVERTADV_CODE\n"
        "#undef NONLIN_FRSURF\n")
    text = stratified.model_options(tmp_path)
    assert "#define SHORTWAVE_HEATING" in text and "#define NONLIN_FRSURF" in text
