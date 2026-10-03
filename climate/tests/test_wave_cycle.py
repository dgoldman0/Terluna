"""Exact reporting boundaries, time occupancy and explicit forcing histories."""
import json

import numpy as np
import pytest

from climate.waves.cycle import CYCLE_HOURS, exceedance, summarize, window
from climate.waves.cycle_checks import read_spectral_series, refinement_window
from climate.waves.weather import drag_coefficient, forcing


def test_spinup_is_excluded_with_fractional_cycle_boundaries():
    # A large spin-up event supplies zero reporting time after hour 2.25.
    times = np.arange(7.)
    heights = np.array([20.,10.,0.,1.,2.,3.,4.])
    t, h = window(times, heights, 2.25, 4.5)
    np.testing.assert_array_equal(t, [2.25,3.,4.,4.5])
    np.testing.assert_array_equal(h, [.25,1.,2.,2.5])
    stats = exceedance(t,h,1.5)
    assert stats["hours"] == 1
    assert stats["longest_episode_hours"] == 1
    assert stats["intervals_hours"] == [[3.5,4.5]]
    assert not stats["begins_above"] and stats["ends_above"]


def test_crossing_durations_and_censored_episodes_are_distinct():
    t = np.array([0.,2.,4.,6.,8.])
    h = np.array([2.,0.,2.,0.,2.])
    stats = exceedance(t,h,1.)
    assert stats["intervals_hours"] == [[0.,1.],[3.,5.],[7.,8.]]
    assert stats["hours"] == 4 and stats["fraction"] == .5
    assert stats["longest_episode_hours"] == 2
    assert stats["episode_count"] == 3
    assert stats["begins_above"] and stats["ends_above"]
    # Endpoint contact carries zero duration; a plateau at threshold counts.
    assert exceedance([0,1,2],[0,1,0],1)["episode_count"] == 0
    assert exceedance([0,1,2,3],[0,1,1,0],1)["hours"] == 1


def test_duration_is_invariant_under_linear_resampling():
    t, h = np.array([0.,1.5,4.,8.]), np.array([0.,3.,1.,0.])
    dense = np.sort(np.r_[t,np.linspace(0,8,257)])
    dense = np.unique(dense)
    for threshold in (0.,.5,1.,2.,4.):
        coarse = exceedance(t,h,threshold)
        fine = exceedance(dense,np.interp(dense,t,h),threshold)
        assert fine["hours"] == pytest.approx(coarse["hours"],abs=1e-12)
        assert fine["longest_episode_hours"] == pytest.approx(coarse["longest_episode_hours"],abs=1e-12)
        assert fine["episode_count"] == coarse["episode_count"]


def test_incomplete_or_corrupt_reporting_histories_are_rejected():
    with pytest.raises(ValueError,match="complete reporting"):
        window([0,1,2],[0,1,2],1,3)
    with pytest.raises(ValueError,match="increasing"):
        window([0,1,1],[0,1,2],0,1)
    with pytest.raises(ValueError,match="finite"):
        exceedance([0,1],[0,np.nan],1)


@pytest.mark.parametrize("external_storage", [False, True], ids=["checkout", "external-drive"])
def test_cycle_product_excludes_spinup_from_every_height_statistic(monkeypatch, tmp_path, external_storage):
    values = np.zeros((7,2,9))
    values[...,0] = np.arange(7)[:,None]*3600
    values[...,1] = [80,81]
    values[...,2] = [0,0]
    values[...,3] = np.array([20,10,0,1,2,3,4])[:,None]
    values[...,4:6] = 10
    values[...,6] = 100
    values[...,7] = np.array([180,180,350,10,30,50,70])[:,None]
    case = dict(schema="terluna.climate.weather-wave-case/1",name="analytic",hours=6,
                stride=4,step_s=150,frequency_intervals=48,values=values.tolist(),
                wet=[[True,True]],longitude_deg=[80,81],latitude_deg=[0],offshore_index=0,
                columns=list(range(9)),atlas_sha256="analytic",run=dict(identity=dict(inputs={}),output_sha256={}),
                forcing=dict(first_snapshot_index=0,atmosphere_sha256="analytic",
                    atmosphere_metadata=dict(first_snapshot_restart_step=5,
                        source_configuration=dict(model=dict(timestep_min=30,sun_clock="synodic"),
                            planet=dict(rotationperiod=CYCLE_HOURS/48,year=CYCLE_HOURS/24)))))
    # Read an analytical history through the same checkout path with either
    # local storage or a symlink to a directory outside that checkout.
    root = tmp_path/"checkout"
    runs = root/"research/runs/waves"
    runs.parent.mkdir(parents=True)
    if external_storage:
        storage = tmp_path/"drive/wave-runs"
        storage.mkdir(parents=True)
        runs.symlink_to(storage, target_is_directory=True)
    else:
        runs.mkdir()
    path = runs/"analytic/product.json"
    path.parent.mkdir()
    path.write_text(json.dumps(case))
    here = root/"climate/waves"
    monkeypatch.setattr("climate.waves.cycle.ROOT", root)
    monkeypatch.setattr("climate.waves.cycle.HERE", here)
    monkeypatch.setattr("climate.waves.cycle.__file__", str(here/"cycle.py"))
    monkeypatch.chdir(root)
    monkeypatch.setattr("climate.waves.cycle.sha256",lambda _:"analytic")
    monkeypatch.setattr("climate.waves.cycle.spectrum_diagnostics",lambda *_:dict(passed=True))
    result=summarize(path.relative_to(root),cycle_hours=2.25)
    assert result["inputs"]["case_path"] == "research/runs/waves/analytic/product.json"
    assert result["statistics"]["max_hs_m"] == 2.5
    assert result["statistics"]["area_time_mean_hs_m"] == pytest.approx(1.375)
    threshold=next(t for t in result["statistics"]["thresholds"] if t["hs_m"]==2)
    assert threshold["node_hours"] == [.5,.5]
    assert threshold["offshore"]["hours"] == .5
    # The reporting boundary lies between directions 350 and 10 degrees.
    assert result["series"]["offshore"][0][7] == pytest.approx(355)
    case["forcing"]["atmosphere_metadata"]["source_configuration"]["model"]["sun_clock"]="rotation"
    path.write_text(json.dumps(case))
    with pytest.raises(ValueError,match="synodic solar clock"):
        summarize(path.relative_to(root),cycle_hours=2.25)


def test_explicit_forcing_uses_the_requested_dates_and_same_overlap(tmp_path):
    shape = (33,3,3)
    tx = np.broadcast_to(np.linspace(.001,.01,33)[:,None,None],shape)
    ty = np.full(shape,-.002)
    product = tmp_path/"atmosphere.npz"
    np.savez(product,time_s=np.arange(33)*10800,
             latitude_deg=[2.,1.,0.],longitude_deg=[0.,1.,2.],
             land_fraction=np.zeros(shape),sea_ice_fraction=np.zeros(shape),
             stress_eastward_Pa=tx,stress_northward_Pa=ty,
             stress_closure_density_kg_m3=np.full(shape,1.4),
             metadata=json.dumps(dict(schema="terluna.climate.gcm-wave-snapshots/1",radius_m=1e6)))
    basin=dict(lon=np.array([.5,1.5]),lat=np.array([.5,1.5]),wet=np.ones((2,2),bool))
    full = forcing(product,basin,1.4,96,start_hour=0)
    later = forcing(product,basin,1.4,48,start_hour=48)
    np.testing.assert_allclose(later["u"],full["u"][192:],atol=1e-14,rtol=0)
    np.testing.assert_allclose(later["v"],full["v"][192:],atol=1e-14,rtol=0)
    assert later["record"]["first_snapshot_index"] == 16
    assert later["record"]["end_day_from_first_snapshot"] == 4
    speed=np.hypot(full["u"][::12],full["v"][::12])
    actual=1.4*drag_coefficient(speed)*speed*full["u"][::12]
    np.testing.assert_allclose(actual[:,0,0],tx[:,0,0],atol=2e-17)
    for start in (-3,1,51):
        with pytest.raises(ValueError,match="Explicit forcing"):
            forcing(product,basin,1.4,48,start_hour=start)


def test_spectral_history_keeps_calm_times_locations_and_variance(tmp_path):
    header = """SWAN 1
TIME
1
LONLAT
2
80 0
81 1
AFREQ
3
.1
.2
.4
QUANT
3
VaDens
m2/Hz
-99
CDIR
degr
-999
DSPRDEGR
degr
-9
"""
    calm="\n".join(["-99 -999 -9"]*3)+"\n"
    waves="1 0 20\n2 10 20\n1 20 20\n"
    path=tmp_path/"history.spc"
    text=header+"20000101.000000\nLOCATION 1\n"+calm+"LOCATION 2\n"+waves
    text+="20000101.030000\nLOCATION 1\n"+waves+"LOCATION 2\n"+calm
    path.write_text(text)
    times,points,freq,density=read_spectral_series(path)
    np.testing.assert_array_equal(times,[0,3])
    np.testing.assert_array_equal(points,[[80,0],[81,1]])
    np.testing.assert_array_equal(density,[[[0,0,0],[1,2,1]],[[1,2,1],[0,0,0]]])
    assert np.trapezoid(density[0,1],freq) == pytest.approx(.45)
    path.write_text(text.replace("LOCATION 2","LOCATION 1",1))
    with pytest.raises(ValueError,match="reordered"):
        read_spectral_series(path)
    path.write_text(text.replace("1 0 20","-1 0 20",1))
    with pytest.raises(ValueError,match="Negative"):
        read_spectral_series(path)


def test_refinement_selection_uses_second_cycle_even_when_spinup_is_stronger():
    hours=np.arange(0,1441,3)
    stress=np.ones(len(hours))
    stress[hours==240]=100
    stress[hours==1008]=5
    selected=refinement_window(dict(monthly_time_days=(hours/24).tolist(),
                                    monthly_basin_mean_stress_Pa=stress.tolist()))
    assert selected["peak_basin_mean_stress_hour"]==1008
    assert CYCLE_HOURS < selected["analysis_start_hour"] < selected["analysis_end_hour"] < 2*CYCLE_HOURS
    assert selected["analysis_start_hour"]-selected["run_start_hour"]==120
