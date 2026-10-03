"""Check atmospheric bytes and numerical sensitivity of the weather-forced sea."""
from __future__ import annotations

import json
from datetime import timedelta
from pathlib import Path
import struct
import re

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot_checks import spectrum_diagnostics
from climate.waves.weather import HERE, RUNS, forcing
from climate.waves.pilot import load_basin


def read_grid_records(path, codes):
    """Independent streaming reader of the pinned model's little-endian records."""
    result = {c: [] for c in codes}
    times = {c: [] for c in codes}
    with path.open("rb") as stream:
        def record():
            marker = stream.read(4)
            if not marker:
                return None
            if len(marker) != 4:
                raise ValueError("Truncated record marker")
            size = struct.unpack("<i", marker)[0]
            if size < 0 or size > 10**7:
                raise ValueError("Unsupported atmospheric record size or byte order")
            data, end = stream.read(size), stream.read(4)
            if len(data) != size or end != marker:
                raise ValueError("Incomplete atmospheric record")
            return data
        while (header := record()) is not None:
            if len(header) != 32:
                raise ValueError("Invalid atmospheric header")
            code, level, _, _, nx, ny, step, _ = struct.unpack("<8i", header)
            values = record()
            if values is None:
                raise ValueError("Missing atmospheric data record")
            if code in codes:
                if level != 0 or len(values) != nx*ny*4:
                    raise ValueError("Expected a surface grid record")
                result[code].append(np.frombuffer(values, dtype="<f4").reshape(ny,nx))
                times[code].append(step)
    return {c: (np.asarray(times[c]), np.asarray(result[c])) for c in codes}


def atmospheric_checks(run_root=RUNS):
    records = [json.loads((run_root / name / "run.json").read_text()) for name in ("snapshots", "control")]
    for record, name in zip(records, ("snapshots", "control")):
        if record["state"] != "complete":
            raise ValueError("Atmospheric run remains incomplete")
        if sha256(run_root / name / "producer.py") != record["producer_sha256"]:
            raise ValueError("Archived atmospheric runner source changed")
        for filename, expected in record["output_sha256"].items():
            if sha256(run_root / name / filename) != expected:
                raise ValueError("Atmospheric run bytes differ from their manifest")
        for filename, expected in record["source_inputs"].items():
            if sha256(Path(filename)) != expected:
                raise ValueError("Source climate input changed")
    preserved = {name: records[0]["output_sha256"][name] == records[1]["output_sha256"][name]
                 for name in ("plasim_status", "plasim_output")}
    if not all(preserved.values()):
        raise ValueError("Snapshot recording changes the atmospheric control")
    codes = {139:"surface_temperature_K",159:"friction_velocity_cubed_m3_s3",172:"land_fraction",
             173:"roughness_m",180:"stress_eastward_Pa",181:"stress_northward_Pa",210:"sea_ice_fraction"}
    raw = read_grid_records(run_root / "snapshots/plasim_snapshot", codes)
    with np.load(run_root / "atmosphere.npz") as d:
        meta = json.loads(d["metadata"].item())
        if sha256(run_root / "snapshots/export_producer.py") != meta["producer"]["source_sha256"]:
            raise ValueError("Archived atmospheric exporter source changed")
        count = 0
        for code, key in codes.items():
            times, values = raw[code]
            np.testing.assert_array_equal(times * records[0]["timestep_s"], d["restart_elapsed_s"])
            np.testing.assert_array_equal(values, d[key])
            count += values.size
        mask = (d["land_fraction"] < .5) & (d["sea_ice_fraction"] < .01)
        tx, ty, u, v = [d[k] for k in ("stress_eastward_Pa","stress_northward_Pa","eastward_lowest_m_s","northward_lowest_m_s")]
        cosine = (tx*u+ty*v)/np.maximum(np.hypot(tx,ty)*np.hypot(u,v),1e-30)
        mask &= np.hypot(tx,ty) > 1e-8
        direction_error = float(np.max(np.abs(cosine[mask]-1)))
        if direction_error > 1e-5:
            raise ValueError("Recovered wind direction differs from the signed surface stress")
    src = ROOT / "climate/gcm/.venv/lib"
    source_paths = sorted(src.glob("python*/site-packages/exoplasim/plasim/src/*.f90"))
    audited = {p.name: sha256(p) for p in source_paths
               if p.name in ("plasim.f90","outmod.f90","fluxmod.f90","seamod.f90",
                             "icemod.f90","oceanmod.f90","glaciermod.f90")}
    optional_restart={name:(Path(records[0]["source_restart"]).parent/name).is_file()
                      for name in ("restart_dsnow","restart_xsnow")}
    return dict(control_equal_bytes=preserved, original_atmospheric_inputs_unchanged=True,
                independent_surface_values_checked=count, independent_surface_values_exact=True,
                wind_stress_direction_cosine_max_error=direction_error,
                snapshot_count=meta["sample_count"], duration_s=meta["duration_s"],
                surface_closure_relative_max=meta["stress_closure_relative_max"],
                source_audit_files=audited, source_audit_scope="Installed source read alongside the archived executable; executable hash retained in each run manifest.",
                optional_source_snow_overrides_present=optional_restart,
                run_manifests={name: sha256(run_root / name / "run.json") for name in ("snapshots","control")})


def load_case(name):
    path = RUNS / name / "product.json"
    data = json.loads(path.read_text())
    if data["schema"] != "terluna.climate.weather-wave-case/1":
        raise ValueError("Unsupported wave case")
    for filename, expected in {**data["run"]["identity"]["inputs"], **data["run"]["output_sha256"]}.items():
        if sha256(path.parent / filename) != expected:
            raise ValueError(f"Changed wave input or output: {name}/{filename}")
    return data


def compare(a, b, start_hour=24, offset_b=0):
    """Matched times/locations; relative error uses the refined or longer-history case."""
    coarse = np.asarray(a["values"])[start_hour:]
    fine = np.asarray(b["values"])[start_hour+offset_b:start_hour+offset_b+len(coarse)]
    if coarse.shape != fine.shape:
        raise ValueError("Comparison shapes differ")
    np.testing.assert_allclose(coarse[...,0]+offset_b*3600,fine[...,0],atol=.1)
    np.testing.assert_array_equal(coarse[...,1:3],fine[...,1:3])
    np.testing.assert_array_equal(coarse[...,6],fine[...,6])
    mask = (coarse[...,3] >= .1) & (fine[...,3] >= .1)
    result = dict(first_compared_hour=start_hour,last_compared_hour=a["hours"],
                  compared_cases=[a.get("name"),b.get("name")],
                  additional_history_hours=offset_b, paired_wet_samples=int(mask.sum()), minimum_hs_m=.1)
    for name, col in (("hs",3),("tm01",5)):
        difference = np.abs(coarse[...,col]-fine[...,col])
        relative = difference[mask]/fine[...,col][mask]
        result[name+"_relative_max"] = float(relative.max())
        result[name+"_relative_p95"] = float(np.percentile(relative,95))
        result[name+("_absolute_max_m" if name=="hs" else "_absolute_max_s")] = float(difference[mask].max())
        all_relative = np.where(mask,difference/np.maximum(fine[...,col],1e-20),-1)
        t,point = np.unravel_index(np.argmax(all_relative),mask.shape)
        result[name+"_worst_pair"] = dict(coarse=coarse[t,point].tolist(),reference=fine[t,point].tolist())
    strong = mask & (coarse[...,3] >= 1) & (fine[...,3] >= 1)
    result["strong_sea_diagnostic"] = dict(paired_samples=int(strong.sum()),minimum_hs_m=1.,
        hs_relative_max=float(np.max(np.abs(coarse[...,3][strong]/fine[...,3][strong]-1))) if strong.any() else None,
        tm01_relative_max=float(np.max(np.abs(coarse[...,5][strong]/fine[...,5][strong]-1))) if strong.any() else None)
    result["within_five_percent"] = result["hs_relative_max"] < .05 and result["tm01_relative_max"] < .05
    return result


def shared_forcing_check(a, b, offset_b=0):
    """Confirm that numerical/history comparisons see the same wind-file values."""
    af,bf=a["forcing"],b["forcing"]
    if (af["atmosphere_sha256"] != bf["atmosphere_sha256"] or
            af["air_density_kg_m3"] != bf["air_density_kg_m3"] or
            af["first_snapshot_index"] != bf["first_snapshot_index"]+offset_b//3):
        raise ValueError("Comparison atmospheric origins differ")
    shape=np.asarray(a["wet"]).shape
    winds=[np.loadtxt(RUNS/c["name"]/"wind.dat").reshape(-1,2,*shape) for c in (a,b)]
    difference=np.max(np.abs(winds[0]-winds[1][offset_b*4:offset_b*4+len(winds[0])]))
    if difference > 1e-8:
        raise ValueError("Numerical comparison changes the forcing")
    return dict(compared_cases=[a["name"],b["name"]],common_forcing_max_absolute_difference_m_s=float(difference),passed=True)


def describe_case(case, display_start=0):
    data = np.asarray(case["values"])[display_start:]
    hs = data[...,3]
    weights = np.cos(np.radians(data[0,:,2]))
    # The first 24 hours of the shorter forcing window are set aside uniformly.
    retained = data[24:]
    peak_time, peak_node = np.unravel_index(np.argmax(retained[...,3]),retained[...,3].shape)
    peak_time += 24
    offshore = case["offshore_index"]
    last = (ORIGIN+timedelta(hours=case["hours"])).strftime("%Y%m%d.%H%M%S")
    spectrum = spectrum_diagnostics(RUNS / case["name"] / "offshore.spc", last)
    basin = load_basin(stride=case["stride"])
    winds = forcing(RUNS / "atmosphere.npz", basin, case["forcing"]["air_density_kg_m3"], case["hours"])
    offshore_wind = np.column_stack((winds["u"][::4,basin["wet"]][:,offshore],
                                    winds["v"][::4,basin["wet"]][:,offshore]))[display_start:]
    fraction = np.average(hs >= 1,axis=1,weights=weights)
    r = dict(name=case["name"],step_s=case["step_s"],run_hours=case["hours"],
             record_time_origin="Raw rows retain seconds since this case's start. Subtract displayed_from_run_hour*3600 to place them within the selected week.",
             displayed_from_run_hour=display_start,analysis_first_display_hour=24,
             max_hs_m=float(retained[...,3].max()),peak_record=data[peak_time,peak_node].tolist(),
             peak_hour=float(data[peak_time,0,0]/3600-display_start),
             max_area_fraction_hs_ge_1m=float(fraction[24:].max()),
             offshore_max_hs_m=float(retained[:,offshore,3].max()),
             offshore_mean_period_range_s=[float(retained[:,offshore,5].min()),float(retained[:,offshore,5].max())],
             elapsed_wall_s=case["run"]["elapsed_wall_s"], source_product_sha256=sha256(RUNS/case["name"]/"product.json"),
             hourly_time_hours=(data[:,0,0]/3600-display_start).tolist(),
             hourly_max_hs_m=hs.max(axis=1).tolist(), hourly_area_fraction_hs_ge_1m=fraction.tolist(),
             offshore_records=data[:,offshore].tolist(),offshore_equivalent_wind_m_s=offshore_wind.tolist(),
             peak_map_records=data[peak_time].tolist(),node_max_hs_m=retained[...,3].max(axis=0).tolist(),
             spectrum_at_final_offshore=spectrum,
             stress_coupling_ustar_absolute_max_m_s=case["stress_coupling_ustar_absolute_max_m_s"],
             run_output_sha256=case["run"]["output_sha256"])
    return r


def main():
    a,b,c,d = [load_case(name) for name in ("basin_dt150","basin_dt75","basin_lead48_dt150","basin_lead96_dt150")]
    dt, lead = compare(a,b), compare(a,c,offset_b=48)
    input_text=(RUNS/a["name"]/"INPUT").read_text()
    gravity,density=map(float,re.search(r"SET GRAV ([\d.eE+-]+) RHO ([\d.eE+-]+)",input_text).groups())
    gamma=float(re.search(r"BREAKING CONSTANT [\d.eE+-]+ ([\d.eE+-]+)",input_text)[1])
    basin=load_basin(stride=a["stride"])
    result = dict(schema="terluna.climate.weather-wave-study/1",
                  evidence="A selected seven-day Smythii–Marginis wave simulation driven by a global atmospheric restart's three-hourly surface stress, with timestep and additional-history comparisons.",
                  reading_rule="Significant heights and periods describe this selected model episode. The stress-equivalent input preserves GCM momentum transfer through SWAN's Wu drag law. The 1-degree basin grid resolves broad exposure; the recorded coastal convergence failures still apply. Occurrence rates across months and seasons require longer, finer ensembles.",
                  units=dict(hs_m="m",period_s="s",direction_deg="degrees counterclockwise from east, waves travelling toward",time_hours="Earth hours since start of selected seven-day window"),
                  columns=a["columns"],
                  undefined_diagnostics="At Hs=0, SWAN's period and direction sentinels (-9 s and -999 degrees) are retained in the raw rows; derived plots mask those entries.",
                  area_weighting="Cos(latitude) weights on the 1-degree wet-node mask; area fractions describe the resolved basin grid.",
                  producer=dict(domain="climate",files={str(p.relative_to(ROOT)):sha256(p) for p in
                    (Path(__file__),HERE/"weather.py",ROOT/"climate/gcm/wave_snapshots.py",ROOT/"shared/constants.json")}),
                  atmosphere=atmospheric_checks(), forcing=a["forcing"],
                  model_settings=dict(gravity_m_s2=gravity,water_density_kg_m3=density,
                      depth_breaking_index=gamma,frequency_intervals=a["frequency_intervals"],
                      executable_sha256=a["run"]["identity"]["executable_sha256"],
                      source_input_sha256=a["run"]["identity"]["inputs"]["INPUT"],
                      physics="Komen/Wu/AGROW with the recorded gravity and air-density patches; DIA, whitecapping and depth breaking"),
                  geography=dict(atlas_sha256=a["atlas_sha256"],atlas_body=874,
                      sea_level_m=basin["metadata"]["sea_level_m"],wave_grid_spacing_deg=float(basin["lon"][1]-basin["lon"][0]),
                      terrain="Flooded LOLA/GRAIL atlas ground; finite dry-node elevations and other disconnected water bodies held dry"),
                  basin_grid=dict(longitude_deg=a["longitude_deg"],latitude_deg=a["latitude_deg"],wet=a["wet"]),
                  numerical_checks=dict(timestep_150_to_75_s=dt,additional_48h_history=lead,
                      prior_history_48_to_96h=compare(c,d,start_hour=72,offset_b=48),
                      shared_forcing=[shared_forcing_check(a,b),shared_forcing_check(a,c,48),shared_forcing_check(c,d,48)],
                      history_windows_48_to_96h={str(h):compare(c,d,start_hour=48+h,offset_b=48) for h in (0,48,72,96)}),
                  displayed_case="basin_lead96_dt150",
                  cases=[describe_case(a),describe_case(b),describe_case(c,48),describe_case(d,96)],
                  limitations=["Selection targets the largest basin-mean stress in one 30-Earth-day record; occurrence statistics need a representative longer sample.",
                               "The atmospheric grid supplies 16 contributing open-water cells. Coastal wind jets, gusts and sub-three-hour events remain unresolved.",
                               "Initial-state sensitivity is assessed with 48 and 96 hours of preceding forcing. The extra 96-hour run was added after the 48-hour extension revealed large initial-state effects.",
                               "Timestep refinement uses the initially calm seven-day episode; combining finer timesteps with the longer histories remains a further check.",
                               "The numerical comparisons cover timesteps and prior history on the 1-degree basin grid. Geographic and spectral refinements for this weather sequence remain open.",
                               "GCM surface stress, SWAN drag, wave growth and dissipation use Earth-derived closures at lunar gravity. Coupling is one-way, with fixed SWAN air density.",
                               "The wave calculation omits ocean currents, time-varying water level, bottom friction and resolved shoreline run-up."])
    path = HERE / "results/weather.json"
    path.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(json.dumps({"cases":[{k:c[k] for k in ("name","max_hs_m","peak_hour","max_area_fraction_hs_ge_1m","offshore_max_hs_m")} for c in result["cases"]],"checks":result["numerical_checks"]},indent=2))


if __name__ == "__main__":
    main()
