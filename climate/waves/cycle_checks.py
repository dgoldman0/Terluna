"""Verify cycle clocks, evolving spectra and matched wave refinements."""
from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
import re

import numpy as np

from climate.waves.cycle import CYCLE_HOURS, RUNS, exceedance, read_case
from climate.waves.model import ORIGIN, sha256
from climate.waves.weather_checks import compare


def refinement_window(forcing):
    """Choose the reporting cycle's strongest forcing week and its preceding days."""
    hours=np.asarray(forcing["monthly_time_days"])*24
    stress=np.asarray(forcing["monthly_basin_mean_stress_Pa"])
    selected=(hours>=CYCLE_HOURS)&(hours<=2*CYCLE_HOURS)
    if not selected.any() or hours[-1] < 2*CYCLE_HOURS:
        raise ValueError("The atmospheric record must cover both solar cycles")
    peak=np.flatnonzero(selected)[np.argmax(stress[selected])]
    week=int(min(max(hours[peak]-24,np.ceil(CYCLE_HOURS/3)*3),
                 np.floor((2*CYCLE_HOURS-168)/3)*3))
    return dict(peak_basin_mean_stress_hour=float(hours[peak]),peak_basin_mean_stress_Pa=float(stress[peak]),
                analysis_start_hour=week+24,analysis_end_hour=week+168,
                run_start_hour=week-96,run_hours=264,preceding_weather_hours=96,
                selection="Seven-day window containing the largest basin-mean stress in the reporting cycle, with 96 preceding hours and the first selected day excluded")


def read_spectral_series(path):
    """Read the recorded SPEC1D variance at ordered times and wet locations."""
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith("$")]
    index = next(i for i,s in enumerate(lines) if s.split()[0] == "LONLAT")
    points = int(lines[index+1].split()[0])
    xy = np.array([[float(v) for v in row.split()] for row in lines[index+2:index+2+points]])
    index = next(i for i,s in enumerate(lines) if s.split()[0] == "AFREQ")
    n = int(lines[index+1].split()[0])
    frequency = np.array([float(row.split()[0]) for row in lines[index+2:index+2+n]])
    index = next(i for i,s in enumerate(lines) if s.split()[0] == "QUANT")
    if int(lines[index+1].split()[0]) != 3 or lines[index+2].split()[0] != "VaDens":
        raise ValueError("Expected variance, mean direction and directional spread")
    exception = float(lines[index+4].split()[0])
    if (xy.shape != (points,2) or n < 3 or frequency[0] <= 0
            or np.any(np.diff(frequency) <= 0) or not np.isfinite(frequency).all() or not np.isfinite(xy).all()):
        raise ValueError("Invalid spectral grid")
    times, density = [], []
    index += 11
    while index < len(lines):
        stamp = lines[index].split()[0]
        if not re.fullmatch(r"\d{8}\.\d{6}",stamp):
            raise ValueError("Expected a spectral time record")
        times.append((datetime.strptime(stamp,"%Y%m%d.%H%M%S")-ORIGIN).total_seconds()/3600)
        index += 1
        rows = []
        for point in range(points):
            if index >= len(lines) or lines[index].split() != ["LOCATION",str(point+1)]:
                raise ValueError("Missing or reordered spectral location")
            index += 1
            values = np.array([[float(v) for v in row.split()] for row in lines[index:index+n]])
            if values.shape != (n,3) or not np.isfinite(values).all():
                raise ValueError("Incomplete spectral density")
            raw = values[:,0]
            if np.any((raw < 0) & (raw != exception)):
                raise ValueError("Negative spectral variance")
            rows.append(np.where(raw == exception,0,raw))
            index += n
        density.append(rows)
    times = np.asarray(times)
    if len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("Expected ordered spectral times")
    return times, xy, frequency, np.asarray(density)


def spectral_checks(case_path, start_hour=CYCLE_HOURS, end_hour=2*CYCLE_HOURS):
    case = read_case(case_path)
    path = case_path.parent/"diagnostics.spc"
    t, xy, frequency, density = read_spectral_series(path)
    np.testing.assert_array_equal(t,np.arange(0,case["hours"]+1,3))
    np.testing.assert_array_equal(xy,np.loadtxt(case_path.parent/"diagnostics.xy"))
    selected = (t >= start_hour) & (t <= end_hour)
    if not selected.any():
        raise ValueError("Empty spectral comparison interval")
    t, density = t[selected], density[selected]
    m0 = np.trapezoid(density,frequency,axis=-1)
    low = np.divide(np.trapezoid(density[...,:3],frequency[:3],axis=-1),m0,
                    out=np.zeros_like(m0),where=m0>0)
    high = np.divide(np.trapezoid(density[...,-3:],frequency[-3:],axis=-1),m0,
                     out=np.zeros_like(m0),where=m0>0)
    peak = np.argmax(density,axis=-1)
    edge = (peak == 0) | (peak == len(frequency)-1)
    hs = 4*np.sqrt(m0)
    included = hs >= .1
    bad = included & ((low >= .001) | (high >= .01) | edge)
    locations = []
    for p in range(len(xy)):
        samples = included[:,p]
        locations.append(dict(longitude_deg=float(xy[p,0]),latitude_deg=float(xy[p,1]),
            included_spectra=int(samples.sum()),failing_spectra=int(bad[:,p].sum()),
            high_edge_fraction_max=float(high[samples,p].max()) if samples.any() else None,
            low_edge_fraction_max=float(low[samples,p].max()) if samples.any() else None,
            first_failure_hour=float(t[np.flatnonzero(bad[:,p])[0]]) if bad[:,p].any() else None))
    return dict(input_sha256=sha256(path),minimum_resolved_hs_m=.1,
                requested_interval_hours=[start_hour,end_hour],sampled_interval_hours=[float(t[0]),float(t[-1])],
                frequency_hz=frequency.tolist(),locations=locations,
                included_spectra=int(included.sum()),failing_spectra=int(bad.sum()),
                all_edges_pass=bool(included.any() and not bad.any()),
                time_hours=t.tolist(),resolved_hs_m=hs.tolist(),low_edge_fraction=low.tolist(),
                high_edge_fraction=high.tolist(),peak_at_edge=edge.tolist(),
                reading_rule="Three-hourly resolved variance at three wet nodes; Hs>=0.1 m enters the criterion. The lowest two intervals must contain <0.1% of variance, the highest two <1%, with the peak inside the band. The printed band excludes SWAN's diagnostic high-frequency tail.")


def matched_history(reference_path, shorter_path, first_archive_hour, last_archive_hour):
    """Compare the same hourly weather dates, retaining each run's provenance."""
    cases = [read_case(p) for p in (reference_path,shorter_path)]
    rows = []
    for case in cases:
        start = case["forcing"]["start_day_from_first_snapshot"]*24
        first, last = first_archive_hour-start, last_archive_hour-start
        if first < 0 or last > case["hours"] or first % 1 or last % 1:
            raise ValueError("Requested comparison lies outside a run's hourly coverage")
        values = np.array(case["values"])[int(first):int(last)+1]
        values[...,0] += start*3600
        rows.append(dict(name=case["name"],hours=last_archive_hour-first_archive_hour,values=values.tolist()))
    if cases[0]["forcing"]["atmosphere_sha256"] != cases[1]["forcing"]["atmosphere_sha256"]:
        raise ValueError("Comparison atmosphere differs")
    # compare uses its second argument as the reference denominator.
    result = compare(rows[1],rows[0],start_hour=0)
    result["archive_hours"]=[first_archive_hour,last_archive_hour]
    result["reference_case"]=cases[0]["name"]
    result["inputs"]={str(p):sha256(p) for p in (reference_path,shorter_path)}
    reference, alternate = [np.asarray(row["values"]) for row in rows]
    times=reference[:,0,0]/3600
    duration=times[-1]-times[0]
    weights=np.cos(np.radians(reference[0,:,2]))
    weights/=weights.sum()
    occupancy=[]
    for level in (.5,1.,2.,3.):
        totals=[np.array([exceedance(times,data[:,point,3],level)["hours"]
                          for point in range(data.shape[1])]) for data in (reference,alternate)]
        fractions=[float(np.dot(weights,total)/duration) for total in totals]
        occupancy.append(dict(hs_m=level,reference_area_time_fraction=fractions[0],
                              compared_area_time_fraction=fractions[1],
                              area_time_fraction_difference=fractions[1]-fractions[0],
                              node_hours_absolute_difference_max=float(np.max(np.abs(totals[1]-totals[0])))))
    result["threshold_occupancy"]=occupancy
    shape=np.asarray(cases[0]["wet"]).shape
    winds=[]
    for case,p in zip(cases,(reference_path,shorter_path)):
        start=case["forcing"]["start_day_from_first_snapshot"]*24
        data=np.loadtxt(p.parent/"wind.dat").reshape(-1,2,*shape)
        first=int((first_archive_hour-start)*4)
        last=int((last_archive_hour-start)*4)
        winds.append(data[first:last+1])
    error=float(np.max(np.abs(winds[0]-winds[1])))
    if error > 1e-8:
        raise ValueError("Matched wave cases have different wind values")
    result["common_forcing_max_absolute_difference_m_s"]=error
    return result


def main():
    baseline=RUNS/"basin_dt150/product.json"
    selection=refinement_window(read_case(baseline)["forcing"])
    start,end=selection["analysis_start_hour"],selection["analysis_end_hour"]
    a,b=RUNS/"lead96_dt150/product.json",RUNS/"lead96_dt75/product.json"
    result=dict(schema="terluna.climate.wave-cycle-checks/1",producer_sha256=sha256(Path(__file__)),
                spectral_edges=spectral_checks(baseline),
                refinement_window=selection,
                timestep=matched_history(b,a,start,end),
                initial_history=matched_history(baseline,a,start,end),
                refined_spectral_edges=spectral_checks(b,start-selection["run_start_hour"],end-selection["run_start_hour"]))
    output=Path(__file__).with_name("results")/"cycle_checks.json"
    output.write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    print(output)


if __name__ == "__main__":
    main()
