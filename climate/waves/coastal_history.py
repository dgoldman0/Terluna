"""Propagate changing basin spectra through regional and coastal terrain."""
from __future__ import annotations

import argparse
from datetime import timedelta
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.spatial import cKDTree

from climate.waves.coastal import load_terrain
from climate.waves.coastal_checks import iter_nesting_spectra, nesting_spectra
from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import array_text, execute, load_basin
from climate.waves.shore import RUNS as SHORE_RUNS, REFERENCE, SITES, coastal_files, stationary_boundary
from climate.waves.shore_history import RUNS, END_HOUR
from climate.waves.shore_region import perimeter, remap_spectra
from climate.waves.weather import forcing


def stamp(hour):
    return (ORIGIN+timedelta(hours=hour)).strftime("%Y%m%d.%H%M%S")


def write_boundary(path, records, xy, wet, spacing):
    """Stream explicitly remapped, time-ordered spectra to SWAN's nest format."""
    count, previous, audit = 0, None, None
    with path.open("x") as output:
        for source in records:
            date = source["dates"][0]
            if date is None or (previous is not None and date <= previous):
                raise ValueError("Boundary spectra require strictly increasing times")
            mapped, current = remap_spectra(source,xy,wet,spacing)
            text = stationary_boundary(mapped,0)
            header, body = text.split("-99\n",1)
            if count == 0:
                output.write(header.replace("SWAN 1\n","SWAN 1\nTIME\n1\n")+"-99\n")
                audit = current
                audit["start_hour"] = (date-ORIGIN).total_seconds()/3600
            else:
                for key in ("maximum_nearest_distance_deg", "maximum_supported_distance_deg",
                            "maximum_nearest_distance_in_source_cells"):
                    audit[key] = max(audit[key],current[key])
            output.write(date.strftime("%Y%m%d.%H%M%S")+"\n"+body)
            previous = date
            count += 1
    if count < 2:
        raise ValueError("An evolving boundary requires at least two spectra")
    audit.update(records=count,end_hour=(previous-ORIGIN).total_seconds()/3600,
                 sha256=sha256(path))
    return audit


def prepare_boundary(level, stride, source_path, *, name=None, start=0, end=END_HOUR):
    source_path=Path(source_path).absolute()
    terrain = load_terrain(SHORE_RUNS/("regional_terrain.npz" if level=="region" else "terrain.npz"))
    source_record = json.loads((source_path.parent/"product.json").read_text())
    expected_schema=("terluna.climate.shore-history-parent/1" if level=="region" else
                     "terluna.climate.coastal-history-case/1")
    if source_record['schema']!=expected_schema or (level=='coast' and source_record['level']!='region'):
        raise ValueError('History boundary requires its completed upstream model')
    if sha256(source_path) != source_record["run"]["output_sha256"][source_path.name]:
        raise ValueError("Boundary source bytes changed")
    xy,wet=perimeter(terrain,stride)
    spacing=(1. if level=="region" else source_record["stride"]/
             load_terrain(ROOT/source_record['terrain_path'])['metadata']['pixels_per_degree'])
    dates={ORIGIN+timedelta(hours=h) for h in range(start,end+1)}
    path=RUNS/(name or f"{level}_s{stride}_{start}_{end}.nest")
    record_path=path.with_suffix(".json")
    identity=dict(source=str(source_path.relative_to(ROOT)),source_sha256=sha256(source_path),
                  source_record_sha256=sha256(source_path.parent/'product.json'),
                  terrain_sha256=sha256(SHORE_RUNS/("regional_terrain.npz" if level=="region" else "terrain.npz")),
                  stride=stride,start=start,end=end)
    if record_path.exists():
        record=json.loads(record_path.read_text())
        if record["identity"] != identity or sha256(path) != record["mapping"]["sha256"]:
            raise ValueError("Prepared boundary identity differs")
        return path,record
    mapping=write_boundary(path,iter_nesting_spectra(source_path,dates),xy,wet,spacing)
    record=dict(schema="terluna.climate.shore-history-boundary/1",identity=identity,mapping=mapping,
                producer_sha256=sha256(Path(__file__)),
                reading_rule="Hourly spectra are interpolated linearly by SWAN. Wet-source interpolation retains the shore-study distance bounds.")
    record_path.write_text(json.dumps(record,indent=2)+"\n")
    return path,record


def history_files(terrain,boundary,*,stride=4,step_s=900,start=0,end=END_HOUR,
                  hotstart=None,field_xy=None,wind_file=REFERENCE/"wind.dat",checkpoints=()):
    if step_s<=0 or 3600%step_s or not 0<=start<end<=END_HOUR:
        raise ValueError("History must fit the archive with timesteps dividing an hour")
    basin=load_basin()
    wind=dict(lon=basin["lon"],lat=basin["lat"],u=np.zeros_like(basin["depth"]),v=np.zeros_like(basin["depth"]))
    settings=next(line for line in (REFERENCE/"INPUT").read_text().splitlines() if line.startswith("SET "))
    files=coastal_files(terrain,"",wind,settings,stride=stride)
    text=files["INPUT"].replace("MODE STATIONARY","MODE NONSTAT")
    text=re.sub(r"NUM ACCUR[^\n]*\n", "NUM ACCUR .005 .005 .005 99 NONSTAT 1\n",text)
    text=re.sub(r"INPGRID WIND([^\n]*)\nREADINP WIND[^\n]*\n",
                rf"INPGRID WIND\1 &\n NONSTAT {stamp(0)} 15 MIN {stamp(END_HOUR)}\n"
                "READINP WIND 1 'wind.dat' 3 0 0 0 FREE\n",text)
    if hotstart:
        text=text.replace("INIT ZERO", "INIT HOTSTART SINGLE 'initial.hot' UNFORMATTED")
        files["initial.hot"]=hotstart
    files["boundary.nest"]=boundary
    files["wind.dat"]=wind_file
    text=text[:text.index("TABLE 'COMPGRID'")]
    text+=f"QUANTITY TSEC REF {stamp(0)}\n"
    text+=("TABLE 'COMPGRID' NOHEAD 'maps.tbl' TSEC XP YP HS RTP TM01 DEP DIR QB TRANSP DISSURF &\n"
           f" OUTPUT {stamp(start)} 6 HR\n")
    text+=("TABLE 'COMPGRID' NOHEAD 'final.tbl' TSEC XP YP HS RTP TM01 DEP DIR QB TRANSP DISSURF &\n"
           f" OUTPUT {stamp(end)} 6 HR\n")
    files["sites.xy"]=array_text(np.array(list(SITES.values())))
    text+=("POINTS 'sites' FILE 'sites.xy'\n"
           "TABLE 'sites' NOHEAD 'sites.tbl' TSEC XP YP HS RTP TM01 DEP DIR QB TRANSP DISSURF &\n"
           f" OUTPUT {stamp(start)} 15 MIN\n"
           f"SPEC 'sites' SPEC2D ABS 'sites.spc' OUTPUT {stamp(start)} 1 HR\n")
    if field_xy is not None:
        files["field.xy"]=array_text(field_xy)
        text+=("POINTS 'field' FILE 'field.xy'\n"
               f"SPEC 'field' SPEC2D ABS 'field.spc' OUTPUT {stamp(start)} 1 HR\n")
    previous=start
    for hour in sorted(set(h for h in (*checkpoints,end) if start<h<=end)):
        text+=f"COMPUTE NONSTAT {stamp(previous)} {step_s} SEC {stamp(hour)}\n"
        text+=f"HOTFILE 'h{hour}.hot' UNFORMATTED\n"
        previous=hour
    files["INPUT"]=text+"STOP\n"
    return files


def history_case(executable,level,*,stride=None,step_s=900,start=0,end=END_HOUR,
                 source=None,initial=None,name=None,timeout=10800):
    producer=Path(__file__).read_bytes()
    if level not in ("region","coast"):
        raise ValueError("Choose regional or coastal propagation")
    if start>0 and initial is None:
        raise ValueError('A later history window requires its recorded checkpoint')
    reference=json.loads((REFERENCE/'run.json').read_text())
    for filename in ('INPUT','wind.dat'):
        if sha256(REFERENCE/filename)!=reference['identity']['inputs'][filename]:
            raise ValueError('The verified basin settings or wind history changed')
    stride=stride or (16 if level=="region" else 4)
    name=name or f"{level}_s{stride}_dt{step_s}_{start}_{end}"
    source=source or (RUNS/"parent/field.spc" if level=="region" else RUNS/"region_s16_dt900_0_1419/field.spc")
    boundary,mapping=prepare_boundary(level,stride,source,start=start,end=end)
    terrain_path=SHORE_RUNS/("regional_terrain.npz" if level=="region" else "terrain.npz")
    terrain=load_terrain(terrain_path)
    field_xy=None
    if level=="region":
        local=load_terrain(SHORE_RUNS/"terrain.npz")
        xy,wet=perimeter(local,4)
        xx,yy=np.meshgrid(terrain["longitude_deg"][::stride],terrain["latitude_deg"][::stride])
        active=terrain["wet"][::stride,::stride] & (terrain["depth_m"][::stride,::stride]>.05)
        points=np.column_stack((xx[active],yy[active]))
        scale=np.array([np.cos(np.deg2rad(np.mean(xy[:,1]))),1.])
        _,idx=cKDTree(points*scale).query(xy[wet]*scale,k=4)
        field_xy=points[np.unique(idx)]
    checkpoints=tuple(h for h in (696,816,1008,1068) if start<h<end)
    files=history_files(terrain,boundary,stride=stride,step_s=step_s,start=start,end=end,
                        hotstart=initial,field_xy=field_xy,checkpoints=checkpoints)
    outputs=["maps.tbl","final.tbl","sites.tbl","sites.spc"]+[f"h{h}.hot" for h in (*checkpoints,end)]
    if field_xy is not None:
        outputs.append("field.spc")
    build=json.loads((executable.parent.parent/"coupled_air.json").read_text())
    if sha256(executable)!=build["build"]["executable_sha256"]:
        raise ValueError("History executable differs from its coupled-air build")
    directory=RUNS/name
    run=execute(executable,directory,files,outputs,timeout_s=timeout,threads=8)
    product_path=directory/'product.json'
    if product_path.exists():
        previous=json.loads(product_path.read_text())
        if (previous['schema']!='terluna.climate.coastal-history-case/1' or
                previous['run']!=run or previous['boundary']!=mapping or
                previous['terrain_sha256']!=sha256(terrain_path) or
                sha256(directory/'producer_source.py')!=previous['producer_sha256']):
            raise ValueError('Completed history metadata differs from its inputs or source record')
        print(json.dumps(dict(name=name,cached=True)),flush=True)
        return previous
    data=np.loadtxt(directory/"sites.tbl",ndmin=2)
    expected=(end-start)*4+1
    # Adjacent COMPUTE commands can repeat their shared output time.
    data=data.reshape(-1,len(SITES),12)
    _,unique=np.unique(data[:,0,0],return_index=True)
    data=data[np.sort(unique)]
    np.testing.assert_allclose(data[:,0,0],np.arange(start,end+.01,.25)*3600,atol=.01,rtol=0)
    if data.shape!=(expected,len(SITES),12) or not np.isfinite(data).all():
        raise ValueError("Invalid evolving station output")
    if level=="coast" and np.any(data[:,:,3]<0):
        raise ValueError("A coastal station became unavailable")
    result=dict(schema="terluna.climate.coastal-history-case/1",name=name,level=level,
        start_hour=start,end_hour=end,stride=stride,step_s=step_s,
        evidence="Nonstationary spectral propagation with changing basin boundary spectra and GCM surface stress.",
        reading_rule="All levels retain the first solar cycle before second-cycle statistics. Boundary spectra are hourly; station wave diagnostics are every 15 minutes. Spatial, temporal and physical limits accompany the analysis.",
        terrain_path=str(terrain_path.relative_to(ROOT)),terrain_sha256=sha256(terrain_path),
        boundary=mapping,initial_sha256=sha256(initial) if initial else None,
        run=run,producer_sha256=hashlib.sha256(producer).hexdigest(),
        units=dict(time="Earth hours from basin origin",height="m",period="s",variance_transport="m3/s"),
        sites=SITES,station_records=int(data.shape[0]),
        columns=["time_s","longitude_deg","latitude_deg","hs_m","peak_period_s","mean_period_s",
                 "depth_m","direction_deg","breaking_fraction","transport_east_m3_s","transport_north_m3_s","breaking_loss_m2_s"])
    (directory/"producer_source.py").write_bytes(producer)
    product_path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(name=name,elapsed_wall_s=run["elapsed_wall_s"])),flush=True)
    return result


def preview(executable,threads=2):
    """Bound run cost using the known eastward state before a full history."""
    terrain=load_terrain(SHORE_RUNS/"terrain.npz")
    source=nesting_spectra(SHORE_RUNS/"region_h852_s16_relaxed/field.spc")
    xy,wet=perimeter(terrain,4)
    boundary,_=remap_spectra(source,xy,wet,16/256)
    basin=load_basin()
    build=json.loads((executable.parent.parent/"coupled_air.json").read_text())
    winds=forcing(ROOT/"research/runs/waves/cycle/atmosphere.npz",basin,build["air_density_kg_m3"],48,start_hour=852)
    wind=dict(lon=basin["lon"],lat=basin["lat"],u=winds["u"][0],v=winds["v"][0])
    settings=next(l for l in (REFERENCE/"INPUT").read_text().splitlines() if l.startswith("SET "))
    files=coastal_files(terrain,stationary_boundary(boundary,0),wind,settings,stride=4,sites=SITES)
    text=files["INPUT"].replace("MODE STATIONARY","MODE NONSTAT")
    text=re.sub(r"NUM ACCUR[^\n]*", "NUM ACCUR .005 .005 .005 99 NONSTAT 1",text)
    text=text.replace("COMPUTE\n",f"COMPUTE NONSTAT {stamp(0)} 900 SEC {stamp(48)}\n")
    files["INPUT"]=text
    run=execute(executable,RUNS/("preview" if threads==2 else f"preview_{threads}threads"),files,
                ["shore.tbl","sites.spc"],timeout_s=900,threads=threads)
    print(json.dumps(dict(preview_elapsed_s=run["elapsed_wall_s"],projected_cycle_s=run["elapsed_wall_s"]*END_HOUR/48)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("level",choices=("region","coast","preview"))
    parser.add_argument("--executable",type=Path,required=True)
    parser.add_argument("--stride",type=int)
    parser.add_argument("--step",type=int,default=900)
    parser.add_argument("--start",type=int,default=0)
    parser.add_argument("--end",type=int,default=END_HOUR)
    parser.add_argument("--source",type=Path)
    parser.add_argument("--initial",type=Path)
    parser.add_argument("--name")
    args=parser.parse_args()
    if args.level=="preview":
        preview(args.executable)
    else:
        history_case(args.executable,args.level,stride=args.stride,step_s=args.step,start=args.start,
                     end=args.end,source=args.source,initial=args.initial,name=args.name)


if __name__=="__main__":
    main()
