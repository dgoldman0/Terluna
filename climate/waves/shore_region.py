"""Propagate second-cycle spectra through a wider eastern Smythii sector."""
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
from climate.waves.coastal_checks import nesting_spectra
from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import array_text, execute, load_basin
from climate.waves.shore import RUNS, REFERENCE, coastal_files, read_shore, stationary_boundary
from climate.waves.weather import forcing

BOUNDS = (90.751953125,.001953125,6.,5.75)


def perimeter(terrain,stride):
    x,y = terrain["longitude_deg"][::stride],terrain["latitude_deg"][::stride]
    indices = ([(0,i) for i in range(1,len(x))]+
               [(j,len(x)-1) for j in range(1,len(y))]+
               [(len(y)-1,i) for i in range(len(x)-2,-1,-1)]+
               [(j,0) for j in range(len(y)-2,-1,-1)])
    ji=np.array(indices)
    active=terrain["wet"] & (terrain["depth_m"]>.05)
    return np.column_stack((x[ji[:,1]],y[ji[:,0]])),active[ji[:,0]*stride,ji[:,1]*stride]


def remap_spectra(source,xy,wet,spacing_deg,*,neighbors=4):
    """Explicit, bounded interpolation from available sea spectra only.

    Inverse-square weights use up to four nearby wet nodes within two source
    spacings. Each wet target must lie within 1.5 spacings of a source node.
    This coastal extrapolation is recorded as a boundary assumption.
    """
    if len(source["spectra"]) != 1 or neighbors not in (1,4):
        raise ValueError("Select one spectrum time and one or four neighbours")
    valid=source["available"][0]
    points=source["locations"][valid]
    scale=np.array([np.cos(np.deg2rad(np.mean(xy[:,1]))),1.])
    distances,indices=cKDTree(points*scale).query(xy[wet]*scale,k=neighbors)
    if neighbors == 1:
        distances,indices=distances[:,None],indices[:,None]
    if np.any(distances[:,0] > 1.5*spacing_deg):
        raise ValueError("A wet boundary lies beyond the bounded coastal spectrum extrapolation")
    weights=np.where(distances<=2*spacing_deg,1/np.maximum(distances,1e-12)**2,0.)
    weights/=weights.sum(axis=1,keepdims=True)
    values=np.zeros((1,len(xy),len(source["frequency"]),len(source["direction"])))
    values[0,wet]=np.einsum('nk,nkfd->nfd',weights,source["spectra"][0,valid][indices])
    mapped=dict(locations=xy,frequency=source["frequency"],direction=source["direction"],
                spectra=values,available=wet[None,:])
    audit=dict(method="Inverse-square interpolation among available wet spectra",neighbors=neighbors,
        source_spacing_deg=spacing_deg,wet_boundary_locations=int(wet.sum()),
        maximum_nearest_distance_deg=float(distances[:,0].max()),
        maximum_supported_distance_deg=float(np.where(weights>0,distances,0).max()),
        maximum_nearest_distance_in_source_cells=float(distances[:,0].max()/spacing_deg),
        reading_rule="The source and receiving grids resolve different shorelines. Wet spectra are extended over the recorded distances; dry locations retain NODATA. Boundary-position and source-grid controls measure sensitivity to this assumption.")
    return mapped,audit


def regional_case(executable,*,hour=852,stride=16,neighbors=4,run_root=RUNS,
                  accuracy=.005,relaxation=.01,name=None):
    run_root=Path(run_root).absolute()
    producer_bytes=Path(__file__).read_bytes()
    builder_path=Path(__file__).with_name("shore.py")
    builder_bytes=builder_path.read_bytes()
    if hour not in (852,1020,1068) or stride not in (8,16):
        raise ValueError("Use a selected second-cycle phase at 948 or 1895 m spacing")
    terrain_path=run_root/"regional_terrain.npz"
    terrain=load_terrain(terrain_path)
    parent_path=run_root/"parent_fields/product.json"
    parent=json.loads(parent_path.read_text())
    field_path=parent_path.parent/("field_b.spc" if hour==1068 else "field_a.spc")
    if sha256(field_path) != parent["run"]["output_sha256"][field_path.name]:
        raise ValueError("Basin directional spectra changed")
    source=nesting_spectra(field_path,{ORIGIN+timedelta(hours=hour)})
    xy,wet=perimeter(terrain,stride)
    boundary,mapping=remap_spectra(source,xy,wet,1.,neighbors=neighbors)
    build=json.loads((executable.parent.parent/"coupled_air.json").read_text())
    if sha256(executable) != build["build"]["executable_sha256"]:
        raise ValueError("Regional executable differs from its build record")
    basin=load_basin()
    winds=forcing(ROOT/"research/runs/waves/cycle/atmosphere.npz",basin,
                  build["air_density_kg_m3"],48,start_hour=hour)
    wind=dict(lon=basin["lon"],lat=basin["lat"],u=winds["u"][0],v=winds["v"][0])
    settings=next(line for line in (REFERENCE/"INPUT").read_text().splitlines() if line.startswith("SET "))
    files=coastal_files(terrain,stationary_boundary(boundary,0),wind,settings,stride=stride,
                        iterations=600,accuracy=accuracy,relaxation=relaxation)
    xx,yy=np.meshgrid(terrain["longitude_deg"][::stride],terrain["latitude_deg"][::stride])
    active=terrain["wet"][::stride,::stride] & (terrain["depth_m"][::stride,::stride]>.05)
    files["field.xy"]=array_text(np.column_stack((xx[active],yy[active])))
    files["INPUT"]=files["INPUT"].replace("COMPUTE\n",
        "POINTS 'field' FILE 'field.xy'\nSPEC 'field' SPEC2D ABS 'field.spc'\nCOMPUTE\n")
    name=name or f"region_h{hour}_s{stride}_relaxed"+("_nearest" if neighbors==1 else "")
    directory=run_root/name
    run=execute(executable,directory,files,["shore.tbl","field.spc"],timeout_s=3600,threads=2)
    data,wet=read_shore(directory/"shore.tbl",terrain,stride)
    matches=re.findall(r"accuracy OK in\s+([\d.]+) % of wet grid points \(\s*([\d.]+) % required\)",
                       (directory/"PRINT").read_text())
    stopping=dict(last_percent=float(matches[-1][0]) if matches else None,
                  passed=bool(matches and float(matches[-1][0])>=99))
    result=dict(schema="terluna.climate.shore-region/1",name=name,hour=hour,stride=stride,
        evidence="Stationary regional propagation of the second-cycle basin spectra across flooded 118 m terrain.",
        reading_rule="The wider sector supplies the small coastal grids. Basin spectra are explicitly interpolated among nearby wet source nodes at the outer boundary. Local winds retain the atmospheric surface-density correction. Regional and boundary refinements remain separate numerical checks.",
        units=dict(height="m",hour="Earth hours from the basin origin",spacing="degrees of longitude and latitude",spectrum="m2/Hz/degree"),
        stationary_relaxation=relaxation,stationary_accuracy=accuracy,stationary_iteration_limit=600,
        parent_path=str(parent_path.relative_to(ROOT)),parent_sha256=sha256(parent_path),
        terrain_path=str(terrain_path.relative_to(ROOT)),terrain_sha256=sha256(terrain_path),
        atmosphere_sha256=winds["record"]["atmosphere_sha256"],boundary_mapping=mapping,
        producer_sha256=hashlib.sha256(producer_bytes).hexdigest(),
        builder_sha256=hashlib.sha256(builder_bytes).hexdigest(),
        run=run,stopping=stopping,wet_nodes=int(wet.sum()),maximum_hs_m=float(data[wet,2].max()))
    (directory/"producer_source.py").write_bytes(producer_bytes)
    (directory/"builder_source.py").write_bytes(builder_bytes)
    (directory/"product.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(name=name,stopping=stopping,elapsed_wall_s=run["elapsed_wall_s"])),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable",type=Path,required=True)
    parser.add_argument("--hour",type=int,default=852)
    parser.add_argument("--stride",type=int,default=16)
    parser.add_argument("--neighbors",type=int,default=4)
    parser.add_argument("--accuracy",type=float,default=.005)
    parser.add_argument("--relaxation",type=float,default=.01)
    parser.add_argument("--name")
    args=parser.parse_args()
    regional_case(args.executable,hour=args.hour,stride=args.stride,neighbors=args.neighbors,
                  accuracy=args.accuracy,relaxation=args.relaxation,name=args.name)


if __name__ == "__main__":
    main()
