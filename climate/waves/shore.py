"""Resolve coastal wave exposure with spectra from the second solar cycle."""
from __future__ import annotations

import argparse
from datetime import timedelta
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import array_text, execute, load_basin
from climate.waves.coastal import load_terrain
from climate.waves.coastal_checks import nesting_spectra
from climate.waves.cycle import CYCLE_HOURS
from climate.waves.weather import forcing
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
RUNS = ROOT / "research/runs/waves/shore"
REFERENCE = ROOT / "research/runs/waves/cycle/basin_dt150"
BOUNDS = (92.751953125, 1.751953125, 1., 1.)
EXPANDED_BOUNDS = (92.501953125, 1.501953125, 1.5, 1.5)
PARENT_HOURS = 1080
EPISODE_HOURS = (852, 1020, 1068)
COLUMNS = ["longitude_deg", "latitude_deg", "hs_m", "peak_period_s", "mean_period_s",
           "depth_m", "direction_deg", "breaking_fraction", "variance_transport_east_m3_s",
           "variance_transport_north_m3_s", "breaking_variance_loss_m2_s"]
SITES = {"west_face": (93.251953125,2.361328125),
         "east_face": (93.408203125,2.330078125),
         "southern_shore": (92.986328125,1.830078125),
         "offshore": (92.845703125,2.345703125)}


def stationary_boundary(boundary, index):
    """Preserve a selected frequency-direction distribution at each nest location."""
    xy, frequency, direction = [boundary[k] for k in ("locations", "frequency", "direction")]
    text = "SWAN 1\nLONLAT\n"+str(len(xy))+"\n"+array_text(xy)
    text += "RFREQ\n"+str(len(frequency))+"\n"+array_text(frequency)
    text += "CDIR\n"+str(len(direction))+"\n"+array_text(direction)
    text += "QUANT\n1\nVaDens\nm2/Hz/degr\n-99\n"
    for available, spectrum in zip(boundary["available"][index], boundary["spectra"][index]):
        if not available:
            text += "NODATA\n"
        elif not spectrum.any():
            text += "ZERO\n"
        else:
            factor = float(spectrum.max()/99999999)
            integers = np.rint(spectrum/factor).astype(np.int64)
            text += f"FACTOR\n{factor:.14e}\n"
            text += "\n".join(" ".join(map(str,row)) for row in integers)+"\n"
    return text


def boundary_audit(terrain, boundary, index):
    x, y = terrain["longitude_deg"], terrain["latitude_deg"]
    ppd = terrain["metadata"]["pixels_per_degree"]
    ij = np.rint((boundary["locations"]-[x[0],y[0]])*ppd).astype(int)
    if np.any(ij < 0) or np.any(ij[:,0] >= len(x)) or np.any(ij[:,1] >= len(y)):
        raise ValueError("Nesting locations lie outside the coastal grid")
    wet = terrain["wet"][ij[:,1],ij[:,0]]
    missing = wet & ~boundary["available"][index]
    if missing.any():
        raise ValueError(f"Wet coastal boundary lacks parent spectra: {boundary['locations'][missing].tolist()}")
    frequency = boundary["frequency"]
    density = boundary["spectra"][index].sum(axis=-1)*(360/len(boundary["direction"]))
    variance = np.trapezoid(density,frequency,axis=-1)
    selected = wet & (variance >= (.1/4)**2)
    if not selected.any():
        raise ValueError("Selected boundary has negligible incident waves")
    low = np.trapezoid(density[:,:3],frequency[:3],axis=-1)[selected]/variance[selected]
    high = np.trapezoid(density[:,-3:],frequency[-3:],axis=-1)[selected]/variance[selected]
    peak = np.argmax(density[selected],axis=-1)
    return dict(wet_locations=int(wet.sum()),missing_parent_spectra=int(missing.sum()),
                frequency_bins=len(frequency),direction_bins=len(boundary["direction"]),
                maximum_low_edge_fraction=float(low.max()),maximum_high_edge_fraction=float(high.max()),
                spectral_edges_pass=bool(low.max()<.001 and high.max()<.01 and np.all((peak>0)&(peak<len(frequency)-1))),
                resolved_hs_range_m=(4*np.sqrt([variance[selected].min(),variance[selected].max()])).tolist())


def coastal_files(terrain, boundary_text, wind, settings, *, stride=4, directions=36,
                  iterations=300, accuracy=.005, extra_frequencies=0, sites=None, relaxation=.01):
    """A stationary response to one second-cycle weather/spectrum snapshot."""
    if stride not in (1,2,4,8,16) or directions not in (36,72) or extra_frequencies not in (0,7):
        raise ValueError("Unsupported coastal refinement")
    xi, yi = terrain["longitude_deg"], terrain["latitude_deg"]
    x, y = xi[::stride], yi[::stride]
    if x[-1] != xi[-1] or y[-1] != yi[-1]:
        raise ValueError("Every coastal grid must reach the same boundary")
    wx, wy = wind["lon"], wind["lat"]
    if wind["u"].shape != (len(wy),len(wx)) or wind["v"].shape != (len(wy),len(wx)):
        raise ValueError("Wind must retain the parent input grid")
    ratio = MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    intervals = 48+extra_frequencies
    high_factor = 100**(extra_frequencies/48)
    # Retain the reference water-density scenario and variance-output convention.
    if not settings.startswith("SET ") or "\n" in settings:
        raise ValueError("Supply the parent SET command")
    text = f"""PROJECT 'Smythii shore' '006'
MODE STATIONARY TWODIMENSIONAL
COORDINATES SPHERICAL {MOON_RADIUS:.12g} CCM
{settings}
OUTPUT OPTIONS TABLE 16
CGRID REG {x[0]:.12g} {y[0]:.12g} 0 {x[-1]-x[0]:.12g} {y[-1]-y[0]:.12g} {len(x)-1} {len(y)-1} &
 CIRCLE {directions} {0.03*ratio:.12g} {3*ratio*high_factor:.12g} {intervals}
INPGRID BOTTOM {xi[0]:.12g} {yi[0]:.12g} 0 {len(xi)-1} {len(yi)-1} {xi[1]-xi[0]:.12g} {yi[1]-yi[0]:.12g}
READINP BOTTOM 1 'depth.bot' 3 0 FREE
BOUNDNEST1 NEST 'boundary.nest' CLOSED
INPGRID WIND {wx[0]:.12g} {wy[0]:.12g} 0 {len(wx)-1} {len(wy)-1} {wx[1]-wx[0]:.12g} {wy[1]-wy[0]:.12g}
READINP WIND 1 'wind.dat' 3 0 0 FREE
GEN3 KOMEN DRAG WU AGROW
BREAKING CONSTANT 1.0 0.73
PROP BSBT
INIT ZERO
NUM ACCUR {accuracy} {accuracy} {accuracy} 99 STAT {iterations} {relaxation}
TABLE 'COMPGRID' NOHEAD 'shore.tbl' XP YP HS RTP TM01 DEP DIR QB TRANSP DISSURF
COMPUTE
STOP
"""
    files = dict(INPUT=text, **{"depth.bot":array_text(terrain["depth_m"]),
                               "boundary.nest":boundary_text,
                               "wind.dat":array_text(wind["u"])+array_text(wind["v"])})
    if sites:
        files["sites.xy"] = array_text(np.array(list(sites.values())))
        output = "POINTS 'sites' FILE 'sites.xy'\nSPEC 'sites' SPEC2D ABS 'sites.spc'\n"
        files["INPUT"] = text.replace("COMPUTE\n",output+"COMPUTE\n")
    return files


def read_shore(path, terrain, stride):
    x, y = terrain["longitude_deg"][::stride], terrain["latitude_deg"][::stride]
    depth = terrain["depth_m"][::stride,::stride]
    wet = terrain["wet"][::stride,::stride] & (depth > .05)
    data = np.loadtxt(path,ndmin=2)
    if data.shape != (len(x)*len(y),len(COLUMNS)) or not np.isfinite(data).all():
        raise ValueError("Incomplete coastal energy output")
    data = data.reshape(len(y),len(x),len(COLUMNS))
    xx, yy = np.meshgrid(x,y)
    np.testing.assert_allclose(data[wet,:2],np.column_stack((xx[wet],yy[wet])),atol=1e-5,rtol=1e-6)
    np.testing.assert_allclose(data[wet,5],depth[wet],atol=.01,rtol=6e-5)
    if np.any(data[wet,2] < 0) or np.any(data[wet,7] < 0) or np.any(data[wet,7] > 1):
        raise ValueError("Invalid coastal wave height or breaking fraction")
    return data, wet


def snapshot_case(executable, *, hour=852, stride=4, directions=36, expanded=False,
                  iterations=300, accuracy=.005, extra_frequencies=0, run_root=RUNS, name=None,
                  regional_stride=16, regional_neighbors=4, regional_name=None):
    run_root = Path(run_root).absolute()
    producer_bytes = Path(__file__).read_bytes()
    if not CYCLE_HOURS <= hour <= 2*CYCLE_HOURS or hour % 3 or hour > PARENT_HOURS:
        raise ValueError("Select a saved snapshot within the second solar cycle")
    terrain_path = run_root/("expanded_terrain.npz" if expanded else "terrain.npz")
    terrain = load_terrain(terrain_path)
    from climate.waves.shore_region import perimeter, remap_spectra
    regional_name = regional_name or f"region_h{hour}_s{regional_stride}_relaxed"+("_nearest" if regional_neighbors==1 else "")
    regional_path = run_root/regional_name/"product.json"
    regional = json.loads(regional_path.read_text())
    if regional["schema"] != "terluna.climate.shore-region/1" or regional["hour"] != hour:
        raise ValueError("The regional source must describe this weather phase")
    if not regional["stopping"]["passed"]:
        raise ValueError("The regional source must satisfy its stopping criterion")
    nest_path = regional_path.parent/"field.spc"
    if sha256(nest_path) != regional["run"]["output_sha256"][nest_path.name]:
        raise ValueError("Regional boundary bytes changed")
    source = nesting_spectra(nest_path)
    xy,wet = perimeter(terrain,stride)
    boundary,mapping = remap_spectra(source,xy,wet,regional["stride"]/terrain["metadata"]["pixels_per_degree"])
    audit = boundary_audit(terrain,boundary,0)
    build = json.loads((executable.parent.parent/"coupled_air.json").read_text())
    if sha256(executable) != build["build"]["executable_sha256"]:
        raise ValueError("Coastal executable differs from its build record")
    basin = load_basin()
    winds = forcing(ROOT/"research/runs/waves/cycle/atmosphere.npz",basin,
                    build["air_density_kg_m3"],48,start_hour=hour)
    settings = next(line for line in (REFERENCE/"INPUT").read_text().splitlines() if line.startswith("SET "))
    wind = dict(lon=basin["lon"],lat=basin["lat"],u=winds["u"][0],v=winds["v"][0])
    files = coastal_files(terrain,stationary_boundary(boundary,0),wind,settings,
                          stride=stride,directions=directions,iterations=iterations,
                          accuracy=accuracy,extra_frequencies=extra_frequencies,sites=SITES)
    name = name or f"h{hour}_s{stride}_d{directions}"+("_expanded" if expanded else "")
    directory = run_root/name
    run = execute(executable,directory,files,["shore.tbl","sites.spc"],timeout_s=3600,threads=2)
    data, wet = read_shore(directory/"shore.tbl",terrain,stride)
    matches = re.findall(r"accuracy OK in\s+([\d.]+) % of wet grid points \(\s*([\d.]+) % required\)",
                         (directory/"PRINT").read_text())
    stopping = dict(criterion_percent=99.,last_percent=float(matches[-1][0]) if matches else None,
                    passed=bool(matches and float(matches[-1][0]) >= float(matches[-1][1])))
    rho = float(re.search(r"\bRHO\s+([\d.eE+-]+)",files["INPUT"])[1])
    result = dict(schema="terluna.climate.shore-case/1",name=name,
        evidence="Stationary coastal response to a directional spectrum and wind field sampled from the second lunar solar cycle.",
        reading_rule="Each case holds its sampled boundary and winds fixed while the coastal solution settles. Spatial comparisons use identical native terrain. This product measures exposure at selected weather phases; event durations require coastal time evolution.",
        hour=hour,day_in_second_cycle=(hour-CYCLE_HOURS)/24,stride=stride,directions=directions,
        frequency_intervals=48+extra_frequencies,stationary_accuracy=accuracy,sites=SITES,
        stationary_relaxation=.01,
        spacing_m=MOON_RADIUS*np.pi/180*stride/terrain["metadata"]["pixels_per_degree"],
        terrain_path=str(terrain_path.relative_to(ROOT)),terrain_sha256=sha256(terrain_path),
        parent_nest_sha256=sha256(nest_path),atmosphere_sha256=winds["record"]["atmosphere_sha256"],
        boundary_source=dict(path=str(regional_path.relative_to(ROOT)),sha256=sha256(regional_path)),
        boundary_mapping=mapping,mapping_producer_sha256=sha256(Path(__file__).with_name("shore_region.py")),
        producer_sha256=hashlib.sha256(producer_bytes).hexdigest(),run=run,boundary=audit,stopping=stopping,
        water_density_kg_m3=rho,gravity_m_s2=MOON_SURFACE_GRAVITY,
        energy_flux_conversion=rho*MOON_SURFACE_GRAVITY,
        columns=COLUMNS,longitude_deg=terrain["longitude_deg"][::stride].tolist(),
        latitude_deg=terrain["latitude_deg"][::stride].tolist(),wet=wet.tolist(),
        values=data.tolist(),maximum_hs_m=float(data[wet,2].max()))
    (directory/"producer_source.py").write_bytes(producer_bytes)
    (directory/"product.json").write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    print(json.dumps(dict(name=name,max_hs_m=result["maximum_hs_m"],stopping=stopping,
                         elapsed_wall_s=run["elapsed_wall_s"])),flush=True)
    return result


def parent_files(reference=REFERENCE, *, fields=False):
    """Replay the verified continuous sea, adding directional coastal outputs."""
    record = json.loads((reference / "run.json").read_text())
    files = {}
    for name, expected in record["identity"]["inputs"].items():
        path = reference / name
        if sha256(path) != expected:
            raise ValueError(f"Parent reference input changed: {name}")
        files[name] = path.read_text()
    start = ORIGIN.strftime("%Y%m%d.%H%M%S")
    old_end = (ORIGIN + timedelta(hours=1419)).strftime("%Y%m%d.%H%M%S")
    end = (ORIGIN + timedelta(hours=PARENT_HOURS)).strftime("%Y%m%d.%H%M%S")
    text = files["INPUT"].replace(old_end, end)
    extra = ""
    for name, bounds in (("shore", BOUNDS), ("expanded", EXPANDED_BOUNDS)):
        west, south, width, height = bounds
        extra += (f"NGRID '{name}' {west:.12g} {south:.12g} 0 {width} {height} "
                  f"{round(width*16)} {round(height*16)}\n"
                  f"NESTOUT '{name}' '{name}.nest' OUTPUT {start} 3 HR\n")
    extra += ("TABLE 'COMPGRID' NOHEAD 'transport.tbl' TSEC XP YP TRANSP &\n"
              f" OUTPUT {start} 3 HR\n")
    if fields:
        files["field.xy"] = array_text(load_basin()["points"])
        extra += "POINTS 'field' FILE 'field.xy'\n"
        for name,hour,interval in (("field_a",852,168),("field_b",1068,2000)):
            date = (ORIGIN+timedelta(hours=hour)).strftime("%Y%m%d.%H%M%S")
            extra += f"SPEC 'field' SPEC2D ABS '{name}.spc' OUTPUT {date} {interval} HR\n"
    files["INPUT"] = text.replace("COMPUTE NONSTAT", extra + "COMPUTE NONSTAT")
    return files


def replay_parent(executable, run_root=RUNS, *, fields=False):
    run_root = Path(run_root).absolute()
    record = json.loads((REFERENCE/"run.json").read_text())
    if sha256(executable) != record["identity"]["executable_sha256"]:
        raise ValueError("Replay requires the recorded continuous-cycle executable")
    directory = run_root/("parent_fields" if fields else "parent")
    run = execute(executable, directory, parent_files(fields=fields),
                  ["basin.tbl", "wind.tbl", "offshore.spc", "diagnostics.spc",
                   "shore.nest", "expanded.nest", "transport.tbl"]+
                  (["field_a.spc","field_b.spc"] if fields else []),
                  timeout_s=7200, threads=4)
    reference = json.loads((REFERENCE/"product.json").read_text())
    from climate.waves.pilot import BasinCase, load_basin, read_basin
    case = BasinCase("shore_parent", stride=4, step_s=150,
                     frequency_intervals=48, hours=PARENT_HOURS)
    actual = read_basin(directory/"basin.tbl", case, load_basin(stride=4))
    expected = np.asarray(reference["values"][:PARENT_HOURS+1])
    np.testing.assert_array_equal(actual, expected)
    result = dict(schema="terluna.climate.shore-parent/1", run=run,
                  evidence="Continuous replay of the established wave history with added coastal spectra and energy-transport diagnostics.",
                  reference_case_sha256=sha256(REFERENCE/"product.json"),
                  reference_values_compared=int(actual.size),
                  reference_values_exact=True, hours=PARENT_HOURS,
                  bounds=dict(shore=BOUNDS, expanded=EXPANDED_BOUNDS),
                  transport_units="m3/s of wave variance transport; multiply by the configured water density and lunar gravity for W/m")
    (directory/"product.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(hours=PARENT_HOURS, values_exact=int(actual.size),
                         elapsed_wall_s=run["elapsed_wall_s"])), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("parent","fields","snapshot"))
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, default=RUNS)
    parser.add_argument("--hour",type=int,default=852)
    parser.add_argument("--stride",type=int,default=4)
    parser.add_argument("--directions",type=int,default=36)
    parser.add_argument("--expanded",action="store_true")
    parser.add_argument("--iterations",type=int,default=300)
    parser.add_argument("--accuracy",type=float,default=.005)
    parser.add_argument("--extra-frequencies",type=int,choices=(0,7),default=0)
    parser.add_argument("--name")
    parser.add_argument("--regional-stride",type=int,choices=(8,16),default=16)
    parser.add_argument("--regional-neighbors",type=int,choices=(1,4),default=4)
    parser.add_argument("--regional-name")
    args = parser.parse_args()
    if args.action in ("parent","fields"):
        replay_parent(args.executable,args.run_root,fields=args.action == "fields")
    else:
        snapshot_case(args.executable,hour=args.hour,stride=args.stride,directions=args.directions,
                      expanded=args.expanded,iterations=args.iterations,accuracy=args.accuracy,
                      extra_frequencies=args.extra_frequencies,run_root=args.run_root,name=args.name,
                      regional_stride=args.regional_stride,regional_neighbors=args.regional_neighbors,
                      regional_name=args.regional_name)


if __name__ == "__main__":
    main()
