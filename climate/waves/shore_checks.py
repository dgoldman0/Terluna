"""Check coastal exposure, directional energy transport and spatial refinement."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from climate.waves.coastal import load_terrain
from climate.waves.coastal_checks import nesting_spectra
from climate.waves.model import ROOT, sha256
from climate.waves.shore import RUNS, SITES, REFERENCE, EPISODE_HOURS
from climate.waves.cycle import CYCLE_HOURS
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY

OUTPUT = Path(__file__).with_name("results")/"shore.json"


def episode_context():
    path=REFERENCE/"product.json"
    values=np.asarray(json.loads(path.read_text())["values"])
    xy=values[0,:,1:3]
    region=(xy[:,0]>=90)&(xy[:,0]<=94)&(xy[:,1]>=0)&(xy[:,1]<=4)
    hours=np.arange(len(values))
    selected=(hours>=CYCLE_HOURS)&(hours<=2*CYCLE_HOURS)
    local=values[:,region]
    proxy=(local[...,3]**2*local[...,5]*np.maximum(np.cos(np.deg2rad(local[...,7])),0)).mean(axis=1)
    peak=int(hours[selected][np.argmax(proxy[selected])])
    phases=[]
    for hour in EPISODE_HOURS:
        row=local[hour]
        angle=np.deg2rad(row[:,7])
        heading=np.arctan2((row[:,3]**2*np.sin(angle)).sum(),(row[:,3]**2*np.cos(angle)).sum())
        phases.append(dict(hour=hour,regional_mean_hs_m=float(row[:,3].mean()),
            regional_mean_tm01_s=float(row[:,5].mean()),
            variance_weighted_travel_direction_deg=float(np.rad2deg(heading)%360)))
    return dict(source_path=str(path.relative_to(ROOT)),source_sha256=sha256(path),
        region_deg=[90,94,0,4],wet_nodes=int(region.sum()),
        selection_proxy="Mean Hs^2 * Tm01 * max(cos(travel_direction), 0) over the region",
        second_cycle_peak_proxy_hour=peak,nearest_saved_spectrum_hour=3*round(peak/3),phases=phases)


def group_velocity(frequency, depth, gravity=MOON_SURFACE_GRAVITY):
    """Linear-wave group speed from the finite-depth dispersion relation."""
    frequency, depth = np.asarray(frequency), np.asarray(depth)
    if np.any(frequency <= 0) or np.any(depth <= 0):
        raise ValueError("Positive frequencies and water depths are required")
    omega = 2*np.pi*frequency
    target = omega**2*depth[..., None]/gravity
    lo, hi = np.zeros_like(target), np.maximum(target, np.sqrt(target))+1
    for _ in range(60):
        middle = (lo+hi)/2
        below = middle*np.tanh(middle) < target
        lo, hi = np.where(below,middle,lo), np.where(below,hi,middle)
    kh = (lo+hi)/2
    ratio = np.zeros_like(kh)
    moderate = kh < 350
    ratio[moderate] = 2*kh[moderate]/np.sinh(2*kh[moderate])
    return .5*omega*depth[..., None]/kh*(1+ratio)


def spectral_transport(frequency, direction, density, depth, rho, gravity, normal=None):
    """Integrate variance spectra in m2/Hz/degree into power per metre."""
    frequency, direction, density = map(np.asarray,(frequency,direction,density))
    if density.shape != (len(frequency),len(direction)) or np.any(density < 0):
        raise ValueError("Expected one nonnegative frequency-direction spectrum")
    if len(direction) < 3 or not np.allclose(np.diff(direction),360/len(direction)):
        raise ValueError("Expected a complete uniformly spaced directional circle")
    unit = np.column_stack((np.cos(np.deg2rad(direction)),np.sin(np.deg2rad(direction))))
    angular_power = np.trapezoid(group_velocity(frequency,depth,gravity)[:,None]*density,
                                 frequency,axis=0)*(360/len(direction))*rho*gravity
    vector = angular_power@unit
    result = dict(vector_w_m=vector.tolist(),magnitude_w_m=float(np.linalg.norm(vector)),
                  all_directions_w_m=float(angular_power.sum()))
    if normal is not None:
        normal = np.asarray(normal,dtype=float)
        if normal.shape != (2,) or not np.isclose(np.linalg.norm(normal),1):
            raise ValueError("Shoreward normal must have unit length")
        projection = unit@normal
        incoming = float(angular_power@np.maximum(projection,0))
        outgoing = float(angular_power@np.maximum(-projection,0))
        result.update(shoreward_w_m=incoming,seaward_w_m=outgoing,
                      net_shoreward_w_m=incoming-outgoing)
    return result


def verify_case(path):
    product = json.loads(path.read_text())
    if product["schema"] != "terluna.climate.shore-case/1":
        raise ValueError("Unsupported shore case")
    run = product["run"]
    for filename,digest in {**run["identity"]["inputs"],**run["output_sha256"]}.items():
        if sha256(path.parent/filename) != digest:
            raise ValueError(f"Changed coastal run file: {path.parent/filename}")
    if sha256(path.parent/"producer_source.py") != product["producer_sha256"]:
        raise ValueError("Recorded coastal producer changed")
    if sha256(ROOT/product["terrain_path"]) != product["terrain_sha256"]:
        raise ValueError("Coastal terrain changed")
    regional=product["boundary_source"]
    if sha256(ROOT/regional["path"]) != regional["sha256"]:
        raise ValueError("Regional spectrum producer record changed")
    data = np.asarray(product["values"])
    np.testing.assert_array_equal(np.loadtxt(path.parent/"shore.tbl").reshape(data.shape),data)
    return product


def shore_normal(terrain,xy):
    """Local direction toward shallower water on the native rock surface."""
    x,y = terrain["longitude_deg"],terrain["latitude_deg"]
    i,j = np.argmin(abs(x-xy[0])),np.argmin(abs(y-xy[1]))
    np.testing.assert_allclose([x[i],y[j]],xy,rtol=0,atol=1e-10)
    physical = terrain["metadata"]["sea_level_m"]-terrain["height_m"]
    dx = MOON_RADIUS*np.deg2rad(x[i+1]-x[i-1])*np.cos(np.deg2rad(y[j]))
    dy = MOON_RADIUS*np.deg2rad(y[j+1]-y[j-1])
    gradient = np.array([(physical[j,i+1]-physical[j,i-1])/dx,
                         (physical[j+1,i]-physical[j-1,i])/dy])
    if np.linalg.norm(gradient) < 1e-8:
        raise ValueError("Shoreward direction is undefined on a level bottom")
    return -gradient/np.linalg.norm(gradient)


def case_summary(case,path):
    values,wet = np.asarray(case["values"]),np.asarray(case["wet"])
    terrain = load_terrain(ROOT/case["terrain_path"])
    spectra = nesting_spectra(path.parent/"sites.spc")
    x,y = np.asarray(case["longitude_deg"]),np.asarray(case["latitude_deg"])
    factor = case["energy_flux_conversion"]
    sites = {}
    for k,(name,xy) in enumerate(SITES.items()):
        np.testing.assert_allclose(spectra["locations"][k],xy,rtol=0,atol=1e-5)
        i,j = np.argmin(abs(x-xy[0])),np.argmin(abs(y-xy[1]))
        # Shared physical stations sit exactly on the 474, 237 and 118 m grids.
        if not np.allclose([x[i],y[j]],xy,rtol=0,atol=1e-10):
            continue
        if not wet[j,i] or not spectra["available"][0,k]:
            raise ValueError(f"Dry comparison station: {name}")
        row = values[j,i]
        normal = None if name == "offshore" else shore_normal(terrain,xy)
        power = spectral_transport(spectra["frequency"],spectra["direction"],
                                   spectra["spectra"][0,k],row[5],case["water_density_kg_m3"],
                                   case["gravity_m_s2"],normal)
        table_vector = row[8:10]*factor
        error = np.linalg.norm(np.asarray(power["vector_w_m"])-table_vector)/max(np.linalg.norm(table_vector),1e-8)
        density = spectra["spectra"][0,k].sum(axis=-1)*(360/len(spectra["direction"]))
        frequency = spectra["frequency"]
        variance = np.trapezoid(density,frequency)
        low=float(np.trapezoid(density[:3],frequency[:3])/variance)
        high=float(np.trapezoid(density[-3:],frequency[-3:])/variance)
        peak=int(np.argmax(density))
        sites[name] = dict(longitude_deg=xy[0],latitude_deg=xy[1],depth_m=float(row[5]),
            hs_m=float(row[2]),peak_period_s=float(row[3]),mean_period_s=float(row[4]),
            direction_deg=float(row[6]),breaking_fraction=float(row[7]),
            shoreward_normal=None if normal is None else normal.tolist(),
            transport=power,table_vector_w_m=table_vector.tolist(),
            transport_relative_difference=float(error),transport_agrees=bool(error<.01),
            low_frequency_edge_fraction=low,high_frequency_edge_fraction=high,
            spectral_edges_pass=bool(low<.001 and high<.01 and 0<peak<len(frequency)-1))
    selected = values[wet]
    shallow = selected[:,5]<10
    return dict(name=case["name"],path=str(path.relative_to(ROOT)),sha256=sha256(path),
        hour=case["hour"],day_in_second_cycle=case["day_in_second_cycle"],
        stride=case["stride"],spacing_m=case["spacing_m"],directions=case["directions"],
        frequency_intervals=case["frequency_intervals"],wet_nodes=int(wet.sum()),
        maximum_hs_m=float(selected[:,2].max()),minimum_depth_m=float(selected[:,5].min()),
        nodes_below_ten_m=int(shallow.sum()),nodes_with_one_percent_breaking=int((selected[:,7]>=.01).sum()),
        maximum_breaking_fraction=float(selected[:,7].max()),
        maximum_breaking_loss_w_m2=float(selected[:,10].max()*factor),
        boundary=case["boundary"],stopping=case["stopping"],sites=sites,
        boundary_source=case["boundary_source"],boundary_mapping=case["boundary_mapping"],
        elapsed_wall_s=case["run"]["elapsed_wall_s"])


def compare(a,b):
    """Compare shared nodes at least 2 km inside the smaller boundary."""
    if a["hour"] != b["hour"]:
        raise ValueError("Spatial checks require the same weather phase")
    ax,ay,bx,by = (np.asarray(case[key]) for case,key in
                    ((a,"longitude_deg"),(a,"latitude_deg"),(b,"longitude_deg"),(b,"latitude_deg")))
    ai,bi = np.nonzero(np.isclose(ax[:,None],bx[None,:],rtol=0,atol=1e-10))
    aj,bj = np.nonzero(np.isclose(ay[:,None],by[None,:],rtol=0,atol=1e-10))
    av,bv = np.asarray(a["values"])[np.ix_(aj,ai)],np.asarray(b["values"])[np.ix_(bj,bi)]
    wet = np.asarray(a["wet"])[np.ix_(aj,ai)] & np.asarray(b["wet"])[np.ix_(bj,bi)]
    xx,yy = np.meshgrid(ax[ai],ay[aj])
    east_west=np.minimum(xx-max(ax[0],bx[0]),min(ax[-1],bx[-1])-xx)*np.cos(np.deg2rad(yy))
    north_south=np.minimum(yy-max(ay[0],by[0]),min(ay[-1],by[-1])-yy)
    interior = np.minimum(east_west,north_south)*MOON_RADIUS*np.pi/180 >= 2000
    selected = wet & interior & (av[...,2]>=.1) & (bv[...,2]>=.1)
    if not selected.any():
        raise ValueError("Empty common energetic interior")
    np.testing.assert_allclose(av[selected,5],bv[selected,5],rtol=1e-5,atol=.01)
    h = abs(av[...,2]-bv[...,2])/np.maximum(bv[...,2],.1)
    period = abs(av[...,4]-bv[...,4])/np.maximum(bv[...,4],1)
    flux_a,flux_b = av[...,8:10]*a["energy_flux_conversion"],bv[...,8:10]*b["energy_flux_conversion"]
    flux = np.linalg.norm(flux_a-flux_b,axis=-1)/np.maximum(np.linalg.norm(flux_b,axis=-1),1)
    def metrics(mask):
        if not mask.any():
            return dict(nodes=0)
        return dict(nodes=int(mask.sum()),
            hs_relative_p95=float(np.quantile(h[mask],.95)),hs_relative_max=float(h[mask].max()),
            hs_absolute_max_m=float(abs(av[...,2]-bv[...,2])[mask].max()),
            mean_period_relative_p95=float(np.quantile(period[mask],.95)),
            flux_vector_relative_p95=float(np.quantile(flux[mask],.95)),
            flux_vector_relative_max=float(flux[mask].max()))
    result = dict(from_case=a["name"],to_case=b["name"],boundary_buffer_m=2000,
        minimum_hs_m=.1,flux_denominator_floor_w_m=1.,period_denominator_floor_s=1.,
        all_depths=metrics(selected),
        depth_bands={f"{lo}_{hi}m":metrics(selected & (bv[...,5]>=lo) & (bv[...,5]<hi))
                     for lo,hi in ((0,10),(10,50),(50,200),(200,100000))})
    worst = np.unravel_index(np.argmax(np.where(selected,h,-1)),h.shape)
    result["largest_height_change"] = dict(longitude_deg=float(xx[worst]),latitude_deg=float(yy[worst]),
        depth_m=float(bv[worst][5]),from_hs_m=float(av[worst][2]),to_hs_m=float(bv[worst][2]))
    m = result["all_depths"]
    result["exposure_tolerance_pass"] = bool(m["hs_relative_p95"]<.05 and m["hs_relative_max"]<.15
        and m["flux_vector_relative_p95"]<.10 and m["flux_vector_relative_max"]<.30)
    return result


def analyse(run_root=RUNS):
    run_root=Path(run_root).absolute()
    paths = sorted(run_root.glob("h*/product.json"))
    cases = {p.parent.name:verify_case(p) for p in paths}
    if not cases:
        raise ValueError("Complete a coastal case first")
    summaries = [case_summary(cases[p.parent.name],p) for p in paths]
    comparisons = []
    for a,b in (("h852_s8_d36","h852_s4_d36"),("h852_s4_d36","h852_s2_d36"),
                ("h852_s2_d36","h852_s1_d36"),("h852_s4_d36","h852_s4_d72"),
                ("h852_s4_d36","h852_s4_d36_expanded"),
                ("h852_s4_d36","h852_s4_d36_band"),
                ("h852_s4_d36","h852_s4_d36_accuracy"),
                ("h852_s4_d36","h852_s4_d36_regional_accuracy"),
                ("h852_s4_d36","h852_s4_d36_r8"),
                ("h852_s4_d36","h852_s4_d36_nearest")):
        if a in cases and b in cases:
            comparisons.append(compare(cases[a],cases[b]))
    parent_path = run_root/"parent/product.json"
    parent = json.loads(parent_path.read_text())
    if sha256(REFERENCE/"product.json") != parent["reference_case_sha256"]:
        raise ValueError("The continuous basin reference changed")
    for filename,digest in {**parent["run"]["identity"]["inputs"],**parent["run"]["output_sha256"]}.items():
        if sha256(parent_path.parent/filename) != digest:
            raise ValueError(f"Parent replay output changed: {filename}")
    regions=[]
    for path in sorted(run_root.glob("region_h*/product.json")):
        region=json.loads(path.read_text())
        if region["schema"] != "terluna.climate.shore-region/1":
            raise ValueError("Unsupported regional product")
        for filename,digest in {**region["run"]["identity"]["inputs"],**region["run"]["output_sha256"]}.items():
            if sha256(path.parent/filename) != digest:
                raise ValueError(f"Regional run changed: {path.parent/filename}")
        for stem in ("terrain","parent"):
            if sha256(ROOT/region[stem+"_path"]) != region[stem+"_sha256"]:
                raise ValueError(f"Regional {stem} product changed")
        for stem in ("producer","builder"):
            if sha256(path.parent/(stem+"_source.py")) != region[stem+"_sha256"]:
                raise ValueError(f"Regional {stem} source changed")
        regions.append(dict(path=str(path.relative_to(ROOT)),sha256=sha256(path),
            name=region["name"],hour=region["hour"],stride=region["stride"],stopping=region["stopping"],
            boundary_mapping=region["boundary_mapping"],wet_nodes=region["wet_nodes"],
            elapsed_wall_s=region["run"]["elapsed_wall_s"]))
    fields_path=run_root/"parent_fields/product.json"
    fields=json.loads(fields_path.read_text())
    if fields["reference_case_sha256"] != parent["reference_case_sha256"]:
        raise ValueError("Basin replays use different reference seas")
    for filename,digest in {**fields["run"]["identity"]["inputs"],**fields["run"]["output_sha256"]}.items():
        if sha256(fields_path.parent/filename) != digest:
            raise ValueError(f"Basin field replay changed: {filename}")
    return dict(schema="terluna.climate.shore-exposure/1",
        evidence="Conditional coastal wave exposure from native 118 m flooded lunar terrain and directional spectra sampled during the second solar cycle of the established basin history.",
        reading_rule="Cases hold the selected wind and boundary spectra fixed. Shoreward power integrates the incoming part of the directional spectrum across the local shallower-water normal. These are station-specific fluxes per metre, with separate grid, boundary, spectral and stopping checks. Coastal event durations and beach run-up require further calculations.",
        units=dict(height="m",period="s",power_per_width="W/m",breaking_loss="W/m2",angles="degrees east and north; wave travel direction counterclockwise from east"),
        producer=dict(path="climate/waves/shore_checks.py",sha256=sha256(Path(__file__))),
        parent=dict(path=str(parent_path.relative_to(ROOT)),sha256=sha256(parent_path),
                    reference_values_exact=parent["reference_values_exact"],
                    reference_values_compared=parent["reference_values_compared"]),
        parent_fields=dict(path=str(fields_path.relative_to(ROOT)),sha256=sha256(fields_path),
                           reference_values_exact=fields["reference_values_exact"]),
        site_selection="Bathymetry-selected points on the 474 m grid: nearest water 10–60 m deep to west (93.26 E, 2.36 N), east (93.42 E, 2.35 N), and south (93 E, 1.84 N) targets; offshore water deeper than 500 m nearest (92.85 E, 2.35 N). All finer grids share these physical points.",
        episode_selection=episode_context(),
        exposure_tolerances=dict(hs_relative_p95=.05,hs_relative_max=.15,
                                 flux_vector_relative_p95=.10,flux_vector_relative_max=.30),
        assumptions=["Static flooded LOLA rock at the atlas sea level; source interpolation and the inherited geoid frame approximation remain.",
            "Inherited gravity-scaled SWAN wind growth and breaking parameters, fixed air-density build and spatial surface-density wind correction.",
            "Current-free waves; bottom friction, diffraction, shoreline reflection, beach sediment and erosion are omitted.",
            "Linear spectral transport and phase-averaged breaking; the narrow swash zone requires a separate model."],
        regional_cases=regions,cases=summaries,comparisons=comparisons)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root",type=Path,default=RUNS)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    args = parser.parse_args()
    result = analyse(args.run_root)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(cases=len(result["cases"]),comparisons=len(result["comparisons"]),
                         output=str(args.output))))


if __name__ == "__main__":
    main()
