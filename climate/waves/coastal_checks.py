"""Audit nested spectra, terrain exposure and completed coastal wave runs."""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import json
from pathlib import Path

import numpy as np

from climate.waves.coastal import RUNS, TERRAIN, load_terrain, read_coastal
from climate.waves.model import ORIGIN, ROOT, sha256


def nesting_spectra(path):
    """Read SWAN's full directional nesting format, including dry locations."""
    lines = iter(line.strip() for line in path.read_text().splitlines()
                 if line.strip() and not line.startswith("$"))

    def line():
        try:
            return next(lines)
        except StopIteration as exc:
            raise ValueError("Truncated nesting spectrum") from exc

    def expect(keyword):
        if line().split()[0] != keyword:
            raise ValueError(f"Expected nesting field {keyword}")

    expect("SWAN")
    expect("TIME")
    if int(line().split()[0]) != 1:
        raise ValueError("Unsupported nesting time code")
    expect("LONLAT")
    nloc = int(line().split()[0])
    locations = np.array([[float(v) for v in line().split()] for _ in range(nloc)])
    expect("RFREQ")
    nf = int(line().split()[0])
    frequency = np.array([float(line().split()[0]) for _ in range(nf)])
    expect("CDIR")
    nd = int(line().split()[0])
    direction = np.array([float(line().split()[0]) for _ in range(nd)])
    expect("QUANT")
    if int(line().split()[0]) != 1:
        raise ValueError("Expected one nesting quantity")
    expect("VaDens")
    expect("m2/Hz/degr")
    exception = float(line().split()[0])
    if (locations.shape != (nloc, 2) or nf < 3 or nd < 3 or
            not all(np.isfinite(a).all() for a in (locations, frequency, direction)) or
            np.any(frequency <= 0) or np.any(np.diff(frequency) <= 0) or
            not np.allclose(np.diff(direction), 360/nd)):
        raise ValueError("Invalid nesting coordinates")
    dates, spectra, available = [], [], []
    for date_line in lines:
        dates.append(datetime.strptime(date_line.split()[0], "%Y%m%d.%H%M%S"))
        values = np.zeros((nloc, nf, nd))
        wet = np.ones(nloc, bool)
        for i in range(nloc):
            kind = line()
            if kind == "NODATA":
                wet[i] = False
            elif kind == "FACTOR":
                scale = float(line())
                table = np.array([[float(v) for v in line().split()] for _ in range(nf)])
                if (table.shape != (nf, nd) or not np.isfinite(table).all() or
                        np.any(table < 0) or not np.isfinite(scale) or scale <= 0):
                    raise ValueError("Invalid nesting variance")
                values[i] = table * scale
            elif kind != "ZERO":
                raise ValueError(f"Unsupported nesting record {kind}")
        spectra.append(values)
        available.append(wet)
    if not dates or any(b <= a for a, b in zip(dates, dates[1:])) or exception >= 0:
        raise ValueError("Invalid nesting times or exception value")
    return dict(dates=dates, locations=locations, frequency=frequency, direction=direction,
                spectra=np.asarray(spectra), available=np.asarray(available))


def audit(product_path, run_root=RUNS):
    product = json.loads(product_path.read_text())
    if product["schema"] != "terluna.climate.coastal-resolution/1":
        raise ValueError("Unsupported coastal product")
    for name, record in product["runs"].items():
        for filename, digest in record["output_sha256"].items():
            if sha256(run_root/name/filename) != digest:
                raise ValueError(f"Changed coastal output: {name}/{filename}")
    if sha256(TERRAIN) != product["inputs"]["terrain_sha256"]:
        raise ValueError("Terrain differs from the coastal run")
    terrain = load_terrain()
    boundary = nesting_spectra(run_root/"parent/coast.nest")
    if boundary["dates"] != [ORIGIN+timedelta(minutes=15*i) for i in range(49)]:
        raise ValueError("Nesting spectra must cover the full twelve hours every fifteen minutes")
    available = boundary["available"]
    if not np.all(available == available[0]):
        raise ValueError("Changing boundary availability")
    x, y = terrain["longitude_deg"], terrain["latitude_deg"]
    ij = np.rint((boundary["locations"]-[x[0], y[0]])*16).astype(int)
    coastal_wet = terrain["wet"][ij[:, 1], ij[:, 0]]
    missing = coastal_wet & ~available[0]
    west = np.isclose(boundary["locations"][:, 0], x[0])
    frequency = boundary["frequency"]
    density = boundary["spectra"][-1].sum(axis=-1)*(360/len(boundary["direction"]))
    variance = np.trapezoid(density, frequency, axis=-1)
    energetic = variance >= (.1/4)**2
    if not energetic.any():
        raise ValueError("Empty final nesting spectrum")
    low = np.trapezoid(density[:, :3], frequency[:3], axis=-1)[energetic]/variance[energetic]
    high = np.trapezoid(density[:, -3:], frequency[-3:], axis=-1)[energetic]/variance[energetic]
    peak = np.argmax(density[energetic], axis=-1)
    details, cubes, masks = {}, {}, {}
    for case in product["cases"]:
        name = case["name"]
        cube, wet = read_coastal(run_root/name/"coastal.tbl", terrain, case["stride"])
        np.testing.assert_array_equal(cube[-1], case["final_values"])
        cubes[name], masks[name] = cube, wet
        final = cube[-1][wet]
        shallow = final[:, 6] < 10
        details[name] = dict(wet_nodes=int(wet.sum()), minimum_depth_m=float(final[:, 6].min()),
                             nodes_shallower_than_ten_m=int(shallow.sum()),
                             nodes_with_one_percent_breaking=int(np.sum(final[:, 8] >= .01)),
                             maximum_final_hs_m=float(final[:, 3].max()))
    locations = {}
    for a, b, factor in (("coast_s4_d36", "coast_s2_d36", 2),
                         ("coast_s2_d36", "coast_s1_d36", 2),
                         ("coast_s2_d36", "coast_s2_d72", 1)):
        coarse, fine = cubes[a][-1], cubes[b][-1, ::factor, ::factor]
        mask = masks[a] & masks[b][::factor, ::factor] & (coarse[..., 3] >= .1) & (fine[..., 3] >= .1)
        mask[:2] = mask[-2:] = False
        mask[:, :2] = mask[:, -2:] = False
        c, f = coarse[mask], fine[mask]
        index = np.argmax(np.abs(c[:, 3]/f[:, 3]-1))
        locations[a+"_to_"+b] = dict(longitude_deg=float(c[index, 1]), latitude_deg=float(c[index, 2]),
                                     depth_m=float(c[index, 6]), coarse_hs_m=float(c[index, 3]),
                                     refined_hs_m=float(f[index, 3]))
    # The extra spectral output should preserve the original parent calculation.
    pilot = json.loads((ROOT/"climate/waves/results/pilot.json").read_text())
    reference = pilot["runs"]["smythii_east_half_grid"]["output_sha256"]["basin.tbl"]
    same_parent = sha256(run_root/"parent/basin.tbl") == reference
    return dict(schema="terluna.climate.coastal-resolution-checks/1",
                evidence="Output integrity and numerical resolution diagnostics for the conditional eastern Smythii experiment.",
                reading_rule="Boundary spectra retain dry parent locations as NODATA. Wet child locations at those points receive zero incident energy through SWAN's nesting interpolation. The western boundary carries the imposed eastward wind sea. Kilometre terrain leaves most of the narrow surf zone between nodes.",
                producer=dict(domain="climate", path="climate/waves/coastal_checks.py", sha256=sha256(Path(__file__))),
                inputs=dict(coastal_sha256=sha256(product_path)),
                units=dict(heights_and_depths="m", locations="degrees east and north", variance_fractions="dimensionless"),
                boundary=dict(snapshots=len(boundary["dates"]), locations=len(available[0]),
                              frequencies=len(boundary["frequency"]), directions=len(boundary["direction"]),
                              wet_child_locations_with_parent_nodata=boundary["locations"][missing].tolist(),
                              western_wet_locations_with_parent_nodata=int(np.sum(west & missing)),
                              variance_nonnegative=True,
                              final_resolved_hs_range_m=(4*np.sqrt([variance[energetic].min(), variance[energetic].max()])).tolist(),
                              maximum_final_low_two_interval_variance_fraction=float(low.max()),
                              maximum_final_high_two_interval_variance_fraction=float(high.max()),
                              final_spectral_edges_pass=bool(low.max() < .001 and high.max() < .01 and np.all((peak > 0) & (peak < len(frequency)-1))),
                              spectral_reading_rule="Resolved printed variance integrated over direction and frequency at final wet locations with Hs >= 0.1 m. SWAN's bulk Hs also includes its diagnostic high-frequency tail."),
                parent_matches_prior_half_degree_table=same_parent,
                cases=details, maximum_height_difference_locations=locations,
                resolution=product["checks"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--product", type=Path, default=ROOT/"climate/waves/results/coastal.json")
    p.add_argument("--run-root", type=Path, default=RUNS)
    p.add_argument("--output", type=Path, default=ROOT/"climate/waves/results/coastal_checks.json")
    args = p.parse_args()
    result = audit(args.product, args.run_root)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
