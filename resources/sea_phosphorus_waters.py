"""Which lakes the land's water crosses on its way to the seas, and how deep the seas lie.

    python -m resources.sea_phosphorus_waters     # results/sea_phosphorus_waters.json (1-3 minutes, 4 GB)

The run repeats the geography drainage routing (geography.drainage.compute with the climatology the committed
drainage product names) and checks that it reproduces geography/results/drainage.json. Tracers then follow the
routing downstream:
- the water that has crossed at least one rain-fed lake, by the lake's residence on its outflow, each lake counted
  once where water leaves it (the drainage product's lakes whose depressions spill, joined across the 0/360 degree
  seam as the routing is);
- phosphorus in solution carried in proportion to land runoff and kept by each lake as Vollenweider's (1976) relation
  keeps total phosphorus, 1 - R = 1 / (1 + sqrt(tau_w));
- phosphorus on particles kept by each lake on Siyam's (2000) form of Brune's (1953) trap-efficiency curves,
  E = exp(-beta I/C), at Earth's settling and with settling slowed by lunar gravity (I/C scaled by g_earth/g_moon).
The atlas grid gives each receiving water body and the seas' depths. The grid inputs are ignored by Git; the product
records their hashes. resources/sea_phosphorus.py reads the product.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
OUT = HERE / "results/sea_phosphorus_waters.json"
SCHEMA = "terluna.resources.sea-phosphorus-waters/1"
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
DRAINAGE = ROOT / "geography/results/drainage.json"
ATLAS = ROOT / "geography/results/atlas.json"
ATLAS_GRID = ROOT / "geography/products/atlas_28pct_4ppd.npz"
IGNORED = ("climate/gcm/products/climatology_A28_dim5_moon.npz", "geography/inputs/ldem_16.img",
           "geography/inputs/ldem_16.lbl", "geography/inputs/ldem_4.img", "geography/inputs/ldem_4.lbl",
           "geography/inputs/jggrx_0420a_sha.tab", "geography/inputs/jggrx_0420a_sha.lbl",
           "geography/products/atlas_28pct_4ppd.npz")
CODE = ("resources/sea_phosphorus_waters.py", "geography/drainage.py", "geography/topography.py",
        "geography/atlas.py")
# Siyam (2000) exp(-beta I/C) with beta for Brune's (1953) upper, median and lower curves (Eizel-Din et al. 2010).
BRUNE_BETA = dict(upper=0.0055, median=0.0079, lower=0.015)
RESIDENCE_THRESHOLDS_YR = (0.0, 0.1, 1.0, 10.0, 100.0)
DEPTH_BINS_M = (0.0, 50.0, 200.0, 1000.0, 3000.0, float("inf"))
EVIDENCE = ("The geography drainage routing rerun with its committed climatology (reproduced to the cubic metre "
            "per second), with tracers carried down the same receivers: water that has crossed a rain-fed lake, "
            "dissolved P kept by each lake as Vollenweider's (1976) relation keeps total P in Earth's lakes, and "
            "particulate P kept on Siyam's (2000) form of Brune's (1953) reservoir curves. P loads follow land runoff "
            "at one concentration. Steady annual flow; lakes do not fill with sediment, and every cell below the "
            "common level is an outlet, as in the drainage product. The atlas grid gives depths below the common "
            "level and the receiving water body.")
READING_RULE = ("Retention bodies are the drainage product's lakes whose depressions spill, 8-connected and joined "
                "across the 0/360 degree seam, with the lakes of closed depressions cleared from them; retention "
                "applies once, where water leaves a body, and tau is a body's volume over its outflow in years. "
                "Closed depressions keep all they receive. Shares are fractions of the land-runoff P load (P "
                "tracers), of the routed discharge to standing water at the common level, rain on lakes included "
                "(lake_passage_share_of_discharge), or of the land's own runoff that reaches it, carried through "
                "spilling lakes (lake_passage_share_of_land_runoff). Particulate survival keys name the Brune curve "
                "(upper, median, lower) and the settling: earth, or moon with I/C scaled by g_earth/g_moon. "
                "Hypsometry areas are km2 of water at the common level by depth band in metres. by_receiving_water "
                "labels outlets by the 4 pixel per degree atlas: a named body, 'small waters' (other common-level "
                "depressions) or 'below the atlas grid' (16 px/deg cells below the common level inside 4 px/deg "
                "land).")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lake_survival_dissolved(tau_yr):
    """Share of inflowing total P that leaves a lake of residence tau_yr (Vollenweider 1976, OECD form)."""
    tau = np.asarray(tau_yr, float)
    if np.any(tau < 0):
        raise ValueError("Negative residence")
    return 1.0 / (1.0 + np.sqrt(tau))


def lake_survival_particulate(tau_yr, beta, gravity_ratio=1.0):
    """Share of inflowing sediment that passes a lake: 1 - exp(-beta I/C) (Siyam 2000), I/C slowed by gravity."""
    tau = np.asarray(tau_yr, float)
    if np.any(tau < 0) or beta <= 0 or not 0 < gravity_ratio <= 1:
        raise ValueError("Invalid trap-efficiency input")
    with np.errstate(divide="ignore"):
        inflow_over_capacity = np.where(tau > 0, 1.0 / np.where(tau > 0, tau, 1.0), np.inf) / gravity_ratio
    return 1.0 - np.exp(-beta * inflow_over_capacity)


def label_wrapped(mask):
    """8-connected components of a latitude-longitude mask, joined across the 0/360 degree seam as the routing's
    neighbours are; rows do not join over the poles, as in the routing. Returns labels, count and seam joins."""
    from scipy import ndimage
    lab, n = ndimage.label(mask, structure=np.ones((3, 3), bool))
    rows, cols = mask.shape
    parent = np.arange(n + 1)

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    joins = 0
    for dr in (-1, 0, 1):
        r0 = np.arange(max(0, -dr), rows - max(0, dr))
        a, b = lab[r0, 0], lab[r0 + dr, cols - 1]
        both = (a > 0) & (b > 0)
        for x, y in zip(a[both].tolist(), b[both].tolist()):
            rx, ry = root(x), root(y)
            if rx != ry:
                parent[max(rx, ry)] = min(rx, ry)
                joins += 1
    roots = np.array([root(i) for i in range(n + 1)])
    _, new = np.unique(roots, return_inverse=True)
    return new[lab], int(new.max()), joins


def lake_bodies(r, shape):
    """Retention bodies: the drainage product's lakes whose depressions spill, joined across the seam, with the lakes
    of closed depressions cleared from them; their volume, outflow and residence, and a census of the lakes.

    The routing chains the depressions inside one merged lake through their spill points, so retention is applied
    once, where water leaves a body, with the body's own volume over its outflow. A closed depression keeps all it
    receives at its floor, so its lake needs no body of its own."""
    lake, label, area = r["lake"], r["label"], r["area"]
    spill_node = r["passed_on"] > 0
    spilling_cell = lake & spill_node[label]
    closed_cell = lake & ~spill_node[label]
    body, n_body, body_joins = label_wrapped(spilling_cell.reshape(shape))
    body = body.ravel()
    final, sinks, outlet = r["receiver"], r["sinks"], r["outlet"]
    out_all = r["discharge"].copy()
    out_all[sinks] = np.maximum(out_all[sinks] + r["lake_net"][1:], 0.0)
    exits = (body > 0) & ~outlet & (body[final] != body)
    outflow = np.bincount(body[exits], weights=out_all[exits], minlength=n_body + 1)
    depth = np.zeros(lake.size)
    depth[spilling_cell] = r["lake_level"][spilling_cell] - r["height"][spilling_cell]
    volume = np.bincount(body[spilling_cell], weights=(depth * area)[spilling_cell], minlength=n_body + 1)
    spills = outflow > 0
    spills[0] = False
    tau = np.where(spills, volume / np.where(spills, outflow * YEAR_S, 1.0), np.inf)
    lakes, n_lakes, lake_joins = label_wrapped(lake.reshape(shape))
    lakes = lakes.ravel()
    has_spill = np.bincount(lakes[spilling_cell], minlength=n_lakes + 1) > 0
    has_closed = np.bincount(lakes[closed_cell], minlength=n_lakes + 1) > 0
    mixed = has_spill & has_closed
    mixed[0] = False
    _, n_closed, _ = label_wrapped(closed_cell.reshape(shape))
    census = dict(lakes_seam_joined=n_lakes, seam_joins=lake_joins, closed_lakes=n_closed,
                  lakes_with_spilling_and_closed_depressions=int(mixed.sum()),
                  closed_lake_km2_cleared_from_spilling_lakes=float(area[closed_cell & mixed[lakes]].sum() / 1e6),
                  retention_bodies=n_body, retention_bodies_spilling=int(spills.sum()),
                  retention_body_seam_joins=body_joins,
                  spilling_lake_km2=float(area[spilling_cell].sum() / 1e6),
                  closed_lake_km2=float(area[closed_cell].sum() / 1e6))
    return body, n_body, outflow, volume, tau, spills, census


def trace(r, runoff_mm):
    """Carry water and P tracers down the drainage receivers, in the routing's own topological order."""
    from geography import drainage as gd
    final, label, outlet, area = r["receiver"], r["label"], r["outlet"], r["area"]
    n = final.size
    hf = r["height"]
    lake_full = ~outlet & (hf < r["level"][label])
    runoff = np.where(outlet | lake_full, 0.0, np.maximum(runoff_mm.ravel(), 0.0) * gd.MM_DAY * area)
    sinks = r["sinks"]
    node_of = np.full(n, -1, np.int64)
    node_of[outlet] = 0
    node_of[sinks] = np.arange(1, sinks.size + 1)
    sink = np.zeros(n, bool)
    sink[sinks] = True
    node_spills = r["passed_on"] > 0
    body, n_body, outflow, volume, tau, spills, census = lake_bodies(r, r["shape"])
    gravity = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    surv_d = np.where(spills, lake_survival_dissolved(np.where(spills, tau, 0.0)), 0.0)
    surv_p = {f"{k}_{g}": np.where(spills, lake_survival_particulate(np.where(spills, tau, 0.0), b, gr), 0.0)
              for k, b in BRUNE_BETA.items() for g, gr in (("earth", 1.0), ("moon", gravity))}
    acc = runoff.copy()
    passed = {t: np.zeros(n) for t in RESIDENCE_THRESHOLDS_YR}
    land = runoff.copy()                     # the land's own runoff, through spilling lakes unchanged
    land_passed = {t: np.zeros(n) for t in RESIDENCE_THRESHOLDS_YR}
    p_d = runoff.copy()
    p_p = {k: runoff.copy() for k in surv_p}
    donor = ~outlet
    indeg = np.bincount(final[donor], minlength=n).astype(np.int64)
    front = np.flatnonzero((indeg == 0) & donor)
    while front.size:
        out = acc[front].copy()
        at_sink = sink[front]
        k = node_of[front[at_sink]]
        if at_sink.any():
            out[at_sink] = np.maximum(out[at_sink] + r["lake_net"][k], 0.0)
        kept = np.ones(front.size)                       # water and P stop at the floor of a closed depression
        if at_sink.any():
            kept[at_sink] = node_spills[k]
        dest = final[front]
        b = body[front]
        leaving = (b > 0) & (body[dest] != b)
        land_out = land[front] * kept
        for t in RESIDENCE_THRESHOLDS_YR:
            mark = leaving & (tau[b] >= t)
            np.add.at(passed[t], dest, np.where(mark, out, np.minimum(passed[t][front], out)))
            np.add.at(land_passed[t], dest, np.where(mark, land_out, np.minimum(land_passed[t][front], land_out)))
        np.add.at(land, dest, land_out)
        d_out = p_d[front] * kept * np.where(leaving, surv_d[b], 1.0)
        np.add.at(p_d, dest, d_out)
        for key in p_p:
            np.add.at(p_p[key], dest, p_p[key][front] * kept * np.where(leaving, surv_p[key][b], 1.0))
        np.add.at(acc, dest, out)
        np.subtract.at(indeg, dest, 1)
        cand = np.unique(dest)
        front = cand[(indeg[cand] == 0) & donor[cand]]
    return dict(acc=acc, runoff=runoff, passed=passed, land=land, land_passed=land_passed, p_d=p_d, p_p=p_p,
                body_tau=tau, body_spills=spills, body_outflow=outflow, body_volume=volume, census=census,
                depression_spills=int(node_spills[1:].sum()))


def receiving_waters(outlet_cells, cols16, weights, atlas):
    """Sum weights at outlet cells by the 4 px/deg atlas water body that holds them."""
    z = np.load(ATLAS_GRID)
    labels, lat4 = z["water_label"], z["lat_deg"]
    if not lat4[0] > lat4[-1]:
        raise ValueError("Atlas grid rows must run from the north")
    row4 = np.clip(outlet_cells // cols16 // 4, 0, labels.shape[0] - 1)
    col4 = np.clip(outlet_cells % cols16 // 4, 0, labels.shape[1] - 1)
    lab = labels[row4, col4].astype(np.int64)
    missing = lab == 0
    if missing.any():
        best = np.zeros(int(missing.sum()), np.int64)
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                rows = np.clip(row4[missing] + dr, 0, labels.shape[0] - 1)
                near = labels[rows, (col4[missing] + dc) % labels.shape[1]]
                best = np.where(best == 0, near, best)
        lab[missing] = best
    named = {b["label"]: (b["water_names"][0] if b["water_names"] else f"unnamed body {b['label']}")
             for b in atlas["bodies"]}
    out = {}
    for key, w in weights.items():
        sums = np.bincount(lab, weights=w, minlength=int(lab.max()) + 1)
        for L in np.flatnonzero(sums != 0):
            name = named.get(int(L), "small waters" if L else "below the atlas grid")
            out.setdefault(name, {})
            out[name][key] = out[name].get(key, 0.0) + float(sums[L])
    return out


def river_mouths(drainage, atlas):
    """The atlas water body at each of the drainage product's largest river mouths, with its area."""
    from geography.topography import Grid
    z = np.load(ATLAS_GRID)
    labels = z["water_label"]
    grid = Grid(z["lat_deg"], z["lon_deg"], 4)
    areas = np.bincount(labels.ravel(), weights=grid.cell_area_m2.ravel() / 1e6)
    named = {b["label"]: (b["water_names"][0] if b["water_names"] else f"unnamed body {b['label']}")
             for b in atlas["bodies"]}
    rows = []
    for river in drainage["rivers"]["largest"]:
        lat, lon = river["mouth"]
        i = int(np.clip(round((90.0 - lat) * 4 - 0.5), 0, labels.shape[0] - 1))
        j = int(round((lon % 360.0) * 4 - 0.5)) % labels.shape[1]
        window = labels[max(i - 2, 0):i + 3, [(j + d) % labels.shape[1] for d in range(-2, 3)]]
        found = sorted({int(x) for x in window.ravel() if x})
        rows.append(dict(mouth=river["mouth"], discharge_m3s=river["discharge_m3s"],
                         waters=[dict(label=L, name=named.get(L, "small water"), area_km2=float(areas[L]))
                                 for L in found]))
    return rows


def hypsometry(atlas):
    """Area of standing water at the common level by depth band, for the named bodies and in all."""
    from geography.topography import Grid
    z = np.load(ATLAS_GRID)
    grid = Grid(z["lat_deg"], z["lon_deg"], 4)
    area_km2 = grid.cell_area_m2 / 1e6
    level = float(atlas["sea_level_m"])
    water = z["water_label"] > 0
    raw = level - z["height_m"].astype(float)
    depth = np.maximum(raw, 0.0)          # polar cells labelled from the 16 px/deg caps can sit above level here
    named = {b["label"]: (b["water_names"][0] if b["water_names"] else f"unnamed body {b['label']}")
             for b in atlas["bodies"]}
    bands = [f"{lo:g}-{hi:g}" for lo, hi in zip(DEPTH_BINS_M[:-1], DEPTH_BINS_M[1:])]

    def table(mask):
        d = depth[mask]
        a = area_km2[mask]
        rows = {b: float(a[(d >= lo) & (d < hi)].sum())
                for b, lo, hi in zip(bands, DEPTH_BINS_M[:-1], DEPTH_BINS_M[1:])}
        return dict(area_km2=float(a.sum()), volume_km3=float((a * d).sum() / 1e3), by_depth_m_km2=rows)

    out = dict(all_common_level_water=table(water),
               water_cells_above_level_km2=float(area_km2[water & (raw < 0)].sum()))
    for L, name in named.items():
        out[name] = table(z["water_label"] == L)
    return out


def run():
    from geography import drainage as gd
    committed = json.loads(DRAINAGE.read_text())
    atlas = json.loads(ATLAS.read_text())
    climatology = ROOT / next(iter(k for k in committed["sources"] if k.startswith("climate/")))
    if digest(climatology) != committed["sources"][str(climatology.relative_to(ROOT))]:
        raise ValueError("Climatology does not match the drainage product")
    atlas_, _, _, level, h, lat, _, runoff_mm, r = gd.compute(climatology)
    summary = gd.summarise(atlas_, level, h, lat, r)
    reproduced = dict(
        discharge_to_sea_m3s=[summary["rivers"]["discharge_to_sea_m3s"], committed["rivers"]["discharge_to_sea_m3s"]],
        lake_area_share=[summary["lakes"]["area_share"], committed["lakes"]["area_share"]],
        lake_volume_km3=[summary["lakes"]["volume_km3"], committed["lakes"]["volume_km3"]],
        lake_bodies=[summary["lakes"]["bodies"], committed["lakes"]["bodies"]],
        largest_river_m3s=[summary["rivers"]["largest"][0]["discharge_m3s"],
                           committed["rivers"]["largest"][0]["discharge_m3s"]])
    for key, (mine, theirs) in reproduced.items():
        if mine != theirs:
            raise ValueError(f"Routing does not reproduce the drainage product: {key} {mine} != {theirs}")
    r = dict(r, height=h.ravel().astype(np.float64), shape=h.shape)
    t = trace(r, runoff_mm)
    outlet = r["outlet"]
    if np.max(np.abs(t["acc"] - r["discharge"])) > 1e-6 * np.max(r["discharge"]):
        raise ValueError("Tracer routing departs from the drainage discharge")
    land = float(t["runoff"].sum())
    to_sea = float(t["acc"][outlet].sum())
    cells = np.flatnonzero(outlet)
    weights = dict(discharge_m3s=t["acc"][cells], lake_passed_m3s=t["passed"][0.0][cells],
                   dissolved_p_share=t["p_d"][cells] / land,
                   particulate_p_share_median_moon=t["p_p"]["median_moon"][cells] / land)
    spill = t["body_spills"]
    w, tau = t["body_outflow"][spill], t["body_tau"][spill]
    order = np.argsort(tau)
    cum = np.cumsum(w[order]) / w.sum()
    producer_files = {p: digest(ROOT / p) for p in CODE}
    return dict(
        schema=SCHEMA,
        producer=dict(files=producer_files,
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in (DRAINAGE, ATLAS)},
                      ignored_inputs={p: digest(ROOT / p) for p in IGNORED},
                      constants=constants_used([ROOT / p for p in CODE])),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        reproduced_drainage=reproduced,
        water_m3s=dict(land_runoff=land, lake_rain_less_evaporation=to_sea - land, to_common_level_water=to_sea),
        lake_passage_share_of_discharge={f"tau_ge_{x:g}_yr": float(t["passed"][x][outlet].sum() / to_sea)
                                         for x in RESIDENCE_THRESHOLDS_YR},
        land_runoff_reaching_common_level_m3s=float(t["land"][outlet].sum()),
        lake_passage_share_of_land_runoff={f"tau_ge_{x:g}_yr": float(t["land_passed"][x][outlet].sum()
                                                                     / t["land"][outlet].sum())
                                           for x in RESIDENCE_THRESHOLDS_YR},
        lakes=dict(drainage_product=summary["lakes"]["bodies"], **t["census"],
                   spilling_depressions=t["depression_spills"]),
        outflow_weighted_tau_yr={f"p{round(q * 100)}": float(tau[order][min(np.searchsorted(cum, q), tau.size - 1)])
                                 for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
        dissolved_p_survival=float(t["p_d"][outlet].sum() / land),
        particulate_p_survival={k: float(v[outlet].sum() / land) for k, v in t["p_p"].items()},
        brune_beta=BRUNE_BETA,
        gravity_ratio=MOON_SURFACE_GRAVITY / STANDARD_GRAVITY,
        by_receiving_water=receiving_waters(cells, h.shape[1], weights, atlas),
        largest_rivers_receiving_water=river_mouths(committed, atlas),
        hypsometry=hypsometry(atlas))


def main():
    product = run()
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(product, indent=2, allow_nan=False) + "\n")
    print(OUT.relative_to(ROOT))


if __name__ == "__main__":
    sys.exit(main())
