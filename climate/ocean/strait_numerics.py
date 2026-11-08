"""The strait study's numerics on its first day: each model's step, and SCHISM's filter and backtracking.

Each case runs the strait study's first day (the smoke case: two half-day chunks, snapshots every three hours, the
tide spinning up from rest) in its own folder, runs/<model>/fecunditatis_test_<label>, and is compared with its
reference on the Fecunditatis basin's mean level, the largest speeds in the western channels and the Crisium
passage, and each water point's largest speed. Over the day the straits' currents grow from nothing to about
0.4-0.7 m/s.

SCHISM's default filter (ishapiro = 1, shapiro0 = 0.5) smooths the velocities once a step at a fixed strength, so the
smoothing per unit time falls as the step grows. Its Smagorinsky-like filter (ishapiro = 2) takes the strength
0.5 tanh(g0 dt mu), mu = sqrt(u_x^2 + v_y^2 + (u_y + v_x)^2 / 2), which is proportional to the step while the tanh
is unsaturated. A filter of strength s acting once a step on a grid of spacing L diffuses momentum like a viscosity
s L^2 / (4 dt), so the unsaturated filter acts as g0 mu L^2 / 8; MITgcm's Smagorinsky viscosity is
(viscC2smag / pi)^2 L^2 |D|, with |D| = sqrt(2) mu for a non-divergent flow. Equal at viscC2smag = 2.2, g0 = 5.6.

`python -m climate.ocean.strait_numerics [label ...]` runs the cases (all by default; about 40 minutes on four
threads) and writes results/ocean/strait_numerics.json under climate/.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import replace

import numpy as np

from climate.ocean import fecunditatis as study

SMAG_G0 = study.SCHISM_FILTER_G0                       # 5.55: SCHISM's Smagorinsky-like filter at MITgcm's strength
SHAPIRO = dict(filter_mode=1, filter_strength=0.5, btrack_max_s=30.0)    # SCHISM's sample settings
SMAG = dict(filter_mode=2, filter_strength=SMAG_G0, btrack_max_s=30.0)
CASES = {                                               # label: (model, settings, reference)
    "mitgcm_dt300": ("mitgcm", dict(deltat=300.0), None),
    "mitgcm_dt150": ("mitgcm", dict(deltat=150.0), "mitgcm_dt300"),
    "schism_dt300_shapiro": ("schism", dict(SHAPIRO, dt=300.0), None),
    "schism_dt150_shapiro": ("schism", dict(SHAPIRO, dt=150.0), "schism_dt300_shapiro"),
    "schism_dt600_shapiro": ("schism", dict(SHAPIRO, dt=600.0), "schism_dt300_shapiro"),
    "schism_dt300_smag": ("schism", dict(SMAG, dt=300.0), None),
    "schism_dt150_smag": ("schism", dict(SMAG, dt=150.0), "schism_dt300_smag"),
    "schism_dt600_smag": ("schism", dict(SMAG, dt=600.0), "schism_dt300_smag"),
    "schism_dt300_smag_bt300": ("schism", dict(SMAG, dt=300.0, btrack_max_s=300.0), "schism_dt300_smag"),
}
PATH = study.RESULTS / "strait_numerics.json"


def experiment(label):
    model, settings, _ = CASES[label]
    if model == "mitgcm":
        return replace(study.mitgcm_case(smoke=True)[0], name=f"fecunditatis_test_{label}", **settings)
    return replace(study.schism_case(smoke=True), name=f"fecunditatis_test_{label}", **settings)


def metrics(label):
    model = CASES[label][0]
    with np.load(study.HERE / "runs" / model / f"fecunditatis_test_{label}" / "reduced.npz", allow_pickle=False) as s:
        r = {k: s[k] for k in s.files}
    return dict(t_h=(r["t"] / 3600.0).tolist(), fecunditatis_level_m=r["fec_mean"].tolist(),
                western_channels_m_s=r["zone_speed_0"].tolist(), crisium_passage_m_s=r["zone_speed_1"].tolist(),
                largest_speed_m_s=float(r["speed_max"].max()), p99_largest_speed_m_s=float(np.percentile(r["speed_max"], 99)),
                speed_max=r["speed_max"])


def compare(a, b):
    """Each series' largest difference over the common times, relative to the reference's largest magnitude."""
    common = np.intersect1d(np.round(a["t_h"], 6), np.round(b["t_h"], 6))
    out = {}
    for key in ("fecunditatis_level_m", "western_channels_m_s", "crisium_passage_m_s"):
        ta, tb = np.round(a["t_h"], 6), np.round(b["t_h"], 6)
        xa = np.array(a[key])[np.isin(ta, common)]
        xb = np.array(b[key])[np.isin(tb, common)]
        out[key] = float(np.max(np.abs(xa - xb)) / np.max(np.abs(xb)))
    out["largest_speed_m_s"] = a["largest_speed_m_s"] / b["largest_speed_m_s"] - 1
    out["p99_largest_speed_m_s"] = a["p99_largest_speed_m_s"] / b["p99_largest_speed_m_s"] - 1
    sa, sb = a["speed_max"], b["speed_max"]
    if sa.shape == sb.shape:
        out["rms_largest_speed"] = float(np.sqrt(np.mean((sa - sb) ** 2)) / np.sqrt(np.mean(sb ** 2)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("labels", nargs="*", help="cases to run (all by default): " + ", ".join(CASES))
    args = parser.parse_args()
    unknown = [x for x in args.labels if x not in CASES]
    if unknown:
        parser.error(f"unknown cases: {unknown}")
    record = json.loads(PATH.read_text()) if PATH.exists() else dict(schema="terluna.climate.ocean-strait-numerics/1",
                                                                      smag_g0=SMAG_G0, cases={})
    for label in args.labels or list(CASES):
        model, settings, reference = CASES[label]
        exp = experiment(label)
        start = time.time()
        (study.run_mitgcm if model == "mitgcm" else study.run_schism)(smoke=True, exp=exp)
        m = metrics(label)
        entry = dict(model=model, settings=settings, reference=reference, wall_minutes=(time.time() - start) / 60.0,
                     **{k: v for k, v in m.items() if k != "speed_max"})
        if reference and reference in record["cases"]:
            entry["against_reference"] = compare(m, metrics(reference))
        record["cases"][label] = entry
        PATH.write_text(json.dumps(record, indent=1) + "\n")
        print(label, json.dumps(entry.get("against_reference", {})), f"{entry['wall_minutes']:.1f} min")
    print(PATH)


if __name__ == "__main__":
    main()
