"""Propagate bounded photon-only feedback and measure shadow service."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from shared import constants as K
from protection.dynamics.ephemeris import Ephemeris
from protection.dynamics.cycling import (CyclingEnvironment, dro_seed, to_inertial,
    propagate_tacker, service_fraction, outgoing_clearance)
from .cycling_search import SCENARIO, CENTRAL, RUN, ROOT
from shared.provenance import constants_used


RAW_SOURCES = [
    "protection/dynamics/cycling.py", "protection/dynamics/ephemeris.py",
    "protection/dynamics/ephemeris.json", "protection/dynamics/model.py",
    "protection/dynamics/optical.py", "research/studies/solar_shield_array/scenario.json",
    "research/studies/solar_shield_array/cycling_search.py",
    "research/studies/solar_shield_array/cycling_feedback.py",
]


def source_identity():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in RAW_SOURCES}


def summarize(trace):
    t = trace["t"]
    avg = lambda x: float(np.trapezoid(x, t)/(t[-1]-t[0]))
    q = trace["state"][:, :3]
    r = np.linalg.norm(q, axis=-1)
    service = service_fraction(trace, SCENARIO["protected_radii"]*K.MOON_RADIUS)
    return dict(completed=bool(trace["completed"]), days=float((t[-1]-t[0])/K.JULIAN_DAY),
                moon_range_km=[float(np.min(r)/1000), float(np.max(r)/1000)],
                useful_projected_area_fraction=avg(service),
                eclipse_fraction=avg(1-trace["visibility"]),
                sail_impulse_m_s=float(np.trapezoid(np.linalg.norm(trace["sail"], axis=-1), t)),
                jacobi_range=[float(np.min(trace["jacobi"])), float(np.max(trace["jacobi"]))],
                minimum_reflected_beam_clearance_rad=float(np.min(outgoing_clearance(q, trace,
                    trace["angles"], SCENARIO["protected_radii"]*K.MOON_RADIUS))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-km", type=float, default=50000.)
    parser.add_argument("--days", type=float, default=90.)
    parser.add_argument("--mode", choices=["tacking", "face_on", "coast", "service_coast"], default="tacking")
    parser.add_argument("--width", type=float, default=.01)
    parser.add_argument("--step", type=float, default=1800.)
    parser.add_argument("--spectrum", choices=["filter_only", "full_off_core"], default="full_off_core")
    parser.add_argument("--max-step", type=float, default=None)
    parser.add_argument("--rtol", type=float, default=2e-9)
    parser.add_argument("--perturb", action="store_true", help="Add 1 m in ICRF x and 1 mm/s in y")
    parser.add_argument("--out-name", help="Stable name under the ignored study run directory")
    args = parser.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    before = time.monotonic()
    identity = source_identity()
    duration = args.days*K.JULIAN_DAY
    eph = Ephemeris(epoch=SCENARIO["epoch_tdb"])
    samples = eph.sample(np.arange(-3600, duration+7200, 1800.))
    eph.close()
    env = CyclingEnvironment(samples, SCENARIO["baseline_areal_mass_kg_m2"], CENTRAL, args.spectrum)
    local, _ = dro_seed(args.radius_km*1000)
    initial = to_inertial(local, env.frames(0.), env.frames(0., 1))
    if args.perturb:
        initial += np.array([1., 0., 0., 0., .001, 0.])
    print(json.dumps(dict(stage="propagating", mode=args.mode, days=args.days,
                          radius_km=args.radius_km, width=args.width, spectrum=args.spectrum)), flush=True)
    trace = propagate_tacker(env, initial, duration, width=args.width, mode=args.mode,
                             step=args.step, protected_radius=SCENARIO["protected_radii"]*K.MOON_RADIUS,
                             max_step=args.max_step, rtol=args.rtol)
    name = args.out_name or f"feedback_{args.radius_km:g}_{args.days:g}_{args.mode}_{args.width:g}_{args.step:g}_{args.spectrum}"
    if Path(name).name != name:
        raise ValueError("Output name must be a filename stem")
    path = RUN/f"{name}.npz"
    np.savez_compressed(path, **trace)
    result = {**summarize(trace), "elapsed_s": time.monotonic()-before}
    if identity != source_identity():
        raise ValueError("Numerical sources changed during the run; do not publish these raw bytes")
    metadata = dict(schema="terluna.research.solar-shield-cycling-raw/1", config=vars(args),
                    producer=dict(source_hashes=identity, constants=constants_used(identity)),
                    raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), summary=result)
    path.with_suffix(".json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
