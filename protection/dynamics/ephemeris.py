"""Pinned DE440s geometric states, without a circular-orbit fallback."""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import urllib.request

import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from shared import constants as K

MANIFEST = Path(__file__).with_name("ephemeris.json")
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_KERNEL = ROOT / "research/runs/solar_shield_array/inputs/de440s.bsp"
BODY_GM = {"sun": K.SUN_GM, "earth": K.EARTH_GM, **K.SOLAR_SYSTEM_GM}
BODY_IDS = {"sun": 10, "earth": 399, "moon": 301, "mercury": 1,
            "venus": 2, "mars": 4, "jupiter": 5, "saturn": 6,
            "uranus": 7, "neptune": 8, "pluto": 9}


def verify_kernel(path: Path) -> dict:
    manifest = json.loads(MANIFEST.read_text())
    data = path.read_bytes()
    if len(data) != manifest["bytes"] or hashlib.sha256(data).hexdigest() != manifest["sha256"]:
        raise ValueError("DE440s kernel size or SHA-256 mismatch; no substitution allowed")
    return manifest


def fetch(path: Path = DEFAULT_KERNEL) -> Path:
    if path.exists():
        verify_kernel(path)
        return path
    manifest = json.loads(MANIFEST.read_text())
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".download")
    try:
        with urllib.request.urlopen(manifest["url"], timeout=45) as source:
            temporary.write_bytes(source.read())
        verify_kernel(temporary)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def calendar_jd_tdb(value: str) -> float:
    """Interpret a calendar date as TDB, explicitly without UTC/leap conversion."""
    from jplephem.calendar import compute_julian_date
    dt = datetime.fromisoformat(value)
    if dt.tzinfo is not None:
        raise ValueError("Specify an unzoned TDB calendar date, rather than UTC")
    return float(compute_julian_date(dt.year, dt.month, dt.day)) + (
        dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6
    ) / K.JULIAN_DAY


class Ephemeris:
    def __init__(self, path: Path = DEFAULT_KERNEL, epoch: str = "2026-10-04"):
        from jplephem.spk import SPK
        self.manifest = verify_kernel(path)
        self.kernel = SPK.open(str(path))
        self.jd0 = calendar_jd_tdb(epoch)
        self.epoch = epoch

    def close(self):
        self.kernel.close()

    def state(self, name: str, t):
        """Barycentric position (m) and velocity (m/s), scalar or vector times."""
        code = BODY_IDS[name]
        if code in (399, 301):
            p0, v0 = self.kernel[0, 3].compute_and_differentiate(self.jd0, np.asarray(t) / K.JULIAN_DAY)
            p1, v1 = self.kernel[3, code].compute_and_differentiate(self.jd0, np.asarray(t) / K.JULIAN_DAY)
            p, v = p0 + p1, v0 + v1
        else:
            p, v = self.kernel[0, code].compute_and_differentiate(self.jd0, np.asarray(t) / K.JULIAN_DAY)
        return np.moveaxis(p, 0, -1) * 1000, np.moveaxis(v, 0, -1) * (1000 / K.JULIAN_DAY)

    def acceleration(self, name: str, t, step: float = 30):
        """Centred velocity difference; convergence checked against 10/90 s."""
        return (self.state(name, np.asarray(t) + step)[1] -
                self.state(name, np.asarray(t) - step)[1]) / (2 * step)

    def sample(self, t, derivative_step: float = 30) -> dict:
        t = np.atleast_1d(t).astype(float)
        mp, mv = self.state("moon", t)
        ma = self.acceleration("moon", t, derivative_step)
        positions, velocities = {}, {}
        for name in BODY_GM:
            p, v = self.state(name, t)
            positions[name], velocities[name] = p - mp, v - mv
        sa = self.acceleration("sun", t, derivative_step) - ma
        ea = self.acceleration("earth", t, derivative_step) - ma
        return {"t": t, "positions": positions, "velocities": velocities,
                "moon_v": mv, "sun_v": mv+velocities["sun"],
                "earth_v": mv+velocities["earth"], "sun_a": ma+sa, "earth_a": ma+ea,
                "moon_a": ma, "sun_relative_a": sa, "earth_relative_a": ea}


class InterpolatedEphemeris:
    """Bounded Hermite state interpolation for trajectory propagation."""
    def __init__(self, samples: dict):
        t = samples["t"]
        self.limits = (t[0], t[-1])
        self.positions = {name: CubicHermiteSpline(t, p, samples["velocities"][name], extrapolate=False)
                          for name, p in samples["positions"].items()}
        self.moon_a = CubicSpline(t, samples["moon_a"], extrapolate=False)
        self.moon_v = CubicSpline(t, samples["moon_v"], extrapolate=False)
        self.sun_v = CubicSpline(t, samples["sun_v"], extrapolate=False)
        self.earth_v = CubicSpline(t, samples["earth_v"], extrapolate=False)
        self.sun_a = CubicSpline(t, samples["sun_a"], extrapolate=False)
        self.earth_a = CubicSpline(t, samples["earth_a"], extrapolate=False)

    def at(self, t):
        if not self.limits[0] <= t <= self.limits[1]:
            raise ValueError("Propagation outside sampled ephemeris")
        return {"positions": {k: s(t) for k, s in self.positions.items()},
                "moon_a": self.moon_a(t), "moon_v": self.moon_v(t), "sun_v": self.sun_v(t),
                "earth_v":self.earth_v(t), "sun_a":self.sun_a(t), "earth_a":self.earth_a(t)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--kernel", type=Path, default=DEFAULT_KERNEL)
    args = parser.parse_args()
    if args.download:
        fetch(args.kernel)
    print(json.dumps(verify_kernel(args.kernel), indent=2))
