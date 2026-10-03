"""Bounded SWAN runs with explicit forcing, provenance and completion checks.

These are uniform-wind, flat-bottom experiments. A wind percentile selects a
forcing scenario; it does not assign a percentile to the resulting waves.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from datetime import datetime, timedelta
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import time

import numpy as np

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "research/runs/waves"
ORIGIN = datetime(2000, 1, 1)  # SWAN clock only; no ephemeris in these cases.
TABLE_COLUMNS = ("time_s", "fetch_m", "hs_m", "peak_period_s", "mean_period_s", "depth_m")


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


@dataclass(frozen=True)
class Case:
    name: str
    gravity_m_s2: float
    wind_m_s: float
    length_m: float = 100_000.0
    depth_m: float = 1000.0
    cells: int = 100
    direction_bins: int = 36
    frequency_intervals: int = 48
    earth_low_hz: float = 0.03
    earth_high_hz: float = 3.0
    steps: int = 576
    step_s: int = 300
    output_every_steps: int = 12
    agrow: bool = True
    seed_height_m: float = 0.0
    seed_period_s: float = 2.0
    water_density_kg_m3: float = 1025.0  # Stated seawater scenario, not measured lunar water.

    def __post_init__(self):
        if not re.fullmatch(r"[a-z0-9_]+", self.name):
            raise ValueError("Case names must use lowercase letters, digits and underscores")
        for name in ("gravity_m_s2", "wind_m_s", "length_m", "depth_m", "earth_low_hz",
                     "earth_high_hz", "seed_period_s", "water_density_kg_m3"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("cells", "direction_bins", "frequency_intervals", "steps", "step_s", "output_every_steps"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if self.direction_bins < 12 or self.direction_bins % 4:
            raise ValueError("At least three directional bins per quadrant are required")
        if self.frequency_intervals < 3 or self.earth_low_hz >= self.earth_high_hz:
            raise ValueError("Invalid frequency grid")
        if self.steps % self.output_every_steps:
            raise ValueError("Final time must be an output time")
        if not math.isfinite(self.seed_height_m) or self.seed_height_m < 0:
            raise ValueError("Initial wave height must be finite and nonnegative")
        if not self.agrow and self.seed_height_m == 0:
            raise ValueError("A zero initial spectrum needs AGROW to seed wind waves")

    @property
    def duration_s(self) -> int:
        return self.steps * self.step_s

    @property
    def frequencies_hz(self) -> tuple[float, float]:
        ratio = self.gravity_m_s2 / STANDARD_GRAVITY
        return self.earth_low_hz * ratio, self.earth_high_hz * ratio

    def input_text(self) -> str:
        low, high = self.frequencies_hz
        start = ORIGIN.strftime("%Y%m%d.%H%M%S")
        end = (ORIGIN + timedelta(seconds=self.duration_s)).strftime("%Y%m%d.%H%M%S")
        seed = "INIT ZERO" if not self.seed_height_m else (
            f"INIT PAR {self.seed_height_m:.12g} {self.seed_period_s:.12g} 0 20")
        return f"""PROJECT 'Terluna waves' '001'
MODE NONSTATIONARY ONEDIMENSIONAL
COORDINATES CARTESIAN
SET GRAV {self.gravity_m_s2:.12g} RHO {self.water_density_kg_m3:.12g}
CGRID REGULAR 0 0 0 {self.length_m:.12g} 0 {self.cells} 0 &
 CIRCLE {self.direction_bins} {low:.12g} {high:.12g} {self.frequency_intervals}
INPGRID BOTTOM 0 0 0 1 0 {self.length_m:.12g} 1
READINP BOTTOM 1. 'depth.bot' 1 0 FREE
WIND {self.wind_m_s:.12g} 0
GEN3 KOMEN DRAG WU{' AGROW' if self.agrow else ''}
BREAKING CONSTANT 1.0 0.73
PROP BSBT
{seed}
QUANTITY TSEC REF {start}
CURVE 'line' 0 0 {self.cells} {self.length_m:.12g} 0
TABLE 'line' NOHEAD 'profile.tbl' TSEC XP HS RTP TM01 DEP &
 OUTPUT {start} {self.output_every_steps * self.step_s} SEC
TABLE 'line' NOHEAD 'breaking.tbl' TSEC XP QB &
 OUTPUT {start} {self.output_every_steps * self.step_s} SEC
POINTS 'end' {self.length_m:.12g} 0
SPEC 'end' SPEC1D ABS 'end.spc' &
 OUTPUT {start} {self.output_every_steps * self.step_s} SEC
COMPUTE NONSTAT {start} {self.step_s} SEC {end}
STOP
"""


def read_profile(path: Path, case: Case) -> np.ndarray:
    """Require every requested output, including final time, even if SWAN exits 0."""
    values = np.loadtxt(path, ndmin=2)
    expected_times = np.arange(0, case.duration_s + 1,
                               case.output_every_steps * case.step_s, dtype=float)
    expected_rows = len(expected_times) * (case.cells + 1)
    if values.shape != (expected_rows, len(TABLE_COLUMNS)) or not np.isfinite(values).all():
        raise RuntimeError(f"Incomplete or non-finite SWAN table: {path}")
    cube = values.reshape(len(expected_times), case.cells + 1, len(TABLE_COLUMNS))
    # SWAN's standard text format keeps about five significant digits.
    if not np.allclose(cube[:, :, 0], expected_times[:, None], rtol=6e-5, atol=0.1):
        raise RuntimeError("SWAN output time coverage differs from the requested run")
    locations = np.linspace(0, case.length_m, case.cells + 1)
    if not np.allclose(cube[:, :, 1], locations[None, :], rtol=6e-5, atol=0.01):
        raise RuntimeError("SWAN output locations differ from the requested grid")
    if np.any(cube[:, :, 2] < 0) or not np.allclose(cube[:, :, 5], case.depth_m, rtol=6e-5):
        raise RuntimeError("Negative heights or unexpected depths in SWAN output")
    if np.any((cube[1:, :, 2:3] > 0) & (cube[1:, :, 3:5] <= 0)):
        raise RuntimeError("Missing or invalid wave periods after initialization")
    return cube


def read_spectral_edges(path: Path, case: Case) -> dict:
    """Read final omnidirectional spectrum; report resolved variance near band edges.

    SWAN's integral wave height includes its diagnostic high-frequency tail;
    these fractions refer only to the explicitly resolved spectral band.
    """
    lines = path.read_text().splitlines()
    dates = [line.split()[0] for line in lines if re.match(r"^\d{8}\.\d{6}\s", line)]
    expected_dates = [(ORIGIN + timedelta(seconds=t)).strftime("%Y%m%d.%H%M%S")
                      for t in range(0, case.duration_s + 1, case.output_every_steps * case.step_s)]
    if dates != expected_dates:
        raise RuntimeError("SWAN spectrum time coverage differs from the requested run")
    marker = next(i for i, line in enumerate(lines)
                  if line.split() and line.split()[0] in ("AFREQ", "RFREQ"))
    count = int(lines[marker + 1].split()[0])
    frequency = np.array([float(line.split()[0]) for line in lines[marker + 2:marker + 2 + count]])
    if (count < 3 or len(frequency) != count or not np.isfinite(frequency).all()
            or np.any(frequency <= 0) or np.any(np.diff(frequency) <= 0)):
        raise RuntimeError("Invalid SWAN spectral frequency grid")
    locations = [i for i, line in enumerate(lines)
                 if line.split() and line.split()[0] == "LOCATION"]
    if len(locations) != len(expected_dates):
        raise RuntimeError("SWAN spectrum must have one LOCATION per requested time")
    first = locations[-1] + 1
    density = np.array([float(line.split()[0]) for line in lines[first:first + count]])
    quantity = next(i for i, line in enumerate(lines) if line.startswith("VaDens"))
    exception = float(lines[quantity + 2].split()[0])
    censored = density == exception
    if len(density) != count or np.any((density < 0) & ~censored) or not np.isfinite(density).all():
        raise RuntimeError("Invalid final variance spectrum")
    # SWCMSP (swanout2.ftn:2388-2403) writes this marker for bins whose
    # density per rad/s is <=1e-12. Bound the omitted energy explicitly.
    density[censored] = 0.0
    omitted_bound = float(np.trapezoid(censored.astype(float) * (2 * math.pi * 1e-12), frequency))
    integral = float(np.trapezoid(density, frequency))
    if integral <= 0:
        raise RuntimeError("Final variance spectrum is empty")
    return dict(resolved_variance_m2=integral,
                below_output_floor_bins=int(censored.sum()),
                below_output_floor_variance_bound_m2=omitted_bound,
                low_two_interval_fraction=float(np.trapezoid(density[:3], frequency[:3]) / integral),
                high_two_interval_fraction=float(np.trapezoid(density[-3:], frequency[-3:]) / integral),
                low_hz=float(frequency[0]), high_hz=float(frequency[-1]),
                peak_at_band_edge=bool(np.argmax(density) in (0, len(density) - 1)))


def breaking_fraction(path: Path, profile: np.ndarray) -> float:
    values = np.loadtxt(path, ndmin=2)
    if values.shape != (profile.shape[0] * profile.shape[1], 3) or not np.isfinite(values).all():
        raise RuntimeError("Incomplete breaking-fraction output")
    cube = values.reshape(profile.shape[0], profile.shape[1], 3)
    if not np.array_equal(cube[..., :2], profile[..., :2]):
        raise RuntimeError("Breaking-fraction output locations/times do not match wave output")
    fraction = cube[1:, :, 2]
    if np.any(fraction < 0) or np.any(fraction > 1):
        raise RuntimeError("Invalid depth-breaking fraction")
    return float(fraction.max())


def run_case(case: Case, executable: Path, directory: Path, *, timeout_s: float = 180) -> dict:
    """Run serially in an isolated directory, caching only identical verified inputs."""
    executable = executable.resolve(strict=True)
    input_text = case.input_text()
    depth_text = f"{case.depth_m:.12g} {case.depth_m:.12g}\n"
    identity = dict(case=asdict(case), executable_sha256=sha256(executable),
                    input_sha256=hashlib.sha256(input_text.encode()).hexdigest(),
                    depth_sha256=hashlib.sha256(depth_text.encode()).hexdigest())
    record_path = directory / "run.json"
    if directory.exists():
        if not record_path.exists():
            raise FileExistsError(f"Incomplete previous run in {directory}; inspect it before retrying")
        record = json.loads(record_path.read_text())
        if record["identity"] != identity:
            raise FileExistsError(f"Different inputs or executable already recorded in {directory}")
        if (sha256(directory / "INPUT") != identity["input_sha256"]
                or sha256(directory / "depth.bot") != identity["depth_sha256"]):
            raise RuntimeError(f"Cached SWAN input hash differs in {directory}")
        for name, digest in record["output_sha256"].items():
            if sha256(directory / name) != digest:
                raise RuntimeError(f"Cached output hash differs: {directory / name}")
        profile = read_profile(directory / "profile.tbl", case)
        record["maximum_breaking_fraction"] = breaking_fraction(directory / "breaking.tbl", profile)
        record["spectral_edges"] = read_spectral_edges(directory / "end.spc", case)
        record["analysis_sha256"] = sha256(Path(__file__))
        record_path.write_text(json.dumps(record, indent=2) + "\n")
        return record
    directory.mkdir(parents=True)
    (directory / "INPUT").write_text(input_text)
    (directory / "depth.bot").write_text(depth_text)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

    def limits():
        resource.setrlimit(resource.RLIMIT_CPU, (math.ceil(timeout_s), math.ceil(timeout_s)))
        resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))

    start = time.monotonic()
    with (directory / "console.log").open("w") as stream:
        process = subprocess.run([str(executable)], cwd=directory, env=environment,
                                 stdout=stream, stderr=subprocess.STDOUT,
                                 timeout=timeout_s, preexec_fn=limits)
    elapsed = time.monotonic() - start
    diagnostic = (directory / "PRINT").read_text() if (directory / "PRINT").exists() else ""
    if (process.returncode or not (directory / "norm_end").exists()
            or re.search(r"\*\*\s+(?:Error|Severe|Fatal)", diagnostic, re.I)):
        raise RuntimeError(f"SWAN failed in {directory}; inspect PRINT and console.log")
    profile = read_profile(directory / "profile.tbl", case)
    edges = read_spectral_edges(directory / "end.spc", case)
    record = dict(schema="terluna.climate.swan-run/1", identity=identity,
                  producer_sha256=sha256(Path(__file__)), analysis_sha256=sha256(Path(__file__)),
                  elapsed_wall_s=elapsed,
                  threads=1, address_space_limit_bytes=2 * 1024**3,
                  spectral_edges=edges,
                  maximum_breaking_fraction=breaking_fraction(directory / "breaking.tbl", profile),
                  warnings=[line.strip() for line in diagnostic.splitlines() if "Warning" in line],
                  output_sha256={name: sha256(directory / name)
                                 for name in ("profile.tbl", "breaking.tbl", "end.spc", "PRINT", "norm_end")})
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    return record


def wind_inputs() -> tuple[dict, dict]:
    path = ROOT / "climate/results/crm/ring_ring_equator.json"
    data = json.loads(path.read_text())
    if data.get("schema") != "terluna.climate.crm-ring/1":
        raise ValueError("Unsupported CM1 ring product schema")
    wind = {name: float(data["wind_10m_m_s"]["water"][name]) for name in ("median", "p90", "p99")}
    if any(not math.isfinite(value) or value <= 0 for value in wind.values()):
        raise ValueError("Invalid water wind quantiles")
    return wind, dict(path=str(path.relative_to(ROOT)), sha256=sha256(path),
                      schema=data["schema"], span_days=data["span_days"], snapshots=data["snapshots"])


def main_cases() -> list[Case]:
    wind, _ = wind_inputs()
    return [Case(f"{body}_{quantile}", gravity, value)
            for body, gravity in (("earth", STANDARD_GRAVITY), ("moon", MOON_SURFACE_GRAVITY))
            for quantile, value in wind.items()]
