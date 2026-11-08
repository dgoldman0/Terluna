"""SCHISM on the Moon's seas: the lunar build, depth-averaged (two vertical levels).

An experiment is a folder (runs/schism/<name>/) with the mesh (hgrid.gr3 and
hgrid.ll in degrees, ics = 2), a two-level vertical grid (SCHISM's 2-D mode),
param.nml from the pinned source's sample with this study's values, an empty
bctides.in (no Earth constituents, no open boundaries), the bed drag
(drag.gr3, quadratic Cd) and the Moon's tide (lunar_tide.th and
lunar_tide_basis.in, ilunartide = 1). Outputs are elevation and depth-averaged
velocity at the nodes (outputs/out2d_*.nc), on one writer rank beside the
compute ranks. A run proceeds in chunks, each ending in a hotstart that the
next chunk starts from, and resumes only under the setup it started with
(setup.sha256: the settings, mesh, initial state, tide files and build).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from climate.ocean import builds, lunar
from climate.ocean.mesh import Mesh

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
VGRID_2D = "2 !ivcor\n2 1 1.e6 !nvrt, kz, h_s\nZ levels\n1 -1.e6\nS levels\n40. 1. 1.e-4 !h_c, theta_b, theta_f\n1 -1.\n2 0.\n"


def write_gr3(path, mesh: Mesh, values, title, boundaries=True):
    """A .gr3 file: nodes with one value each, the elements, and (for hgrid) the land boundaries."""
    with open(path, "w") as out:
        out.write(f"{title}\n{mesh.elements} {mesh.nodes}\n")
        nodes = np.column_stack([np.arange(1, mesh.nodes + 1), mesh.lon, mesh.lat, values])
        np.savetxt(out, nodes, fmt=["%d", "%.8f", "%.8f", "%.6f"])
        elements = np.column_stack([np.arange(1, mesh.elements + 1), np.full(mesh.elements, 3), mesh.triangles + 1])
        np.savetxt(out, elements, fmt="%d")
        if boundaries:
            loops = mesh.boundary_loops()
            out.write("0 = Number of open boundaries\n0 = Total number of open boundary nodes\n")
            out.write(f"{len(loops)} = Number of land boundaries\n{sum(len(l) for l in loops)} = Total number of land boundary nodes\n")
            for k, loop in enumerate(loops):
                out.write(f"{len(loop)} {0 if k == 0 else 1} = Number of nodes for land boundary {k + 1}\n")
                for n in loop:
                    out.write(f"{n + 1}\n")


def chunk_steps(start_s, span_s, dt, output_s, file_s):
    """A chunk's steps in each output file (ihfskip) and between hotstarts (nhot_write).

    A full chunk takes file_s's files and one hotstart at its end. A shorter last chunk takes files of the largest
    number of steps dividing both its length and its first step, so that a continued run's files line up with its
    start, and, since SCHISM writes a hotstart whenever the run's own step count is a multiple of nhot_write, a
    hotstart interval equal to its last step."""
    steps, per_file = int(round(span_s / dt)), int(round(file_s / dt))
    if steps % per_file == 0:
        return per_file, steps
    first = int(round(start_s / dt))
    per_file = math.gcd(steps, first)
    if per_file % int(round(output_s / dt)):
        raise ValueError("a chunk must hold a whole number of output files (nhot_write a multiple of ihfskip)")
    return per_file, first + steps


def param_nml(values: dict, opt_extra: dict) -> str:
    """The pinned source's sample param.nml with `values` set (key = value lines) and extra &OPT entries."""
    text = (builds.checkout("schism") / "sample_inputs" / "param.nml").read_text()
    for key, value in values.items():
        if isinstance(value, bool):
            value = 1 if value else 0
        pattern = re.compile(rf"^(\s*){re.escape(key)}\s*=\s*[^!\n]*", re.M)
        if not pattern.search(text):
            raise KeyError(f"{key} is not in the sample param.nml")
        text = pattern.sub(lambda m: f"{m.group(1)}{key} = {value} ", text, count=1)
    extra = "".join(f"  {k} = {v}\n" for k, v in opt_extra.items())
    text = re.sub(r"(&OPT\s*\n)", lambda m: m.group(1) + extra, text, count=1)
    # v5.13.0 declares veg_lai and veg_cw INTEGER while its sample gives them reals (1. and 1.5), which gfortran
    # cannot read; vegetation is off here, so they keep the code's defaults. Comments are dropped as well.
    lines = [re.sub(r"\s*!.*$", "", line) for line in text.splitlines()]
    lines = [line for line in lines if line.strip() and not re.match(r"\s*veg_(lai|cw)\s*=", line)]
    return "\n".join(lines) + "\n"


@dataclass
class Experiment:
    name: str
    mesh: Mesh
    dt: float                          # s
    drag: float = lunar.DRAG           # quadratic Cd (0 for none)
    advection: bool = True
    rotation: bool = True
    output_s: float = 10800.0
    file_s: float = 86400.0 * 10       # one output file per this span
    ramp_days: float = 0.0             # SCHISM's own tanh ramp of the tide (dramp); 0 when the series carries it
    initial_eta: np.ndarray | None = None
    compute_ranks: int = 3
    tide: bool = False                 # set by write() when the Moon's tide files are written
    filter_mode: int = 1               # the Shapiro filter with advection (ishapiro): 1 a fixed strength each step,
                                       # 2 Smagorinsky-like, 0.5 tanh(coefficient dt mu) with mu the local deformation
    filter_strength: float = 0.5       # shapiro0 (mode 1) or the Smagorinsky-like coefficient in shapiro.gr3 (mode 2)
    btrack_max_s: float = 30.0         # the backtracking's longest and shortest sub-steps (dtb_max, dtb_min): each step
    btrack_min_s: float = 10.0         # takes at least dt / dtb_max of them, more where the flow's gradients ask

    @property
    def folder(self) -> Path:
        return RUNS / "schism" / self.name

    def write(self, tide_seconds=None, tide_coefs=None, attraction=True, fresh=False):
        """Write the run's inputs into its folder; `fresh` first clears a folder left by an earlier run."""
        f = self.folder
        if fresh and f.exists():
            if not f.resolve().is_relative_to((RUNS / "schism").resolve()):
                raise RuntimeError(f"{f} is outside the runs folder")
            shutil.rmtree(f)
        (f / "outputs").mkdir(parents=True, exist_ok=True)
        m = self.mesh
        write_gr3(f / "hgrid.gr3", m, m.depth, "hgrid (lon, lat, depth)")
        shutil.copy(f / "hgrid.gr3", f / "hgrid.ll")
        write_gr3(f / "drag.gr3", m, np.full(m.nodes, self.drag), "quadratic bed drag Cd", boundaries=False)
        (f / "vgrid.in").write_text(VGRID_2D)
        if not self.advection:                  # nadv = 0 reads momentum advection on (1) or off (0) per node
            write_gr3(f / "adv.gr3", m, np.zeros(m.nodes), "momentum advection off", boundaries=False)
        (f / "bctides.in").write_text("01/01/2026 00:00:00 TT (model time zero)\n0 0. !ntip, tip_dp\n0 !nbfr\n0 !nope\n")
        if self.initial_eta is not None:
            write_gr3(f / "elev.ic", m, self.initial_eta, "initial elevation", boundaries=False)
        if self.advection and self.filter_mode in (-1, 2):
            write_gr3(f / "shapiro.gr3", m, np.full(m.nodes, self.filter_strength), "Shapiro filter coefficient",
                      boundaries=False)
        if tide_coefs is not None:
            from climate.ocean import forcing
            forcing.write_schism(f, m.lat, m.lon, tide_seconds, tide_coefs, attraction=attraction)
            self.tide = True
        (f / "experiment.json").write_text(json.dumps(dict(
            name=self.name, nodes=m.nodes, elements=m.elements, dt=self.dt, drag=self.drag, advection=self.advection,
            filter_mode=self.filter_mode, filter_strength=self.filter_strength, btrack_max_s=self.btrack_max_s,
            btrack_min_s=self.btrack_min_s,
            rotation=self.rotation, output_s=self.output_s, ramp_days=self.ramp_days, compute_ranks=self.compute_ranks,
            tide=self.tide or tide_coefs is not None, gravity=lunar.GRAVITY, radius=lunar.RADIUS, rotation_period_s=lunar.ROTATION_PERIOD_S),
            indent=1) + "\n")
        np.savez_compressed(f / "mesh.npz", lon=m.lon, lat=m.lat, depth=m.depth, triangles=m.triangles)

    def _param(self, start_s, span_s, hot):
        nspool = int(round(self.output_s / self.dt))
        ihfskip, nhot_write = chunk_steps(start_s, span_s, self.dt, self.output_s, self.file_s)
        values = dict(ipre=0, ibc=1, ibtp=0, rnday=(start_s + span_s) / 86400.0, dt=self.dt, nspool=nspool,
                      ihfskip=ihfskip, ics=2, ihot=2 if hot else 0, dramp=self.ramp_days,
                      nchi=0, ncor=1 if self.rotation else 0, coricoef=0.0, nadv=1 if self.advection else 0,
                      ic_elev=1 if (self.initial_eta is not None and not hot) else 0, h0=0.01, thetai=0.6,
                      itur=0, ihorcon=0, ishapiro=self.filter_mode if self.advection else 0, indvel=0, inunfl=0, nws=0,
                      ihconsv=0, isconsv=0, iloadtide=0, nhot=1, nhot_write=nhot_write,
                      rearth_pole=lunar.RADIUS, rearth_eq=lunar.RADIUS,
                      dtb_max=self.btrack_max_s, dtb_min=self.btrack_min_s,
                      shapiro0=self.filter_strength if self.filter_mode == 1 else 0.5,
                      start_year=2026, start_month=1, start_day=1, start_hour=0, utc_start=0)
        if values["nhot_write"] % ihfskip:
            raise ValueError("a chunk must hold a whole number of output files (nhot_write a multiple of ihfskip)")
        text = param_nml(values, dict(ilunartide=1 if (self.folder / "lunar_tide.th").exists() else 0))
        # outputs: elevation and depth-averaged velocity; no 3-D velocity
        text = re.sub(r"^(\s*)iof_hydro\(26\)\s*=\s*\d", r"\g<1>iof_hydro(26) = 0", text, flags=re.M)
        # the code switches on 3-D layer heights by default (iof_hydro(25)); a 2-D run has no use for them
        text = re.sub(r"(&SCHOUT\s*\n)", lambda m: m.group(1) + "  iof_hydro(25) = 0\n", text, count=1)
        text = re.sub(r"^(\s*)iof_hydro\(16\)\s*=\s*\d", r"\g<1>iof_hydro(16) = 1", text, flags=re.M)
        return text

    def chunks_done(self):
        p = self.folder / "chunks.json"
        return json.loads(p.read_text()) if p.exists() else []

    def setup_digest(self, build: str) -> str:
        m = self.mesh
        h = hashlib.sha256(json.dumps(dict(
            name=self.name, dt=self.dt, drag=self.drag, advection=self.advection, rotation=self.rotation,
            output_s=self.output_s, file_s=self.file_s, ramp_days=self.ramp_days, compute_ranks=self.compute_ranks,
            param=self._param(0.0, self.file_s, False), build=build), sort_keys=True).encode())
        for values in (m.lon, m.lat, m.depth, m.triangles, self.initial_eta):
            if values is not None:
                h.update(np.ascontiguousarray(values, dtype=float).tobytes())
        for name in ("lunar_tide.th", "lunar_tide_basis.in"):
            if (self.folder / name).exists():
                h.update(name.encode() + (self.folder / name).read_bytes())
        return h.hexdigest()

    def run(self, until_s, chunk_s, *, after_chunk=None):
        exe = builds.schism()
        combine = exe.parent / "combine_hotstart7"
        f = self.folder
        digest = self.setup_digest(exe.parent.name)
        stamp = f / "setup.sha256"
        if not self.chunks_done():
            stamp.write_text(digest + "\n")
        elif not stamp.exists() or stamp.read_text().split()[0] != digest:
            raise RuntimeError(f"{f} holds chunks from another setup; write(fresh=True) starts the run again")
        while True:
            done = self.chunks_done()
            start = done[-1]["end_s"] if done else 0.0
            if start >= until_s - 1e-6:
                return start
            span = min(chunk_s, until_s - start)
            hot = start > 0
            if hot:
                shutil.copy(f / "hotstart_next.nc", f / "hotstart.nc")
            (f / "param.nml").write_text(self._param(start, span, hot))
            out = f / f"outputs_{len(done):03d}"
            if (f / "outputs").exists():
                shutil.rmtree(f / "outputs")
            (f / "outputs").mkdir()
            if hot:                             # a continued run appends to the previous chunk's running logs
                for log_file in (f / done[-1]["folder"]).glob("*.out"):
                    shutil.copy(log_file, f / "outputs" / log_file.name)
            env = dict(os.environ, OMP_NUM_THREADS="1")
            with open(f / "run.log", "a") as log:
                log.write(f"\n== chunk {len(done)} from {start / 86400:.3f} d for {span / 86400:.3f} d\n")
                log.flush()
                result = subprocess.run(["nice", "-n", "10", "mpirun", "--bind-to", "none", "--mca", "mpi_yield_when_idle", "1",
                                         "-np", str(self.compute_ranks + 1),
                                         str(exe), "1"], cwd=f, stdout=log, stderr=subprocess.STDOUT, env=env)
            steps = int(round((start + span) / self.dt))
            hotfiles = sorted((f / "outputs").glob(f"hotstart_*_{steps}.nc"))
            if result.returncode != 0 or not hotfiles:
                raise RuntimeError(f"SCHISM stopped in {f} (see run.log and outputs/mirror.out)")
            subprocess.run([str(combine), "-i", str(steps)], cwd=f / "outputs", check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            shutil.move(f / "outputs" / f"hotstart_it={steps}.nc", f / "hotstart_next.nc")
            for p in hotfiles:
                p.unlink()
            (f / "outputs").rename(out)
            done.append(dict(start_s=start, end_s=start + span, folder=out.name))
            (f / "chunks.json").write_text(json.dumps(done, indent=1) + "\n")
            if after_chunk:
                after_chunk(self, out, start, start + span)

    def series(self, variables=("elevation",), chunks=None):
        """Model times (s) and node fields of the chosen variables from the saved chunks' out2d files."""
        import netCDF4
        times, fields = [], {v: [] for v in variables}
        for chunk in self.chunks_done() if chunks is None else chunks:
            for path in sorted((self.folder / chunk["folder"]).glob("out2d_*.nc"), key=lambda p: int(p.stem.split("_")[1])):
                with netCDF4.Dataset(path) as nc:
                    t = np.array(nc["time"][:], dtype=float)        # seconds from model time zero
                    times.append(t)
                    for v in variables:
                        fields[v].append(np.array(nc[v][:], dtype=float))
        return np.concatenate(times), {v: np.concatenate(x) for v, x in fields.items()}
