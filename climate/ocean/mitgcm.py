"""MITgcm on the Moon's seas: one depth-averaged layer under the Moon's tide.

An experiment is a folder (runs/<name>/) holding MITgcm's namelists, the
bathymetry, the grid spacing, the tide's geopotential and the output. The
lunar settings are all run-time values: gravity, the reference sphere's
radius, the rotation period (the sidereal month) and the Coriolis map
2 Omega sin(latitude). The tide enters through the forcing package's tidal
geopotential, which MITgcm applies as -grad(phi) in the momentum equations,
so the file carries phi = -g F for an equilibrium height F.

The free surface is linear for the checks, as in geography's solver; the
strait study carries each column's full thickness (z-star). Bed drag is solved
implicitly, together with the surface pressure
(selectImplicitDrag = 2, compiled in through CPP_OPTIONS.h): MITgcm's default
explicit drag, stepped by Adams-Bashforth, holds only while r dt / H stays
below about 1, and shallow cells under strong friction pass that.

The grid is spherical-polar with any spacing in longitude and latitude, so a
fine region can sit inside a coarse sea. A run proceeds in chunks, each
ending in a permanent pickup; a stopped run resumes from its last pickup, and
only under the setup it started with (setup.sha256: the namelists, the fields
and the build).
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

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"


# ---------------------------------------------------------------- grids

def uniform_edges(start, stop, step):
    n = int(round((stop - start) / step))
    return start + step * np.arange(n + 1)


def stretched_edges(start, stop, fine, fine_step, coarse_step, ratio=1.1):
    """Cell edges from `start` to `stop`: `fine_step` inside fine=(a, b), growing by `ratio` per cell to
    `coarse_step` outside. Edges fall on the fine region's bounds exactly."""
    a, b = fine
    inner = uniform_edges(a, b, fine_step)
    def outward(origin, limit, sign):
        edges, step, x = [], fine_step, origin
        while sign * (limit - x) > 1e-9:
            step = min(step * ratio, coarse_step)
            x = x + sign * step
            if sign * (x - limit) > 0 or sign * (limit - x) < 0.5 * coarse_step:
                x = limit
            edges.append(x)
        return edges
    left = outward(a, start, -1)[::-1]
    right = outward(b, stop, +1)
    return np.array(left + list(inner) + right)


@dataclass
class Grid:
    """A spherical-polar grid: edges in degrees (longitudes east, latitudes south to north) and depth (m)."""
    xg: np.ndarray
    yg: np.ndarray
    depth: np.ndarray                  # (ny, nx), positive metres, 0 on land
    name: str = ""

    @property
    def nx(self):
        return len(self.xg) - 1

    @property
    def ny(self):
        return len(self.yg) - 1

    @property
    def xc(self):
        return 0.5 * (self.xg[1:] + self.xg[:-1])

    @property
    def yc(self):
        return 0.5 * (self.yg[1:] + self.yg[:-1])

    @property
    def area(self):
        """Cell areas on the Moon's sphere (m^2), (ny, nx)."""
        r = lunar.RADIUS
        band = r * r * (np.sin(np.radians(self.yg[1:])) - np.sin(np.radians(self.yg[:-1])))
        return band[:, None] * np.radians(np.diff(self.xg))[None, :]

    @property
    def wet(self):
        return self.depth > 0


def aggregate(terrain_box, mask, xg, yg, min_wet_share=0.5, min_depth=5.0):
    """Depth of each cell of a grid with edges (xg, yg) from a finer terrain: wet where at least
    `min_wet_share` of the fine cells whose centres fall inside it are water of `mask`, with the mean depth
    of those water cells; cells shallower than `min_depth` count as land."""
    lat, lon = terrain_box["lat"], terrain_box["lon"]
    d = np.where(mask, terrain_box["level"] - terrain_box["height"], 0.0)
    jj = np.searchsorted(yg, lat) - 1
    ii = np.searchsorted(xg, lon) - 1
    rows = (jj >= 0) & (jj < len(yg) - 1)
    cols = (ii >= 0) & (ii < len(xg) - 1)
    jr, ic = jj[rows], ii[cols]
    sub_d = d[np.ix_(rows, cols)]
    sub_w = (mask & (d > 0))[np.ix_(rows, cols)]
    ny, nx = len(yg) - 1, len(xg) - 1
    flat = (jr[:, None] * nx + ic[None, :]).ravel()
    total = np.bincount(flat, minlength=ny * nx).reshape(ny, nx)
    wet = np.bincount(flat, weights=sub_w.ravel().astype(float), minlength=ny * nx).reshape(ny, nx)
    deep = np.bincount(flat, weights=(sub_d * sub_w).ravel(), minlength=ny * nx).reshape(ny, nx)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(wet > 0, deep / wet, 0.0)
        share = np.where(total > 0, wet / total, 0.0)
    out = np.where((share >= min_wet_share) & (mean >= min_depth), mean, 0.0)
    return out


def keep_connected(depth, seeds_ji):
    """Zero every wet region of `depth` (joined through faces) that holds none of the seed cells."""
    from climate.ocean.terrain import components
    labels, _ = components(depth > 0)
    keep = {int(labels[j, i]) for j, i in seeds_ji if labels[j, i] > 0}
    return np.where(np.isin(labels, list(keep)), depth, 0.0)


# ---------------------------------------------------------------- binary files

def write_field(path, array, precision=64):
    """A big-endian MITgcm input file (the gfortran build reads big-endian)."""
    np.asarray(array, dtype=f">f{precision // 8}").tofile(path)


def read_meta(path):
    text = Path(path).read_text()
    dims = [int(x) for x in re.search(r"dimList\s*=\s*\[(.*?)\]", text, re.S).group(1).replace(",", " ").split()]
    prec = re.search(r"dataprec\s*=\s*\[\s*'(\w+)'", text).group(1)
    nrec = int(re.search(r"nrecords\s*=\s*\[\s*(\d+)", text).group(1))
    names = re.findall(r"'(\w+)\s*'", re.search(r"fldList\s*=\s*\{(.*?)\}", text, re.S).group(1)) if "fldList" in text else []
    shape = tuple(dims[0::3][::-1])        # dimList holds (global size, start, end) per dimension, x first
    return dict(shape=shape, dtype=">f4" if prec == "float32" else ">f8", nrecords=nrec, fields=names)


def read_mds(stem):
    """An MDS output (stem.meta and stem.data) as {field: array} (or one array when it holds no field names)."""
    meta = read_meta(f"{stem}.meta")
    data = np.fromfile(f"{stem}.data", dtype=meta["dtype"]).astype(float)
    data = data.reshape((meta["nrecords"],) + meta["shape"])
    if meta["fields"]:
        return {name: data[k] for k, name in enumerate(meta["fields"])}
    return data[0] if meta["nrecords"] == 1 else data


# ---------------------------------------------------------------- the experiment

def size_h(nx, ny, npx, npy, olx=3):
    if nx % npx or ny % npy:
        raise ValueError(f"{nx}x{ny} does not divide into {npx}x{npy} tiles")
    return f"""C SIZE.h for a one-layer Terluna sea: {nx} x {ny} cells on {npx} x {npy} processes
      INTEGER sNx, sNy, OLx, OLy, nSx, nSy, nPx, nPy, Nx, Ny, Nr
      PARAMETER (
     &           sNx = {nx // npx},
     &           sNy = {ny // npy},
     &           OLx = {olx},
     &           OLy = {olx},
     &           nSx = 1,
     &           nSy = 1,
     &           nPx = {npx},
     &           nPy = {npy},
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  = 1)
      INTEGER MAX_OLX, MAX_OLY
      PARAMETER ( MAX_OLX = OLx, MAX_OLY = OLy )
"""


PACKAGES_CONF = "mom_common\nmom_fluxform\nmom_vecinv\ngeneric_advdiff\ncal\nexf\ndiagnostics\nmdsio\nrw\nmonitor\n"


def exf_options(source: Path) -> str:
    """The forcing package's options with the tidal geopotential and input interpolation switched on."""
    text = (source / "pkg" / "exf" / "EXF_OPTIONS.h").read_text()
    for flag in ("EXF_ALLOW_TIDES", "USE_EXF_INTERPOLATION"):
        text = re.sub(rf"#undef\s+{flag}\b", f"#define {flag}", text)
        if f"#define {flag}" not in text:
            text = text.replace("#endif /* ALLOW_EXF */", f"#define {flag}\n#endif /* ALLOW_EXF */")
    return text


def cpp_options(source: Path, nonlinear_surface=False) -> str:
    """The model's options with bed drag solvable together with the surface pressure (selectImplicitDrag = 2) and,
    when asked, the nonlinear free surface (z-star)."""
    text = (source / "model" / "inc" / "CPP_OPTIONS.h").read_text()
    flags = ["ALLOW_SOLVE4_PS_AND_DRAG"] + (["NONLIN_FRSURF"] if nonlinear_surface else [])
    for flag in flags:
        text, n = re.subn(rf"#undef\s+{flag}\b", f"#define {flag}", text)
        if n != 1:
            raise RuntimeError(f"CPP_OPTIONS.h lacks the {flag} switch this setup expects")
    if (not re.search(r"^#define\s+INCLUDE_IMPLVERTADV_CODE\b", text, re.M)
            or re.search(r"^#define\s+DISABLE_RSTAR_CODE\b", text, re.M)):
        raise RuntimeError("CPP_OPTIONS.h no longer has the implicit-drag and z-star code this setup expects")
    return text


def namelist(group, values):
    lines = [f" &{group}"]
    for key, value in values.items():
        if isinstance(value, bool):
            value = ".TRUE." if value else ".FALSE."
        elif isinstance(value, str) and not value.startswith(("'", '"')) and "*" not in value:
            value = f"'{value}'"
        elif isinstance(value, (list, tuple, np.ndarray)):
            items = [f"{float(v):.10g}" for v in value]
            value = ",\n   ".join(", ".join(items[k:k + 8]) for k in range(0, len(items), 8))
        elif isinstance(value, float):
            value = f"{value:.12g}"
        lines.append(f" {key}={value},")
    lines.append(" &")
    return "\n".join(lines) + "\n"


@dataclass
class Experiment:
    """One MITgcm configuration on a Grid, in runs/<name>/."""
    name: str
    grid: Grid
    deltat: float                       # s
    nonlinear: bool = False             # momentum advection and quadratic bed drag
    nonlinear_surface: bool = False     # z-star free surface: the column's thickness follows the level
    drag: float = lunar.DRAG            # quadratic coefficient when nonlinear
    linear_friction: float = 0.0        # m/s, r in r u / H, when linear
    rotation: bool = True
    tide: dict | None = None            # forcing grid: dict(lon0, lat0, step, nlon, nlat, period_s, file)
    output_s: float = 10800.0           # snapshot interval
    output_fields: tuple = ("ETAN", "UVEL", "VVEL")
    output_from_s: float = 0.0          # chunks ending before this write no snapshots
    implicit: float = 0.5               # implicSurfPress = implicDiv2DFlow (0.5 neutral, 1 damping)
    viscosity_smag: float = 0.0         # Smagorinsky coefficient for the nonlinear runs
    tiles: tuple = (1, 1)
    cg2d_residual: float = 1e-10        # the barotropic solver's target (relative to the right-hand side's largest)
    initial: dict = field(default_factory=dict)   # optional fields: eta, u, v (arrays on the grid)

    @property
    def folder(self) -> Path:
        return RUNS / "mitgcm" / self.name

    def code(self) -> dict:
        source = builds.checkout("mitgcm")
        return {"SIZE.h": size_h(self.grid.nx, self.grid.ny, *self.tiles), "packages.conf": PACKAGES_CONF,
                "CPP_OPTIONS.h": cpp_options(source, self.nonlinear_surface), "EXF_OPTIONS.h": exf_options(source)}

    def binary(self, jobs=3) -> Path:
        return builds.mitgcm(self.code(), mpi=self.tiles != (1, 1), jobs=jobs)

    def inputs(self) -> tuple[dict, dict]:
        """The run's namelists (file name to text) and binary fields (file name to array)."""
        g = self.grid
        fields = {"bathy.bin": -g.depth}
        files = dict(bathyFile="bathy.bin")
        for key, name in (("eta", "pSurfInitFile"), ("u", "uVelInitFile"), ("v", "vVelInitFile")):
            if key in self.initial:
                fields[f"init_{key}.bin"] = self.initial[key]          # eta in metres, flows in m/s
                files[name] = f"init_{key}.bin"
        hmax = float(g.depth.max()) * 1.01 + 1.0
        parm01 = dict(gravity=lunar.GRAVITY, gBaro=lunar.GRAVITY, rhoConst=lunar.WATER_DENSITY, rhoNil=lunar.WATER_DENSITY,
                      rotationPeriod=lunar.ROTATION_PERIOD_S, selectCoriMap=2 if self.rotation else 0,
                      f0=0.0, beta=0.0,
                      implicitFreeSurface=True, rigidLid=False, nonlinFreeSurf=4 if self.nonlinear_surface else 0,
                      implicSurfPress=self.implicit, implicDiv2DFlow=self.implicit, exactConserv=True,
                      momAdvection=self.nonlinear,
                      momViscosity=True, viscC2smag=self.viscosity_smag, useFullLeith=False,
                      bottomDragQuadratic=self.drag if self.nonlinear else 0.0,
                      bottomDragLinear=0.0 if self.nonlinear else self.linear_friction,
                      selectImplicitDrag=2,        # bed drag solved with the surface pressure
                      no_slip_sides=False, no_slip_bottom=False,
                      tempStepping=False, saltStepping=False, momTidalForcing=self.tide is not None,
                      eosType="LINEAR", tAlpha=0.0, sBeta=0.0,
                      hFacMin=1e-6, hFacMinDr=1.0,
                      readBinaryPrec=64, writeBinaryPrec=32, useSingleCpuIO=True, globalFiles=True)
        if not self.rotation:
            parm01.pop("rotationPeriod")
        if self.viscosity_smag > 0:     # Smagorinsky held below L^2 / (8 dt), the explicit limit under Adams-Bashforth
            parm01.update(viscAhGridMax=0.5)
        if self.nonlinear_surface:      # MITgcm cannot dry a cell: it stops when a column thins below hFacInf of its depth
            parm01.update(select_rStar=2, hFacInf=0.05, hFacSup=3.0)
        text = namelist("PARM01", parm01)
        text += namelist("PARM02", dict(cg2dTargetResidual=self.cg2d_residual, cg2dMaxIters=2000, cg2dUseMinResSol=1))
        text += "PARM03_PLACEHOLDER"
        text += namelist("PARM04", dict(usingSphericalPolarGrid=True, rSphere=lunar.RADIUS,
                                        xgOrigin=float(g.xg[0]), ygOrigin=float(g.yg[0]),
                                        delX=np.diff(g.xg), delY=np.diff(g.yg), delR=[hmax]))
        text += namelist("PARM05", files)
        texts = {"data.template": text,
                 "data.pkg": namelist("PACKAGES", dict(useEXF=True, useCAL=True, useDiagnostics=True)),
                 "eedata": " &EEPARMS\n nTx=1,\n nTy=1,\n &\n",
                 "data.cal": namelist("CAL_NML", dict(TheCalendar="gregorian", startDate_1=20260101, startDate_2=0))}
        exf = namelist("EXF_NML_01", dict(useExfCheckRange=False, exf_iprec=64, useExfYearlyFields=False))
        if self.tide:
            t = self.tide
            exf += namelist("EXF_NML_02", dict(tidePotFile=t["file"], tidePotstartdate1=20260101, tidePotstartdate2=0,
                                               tidePotPeriod=t["period_s"]))   # record n at (n - 1) periods from model time zero
            exf += namelist("EXF_NML_03", dict(exf_inscal_tidePot=-lunar.GRAVITY))
            exf += namelist("EXF_NML_04", dict(tidePot_lon0=t["lon0"], tidePot_lon_inc=t["step"],
                                               tidePot_lat0=t["lat0"], tidePot_lat_inc=f"{t['nlat'] - 1}*{t['step']}",
                                               tidePot_nlon=t["nlon"], tidePot_nlat=t["nlat"], tidePot_interpMethod=1))
        else:
            exf += namelist("EXF_NML_02", {}) + namelist("EXF_NML_03", {}) + namelist("EXF_NML_04", {})
        texts["data.exf"] = exf
        texts["experiment.json"] = json.dumps(dict(
            name=self.name, nx=g.nx, ny=g.ny, wet_cells=int(g.wet.sum()), deltat=self.deltat, nonlinear=self.nonlinear,
            drag=self.drag, linear_friction=self.linear_friction, rotation=self.rotation, implicit=self.implicit,
            nonlinear_surface=self.nonlinear_surface,
            viscosity_smag=self.viscosity_smag, tiles=list(self.tiles), output_s=self.output_s, tide=self.tide,
            implicit_drag=True, gravity=lunar.GRAVITY, radius=lunar.RADIUS, rotation_period_s=lunar.ROTATION_PERIOD_S,
            water_density=lunar.WATER_DENSITY), indent=1) + "\n"
        return texts, fields

    def setup_digest(self, build: str) -> str:
        texts, fields = self.inputs()
        h = hashlib.sha256(json.dumps(dict(texts=texts, build=build), sort_keys=True).encode())
        for name in sorted(fields):
            h.update(name.encode() + np.ascontiguousarray(fields[name], dtype=">f8").tobytes())
        return h.hexdigest()

    def write(self, fresh=False):
        """Write the run's inputs into its folder; `fresh` first clears a folder left by an earlier run."""
        f = self.folder
        if fresh and f.exists():
            if not f.resolve().is_relative_to((RUNS / "mitgcm").resolve()):
                raise RuntimeError(f"{f} is outside the runs folder")
            shutil.rmtree(f)
        f.mkdir(parents=True, exist_ok=True)
        texts, fields = self.inputs()
        for name, values in fields.items():
            write_field(f / name, values)
        for name, text in texts.items():
            (f / name).write_text(text)
        g = self.grid
        np.savez_compressed(f / "grid.npz", xg=g.xg, yg=g.yg, depth=g.depth)

    def write_diagnostics(self, active=True):
        fields = ",".join(f"'{name:<8s}'" for name in self.output_fields)
        frequency = -self.output_s if active else -1.0e12
        (self.folder / "data.diagnostics").write_text(
            " &DIAGNOSTICS_LIST\n"
            f"  fields(1:{len(self.output_fields)},1) = {fields},\n"
            "  fileName(1) = 'state',\n"
            f"  frequency(1) = {frequency:.1f},\n"
            "  timePhase(1) = 0.,\n"
            " &\n &DIAG_STATIS_PARMS\n &\n")

    def pickups(self):
        return sorted(int(p.name.split(".")[1]) for p in self.folder.glob("pickup.*.data") if p.name.split(".")[1].isdigit())

    def run(self, until_s: float, chunk_s: float, *, threads=None, log_name="run.log", after_chunk=None):
        """Advance to model time `until_s` in chunks of `chunk_s`, each from the latest pickup."""
        exe = self.binary()
        f = self.folder
        steps_per_chunk = int(round(chunk_s / self.deltat))
        total = int(round(until_s / self.deltat))
        ranks = self.tiles[0] * self.tiles[1]
        digest = self.setup_digest(exe.parent.name)
        stamp = f / "setup.sha256"
        if not self.pickups():
            stamp.write_text(digest + "\n")
        elif not stamp.exists() or stamp.read_text().split()[0] != digest:
            raise RuntimeError(f"{f} holds pickups from another setup; write(fresh=True) starts the run again")
        while True:
            done = max([0] + self.pickups())
            if done >= total:
                return done
            n = min(steps_per_chunk, total - done)
            # permanent pickups fall on multiples of pChkptFreq from time zero: the chunk's end time makes its end
            # the only one (an end-of-run pickup alone would take the rolling checkpoint's name)
            parm03 = dict(nIter0=done, nTimeSteps=n, deltaT=self.deltat, pChkptFreq=(done + n) * self.deltat,
                          chkptFreq=0.0, dumpFreq=0.0, monitorFreq=max(self.deltat, n * self.deltat / 4),
                          monitorSelect=2, pickupStrictlyMatch=False, writePickupAtEnd=True)
            text = (f / "data.template").read_text().replace("PARM03_PLACEHOLDER", namelist("PARM03", parm03))
            (f / "data").write_text(text)
            self.write_diagnostics(active=(done + n) * self.deltat > self.output_from_s + 1e-6)
            command = ([str(exe)] if ranks == 1 else ["mpirun", "--bind-to", "none", "--mca", "mpi_yield_when_idle", "1",
                                                      "-np", str(ranks), str(exe)])
            env = dict(os.environ, OMP_NUM_THREADS="1")
            with open(f / log_name, "a") as out:
                out.write(f"\n== chunk from iteration {done}, {n} steps\n")
                out.flush()
                result = subprocess.run(["nice", "-n", "10"] + command, cwd=f, stdout=out, stderr=subprocess.STDOUT, env=env)
            if result.returncode != 0 or (done + n) not in self.pickups():
                raise RuntimeError(f"MITgcm stopped in {f} (see {log_name} and STDOUT.0000)")
            if after_chunk:
                after_chunk(self, done, done + n)

    def snapshots(self, first=0, last=None):
        """Output iterations available (state.<iteration>.data), in order."""
        its = sorted(int(p.name.split(".")[1]) for p in self.folder.glob("state.*.data"))
        return [i for i in its if i >= first and (last is None or i <= last)]

    def state(self, iteration):
        return read_mds(self.folder / f"state.{iteration:010d}")


def pad_even(grid: Grid) -> Grid:
    """Add a land column and row where needed so the grid splits into 2 x 2 tiles."""
    xg, yg, depth = grid.xg, grid.yg, grid.depth
    if grid.nx % 2:
        xg = np.r_[xg, 2 * xg[-1] - xg[-2]]
        depth = np.c_[depth, np.zeros(depth.shape[0])]
    if grid.ny % 2:
        yg = np.r_[yg, 2 * yg[-1] - yg[-2]]
        depth = np.r_[depth, np.zeros((1, depth.shape[1]))]
    return Grid(xg=xg, yg=yg, depth=depth, name=grid.name)


def tide_forcing_grid(grid: Grid, step=1.0, margin=3.0):
    """A regular forcing grid covering `grid` with a margin (degrees)."""
    lon0 = math.floor(grid.xg[0] - margin)
    lat0 = math.floor(grid.yg[0] - margin)
    nlon = int(math.ceil((grid.xg[-1] + margin - lon0) / step)) + 1
    nlat = int(math.ceil((grid.yg[-1] + margin - lat0) / step)) + 1
    return dict(lon0=float(lon0), lat0=float(lat0), step=float(step), nlon=nlon, nlat=nlat)
