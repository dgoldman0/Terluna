"""MITgcm on one of the Moon's seas with temperature and salt: its circulation under wind, heat, evaporation, rain and
rivers.

The tide runs (mitgcm.py) carry one depth-averaged layer of water at one density. This module adds what a sea's own
circulation needs: levels in the vertical, temperature and salinity with seawater's equation of state, a turbulence
closure, and the surface's fluxes read from files. It reuses the tide runs' grid, files, builds and chunked runs with
pickups.

The grid. A sea that reaches high latitudes is run on a spherical grid whose own equator passes through the sea's
centre (`Frame`), so that its cells keep one size; the Coriolis parameter then comes from files holding
2 Omega sin(true latitude) at the cells' centres and corners (selectCoriMap = 3), and winds are turned into the
grid's own east and north before they are written.

The lunar settings. Gravity, radius and rotation are run-time values, as in the tide runs. The equation of state
(Jackett and McDougall 1995) takes its pressure as rho g z with the Moon's gravity. The closure (pkg/ggl90, Gaspar et
al. 1990) builds its buoyancy frequency and its mixing lengths from that gravity; its constants are dimensionless
and were fitted to Earth's ocean, and it handles convection itself. Bed drag is quadratic with the tide runs'
coefficient.

The forcing. The forcing package reads wind stress, the net heat flux, the sunlight within it, the net fresh-water
flux and the rivers' inflow on the model's own grid; its bulk formulas, which hold Earth's constants, are compiled
out. Sunlight heats below the top layer by the model's two-band profile for clear ocean water (Jerlov's type IA,
which the model's source fixes), an Earth water's clarity taken as it stands. The free surface
carries the water the rivers and the rain add and the evaporation takes (z-star with real fresh-water flux), so a
river dilutes the cell it enters to its own strength and the sea's salt is conserved. The sea surface can be drawn
towards a prescribed temperature with a stated time scale.

A patch of fine cells. The grid's spacing may vary (mitgcm.stretched_edges), so a shore can carry fine cells inside
a coarse sea. The lateral viscosity then follows each cell's own size (viscAhGrid), and a run can start from fields
carried over from a coarser run in place of one temperature and salinity.

The old water. A dye starts at 1 everywhere and is drawn to 0 in the top layer within a day (pkg/ptracers with
pkg/rbcs), and is otherwise carried and mixed like salt. What is left of it at a depth is the share of that water
that has not been at the surface since the run began, so its fall measures how fast the deeps are renewed.

Other tracers (`Tracer`). A run may carry several tracers in place of the dye, each drawn towards its own value in
the top layer and given its own source below it. The relaxation package carries the source as a pull towards a
target far beyond the tracer's range (BEYOND), which adds the source to within the tracer's value over BEYOND.

The Moon's own interior mixing (`idemix`). MITgcm's IDEMIX (Olbers and Eden 2013, in pkg/ggl90 and compiled in only
when asked) follows the energy of the internal waves: what the tide puts in at the floor and the wind at the surface,
both given as maps in W/m2, propagates and dissipates by the model's own buoyancy frequency and Coriolis parameter,
and the dissipation feeds the closure's turbulence. Its constants are fitted to Earth's internal waves; it takes a
fifth of the wind's inertial power into the waves and divides both inputs by 1024 kg/m3. With it the run's background
diffusivity can fall to a molecular value. The floor's own heat (`geothermal`, W/m2) enters the deepest cell of each
column.

Restarts (`restart_from`). A run can begin from another run's pickup on the same grid and levels: its velocities,
temperature, salt, surface and turbulence carry over, and its tracers start from files at that iteration.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from climate.ocean import builds, lunar
from climate.ocean.mitgcm import Experiment, Grid, cpp_options, diagnostics_size, namelist

PACKAGES_CONF = ("mom_common\nmom_fluxform\nmom_vecinv\ngeneric_advdiff\ncal\nexf\nggl90\nptracers\nrbcs\ndiagnostics\n"
                 "mdsio\nrw\nmonitor\n")
FORCING = ("ustress", "vstress", "hflux", "swflux", "sflux", "climsst")   # fields that vary through the forcing's cycle
BEYOND = 1.0e6                  # a source's target: the source comes out low by the tracer's value over this


@dataclass
class Tracer:
    """A passive tracer carried and mixed like salt, drawn towards `top` in the top layer with the time scale
    `top_s`, and given `source` (tracer units a second, on every level below the top; positive adds)."""
    name: str
    start: np.ndarray | float = 1.0     # (levels, ny, nx), or one value everywhere
    top: np.ndarray | float = 0.0       # the top layer's value, (ny, nx) or one value
    top_s: float = 86400.0
    source: np.ndarray | None = None    # (levels, ny, nx)

    def relaxation(self, wet3d):
        """The relaxation package's mask and target, (levels, ny, nx) each: the top layer drawn to `top` at the rate
        1/top_s, and below it a pull at |source| top_s / BEYOND towards plus or minus BEYOND, which adds `source`."""
        mask = np.zeros(wet3d.shape)
        target = np.zeros(wet3d.shape)
        mask[0] = wet3d[0]
        target[0] = np.broadcast_to(self.top, wet3d.shape[1:])
        if self.source is not None:
            s = np.where(wet3d[1:], self.source[1:], 0.0)
            mask[1:] = np.abs(s) * self.top_s / BEYOND
            target[1:] = np.where(s < 0, -BEYOND, BEYOND)
        return mask, target


def options_with(text: str, flag: str, where: str) -> str:
    """An options file's text with `flag` defined in place of its #undef."""
    text, n = re.subn(rf"^#undef\s+{flag}\b", f"#define {flag}", text, flags=re.M)
    if n != 1:
        raise RuntimeError(f"{where} lacks the {flag} switch this setup expects")
    return text


# ---------------------------------------------------------------- the sea's own frame

def unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)


def latlon(v):
    return np.degrees(np.arcsin(np.clip(v[..., 2], -1.0, 1.0))), np.degrees(np.arctan2(v[..., 1], v[..., 0]))


def tangents(lat_deg, lon_deg):
    """Unit vectors to the east and to the north at places on the sphere."""
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    east = np.stack((-np.sin(lo), np.cos(lo), np.zeros_like(lo)), axis=-1)
    north = np.stack((-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)), axis=-1)
    return east, north


@dataclass(frozen=True)
class Frame:
    """A spherical frame whose equator and prime meridian cross at a chosen place, with its north towards true north
    there. `matrix` turns true unit vectors into the frame's."""
    centre_lat: float
    centre_lon: float

    @property
    def matrix(self):
        x = unit(self.centre_lat, self.centre_lon)
        z = np.array([0.0, 0.0, 1.0]) - x[2] * x
        z = z / np.linalg.norm(z)
        return np.stack((x, np.cross(z, x), z))

    def to_grid(self, lat_deg, lon_deg):
        return latlon(unit(np.asarray(lat_deg, float), np.asarray(lon_deg, float)) @ self.matrix.T)

    def to_true(self, glat_deg, glon_deg):
        return latlon(unit(np.asarray(glat_deg, float), np.asarray(glon_deg, float)) @ self.matrix)

    def vector_to_grid(self, glat_deg, glon_deg, east, north):
        """True eastward and northward components at grid places, as the grid's own eastward and northward ones."""
        lat, lon = self.to_true(glat_deg, glon_deg)
        e_true, n_true = tangents(lat, lon)
        e_grid, n_grid = tangents(np.asarray(glat_deg, float), np.asarray(glon_deg, float))
        e_grid, n_grid = e_grid @ self.matrix, n_grid @ self.matrix          # back in the true frame
        ee, ne = (e_true * e_grid).sum(-1), (n_true * e_grid).sum(-1)
        en, nn = (e_true * n_grid).sum(-1), (n_true * n_grid).sum(-1)
        return east * ee + north * ne, east * en + north * nn


def levels(depth_m, top=5.0, growth=1.2, largest=500.0):
    """Layer thicknesses (m) from the surface down: `top` at the surface, growing by `growth` a layer up to
    `largest`, until they hold `depth_m`."""
    out, dz = [], top
    while sum(out) < depth_m:
        out.append(round(min(dz, largest), 2))
        dz *= growth
    return np.array(out)


def size_h(nx, ny, nr, npx, npy, olx=4):
    if nx % npx or ny % npy:
        raise ValueError(f"{nx}x{ny} does not divide into {npx}x{npy} tiles")
    return f"""C SIZE.h for a Terluna sea with levels: {nx} x {ny} x {nr} cells on {npx} x {npy} processes
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
     &           Nr  = {nr})
      INTEGER MAX_OLX, MAX_OLY
      PARAMETER ( MAX_OLX = OLx, MAX_OLY = OLy )
"""


def model_options(source: Path) -> str:
    """The model's options with z-star, the implicit bed drag's code and sunlight that heats below the top layer."""
    text, n = re.subn(r"#undef\s+SHORTWAVE_HEATING\b", "#define SHORTWAVE_HEATING", cpp_options(source, True))
    if n != 1:
        raise RuntimeError("CPP_OPTIONS.h lacks the SHORTWAVE_HEATING switch this setup expects")
    return text


def exf_options(source: Path, tides=False) -> str:
    """The forcing package's options with its bulk formulas compiled out: stress, heat and fresh water are read. With
    `tides`, its reader of the tide's geopotential is compiled in."""
    text = (source / "pkg" / "exf" / "EXF_OPTIONS.h").read_text()
    if tides:
        text, n = re.subn(r"^#undef\s+EXF_ALLOW_TIDES\b", "#define EXF_ALLOW_TIDES", text, flags=re.M)
        if n != 1:
            raise RuntimeError("EXF_OPTIONS.h lacks the EXF_ALLOW_TIDES switch this setup expects")
    for flag in ("ALLOW_ATM_TEMP", "ALLOW_ATM_WIND", "ALLOW_DOWNWARD_RADIATION"):
        text, n = re.subn(rf"^#define\s+{flag}\b", f"#undef {flag}", text, flags=re.M)
        if n != 1:
            raise RuntimeError(f"EXF_OPTIONS.h lacks the {flag} switch this setup expects")
    for flag in ("ALLOW_RUNOFF", "ALLOW_CLIMSST_RELAXATION"):
        if not re.search(rf"^#define\s+{flag}\b", text, re.M):
            raise RuntimeError(f"EXF_OPTIONS.h no longer defines {flag}")
    return text


@dataclass
class Sea(Experiment):
    """One MITgcm configuration of a sea with levels, temperature and salt, in runs/mitgcm/<name>/."""
    delr: np.ndarray = None             # layer thicknesses from the surface down, m
    coriolis: dict = None               # centre, corner (south-west corners) and cosine fields, 1/s, on the grid
    theta0: float = 27.0                # starting temperature, C, everywhere
    salt0: float = 35.0                 # starting salinity, g/kg, everywhere
    rho_const: float = 1022.0           # reference density, kg/m3
    forcing: dict = field(default_factory=dict)   # any of FORCING: arrays (records, ny, nx); fluxes positive upward
    runoff: np.ndarray = None           # m/s into each cell, constant
    forcing_period_s: float = 0.0       # interval between the forcing's records
    sst_relax_s: float = 0.0            # time scale drawing the surface layer to climsst; 0 for none
    visc_h: float = 300.0               # background lateral viscosity, m2/s
    diff_r: float = 1.0e-5              # background vertical diffusivity of heat and salt, m2/s
    visc_r: float = 1.0e-4              # background vertical viscosity, m2/s
    surface_fields: tuple = ("ETAN", "MXLDEPTH")
    top_fields: tuple = ("SALT", "THETA", "UVEL", "VVEL")
    deep_fields: tuple = ("SALT", "THETA", "UVEL", "VVEL", "WVEL")
    deep_output_s: float = 0.0          # interval of the whole-depth means; 0 for none
    visc_grid: float = 0.0              # lateral viscosity as a share of each cell's L^2 / (4 dt); for mixed cell sizes
    start: dict = field(default_factory=dict)     # optional starting fields (levels, ny, nx): theta, salt, dye
    tide_m: np.ndarray = None           # the tide's equilibrium height on the grid (records, ny, nx), m; None for no tide
    tide_period_s: float = 86400.0      # interval between the tide's records; the series repeats after its last
    dye: bool = True                    # carry the old-water dye (TRAC01)
    dye_surface_s: float = 86400.0      # time scale drawing the top layer's dye to 0
    tracers: tuple = ()                 # Tracer objects carried in place of the dye (TRAC01, TRAC02, ...)
    tracers_from: int = 0               # the iteration at which the tracers start from their `start`
    idemix: dict = None                 # tidal and wind: the power put into internal waves (ny, nx), W/m2
    geothermal: np.ndarray = None       # the floor's heat into the sea (ny, nx), W/m2
    extra_lists: tuple = ()             # further outputs: (field names, file stem, interval s, level or None)

    def __post_init__(self):
        self.nonlinear = self.nonlinear_surface = True
        if self.tracers and self.dye:
            raise ValueError("a sea carries either the dye or its own tracers")

    def tracer_names(self):
        n = len(self.tracers) if self.tracers else int(self.dye)
        return tuple(f"TRAC{i:02d}" for i in range(1, n + 1))

    def code(self) -> dict:
        source = builds.checkout("mitgcm")
        nr = len(self.delr)
        options = model_options(source)
        if self.geothermal is not None:
            options = options_with(options, "ALLOW_GEOTHERMAL_FLUX", "CPP_OPTIONS.h")
        extra = sum(len(names) for names, *_ in self.extra_lists)
        code = {"SIZE.h": size_h(self.grid.nx, self.grid.ny, nr, *self.tiles), "packages.conf": PACKAGES_CONF,
                "CPP_OPTIONS.h": options,
                "EXF_OPTIONS.h": exf_options(source, tides=self.tide_m is not None),
                # lists written at different intervals each hold their own levels, the top layer's a whole column
                "DIAGNOSTICS_SIZE.h": diagnostics_size(source, len(self.top_fields) + len(self.deep_fields) + 2
                                                       + (len(self.tracers) + extra if self.tracers or extra else 0))}
        if self.tracers:
            text = (source / "pkg" / "rbcs" / "RBCS_SIZE.h").read_text()
            text, n = re.subn(r"PARAMETER\(\s*maskLEN\s*=\s*3\s*\)", f"PARAMETER( maskLEN = {2 + len(self.tracers)} )", text)
            if n != 1:
                raise RuntimeError("RBCS_SIZE.h no longer sets maskLEN as this setup expects")
            code["RBCS_SIZE.h"] = text
            text = (source / "pkg" / "ptracers" / "PTRACERS_SIZE.h").read_text()
            text, n = re.subn(r"PARAMETER\(\s*PTRACERS_num\s*=\s*1\s*\)", f"PARAMETER(PTRACERS_num = {len(self.tracers)} )",
                              text)
            if n != 1:
                raise RuntimeError("PTRACERS_SIZE.h no longer sets PTRACERS_num as this setup expects")
            code["PTRACERS_SIZE.h"] = text
        if self.idemix is not None:
            code["GGL90_OPTIONS.h"] = options_with((source / "pkg" / "ggl90" / "GGL90_OPTIONS.h").read_text(),
                                                   "ALLOW_GGL90_IDEMIX", "GGL90_OPTIONS.h")
        return code

    def inputs(self) -> tuple[dict, dict]:
        g, nr = self.grid, len(self.delr)
        fields = {"bathy.bin": -g.depth, "fCoriC.bin": self.coriolis["centre"], "fCoriG.bin": self.coriolis["corner"],
                  "fCorCs.bin": self.coriolis["cosine"]}
        parm01 = dict(
            gravity=lunar.GRAVITY, gBaro=lunar.GRAVITY, rhoConst=self.rho_const, rhoNil=self.rho_const,
            rotationPeriod=lunar.ROTATION_PERIOD_S, selectCoriMap=3,
            tRef=f"{nr}*{self.theta0:.6g}", sRef=f"{nr}*{self.salt0:.6g}",
            eosType="JMD95Z", selectP_inEOS_Zc=0,
            implicitFreeSurface=True, rigidLid=False, nonlinFreeSurf=4, select_rStar=2, hFacInf=0.1, hFacSup=5.0,
            exactConserv=True, useRealFreshWaterFlux=True,
            implicSurfPress=1.0, implicDiv2DFlow=1.0,
            momAdvection=True, momViscosity=True, viscAh=self.visc_h, viscC2smag=self.viscosity_smag,
            viscAhGrid=self.visc_grid, viscAhGridMax=0.5, useFullLeith=False,
            viscAr=self.visc_r, diffKrT=self.diff_r, diffKrS=self.diff_r,
            implicitDiffusion=True, implicitViscosity=True,
            bottomDragQuadratic=self.drag, no_slip_sides=False, no_slip_bottom=True,
            tempStepping=True, saltStepping=True, tempAdvScheme=33, saltAdvScheme=33, staggerTimeStep=True,
            hFacMin=0.2, hFacMinDr=float(self.delr[0]),
            readBinaryPrec=64, writeBinaryPrec=32, useSingleCpuIO=True, globalFiles=True)
        if self.tide_m is not None:
            parm01["momTidalForcing"] = True
        text = namelist("PARM01", parm01)
        text += namelist("PARM02", dict(cg2dTargetResidual=self.cg2d_residual, cg2dMaxIters=2000))
        text += "PARM03_PLACEHOLDER"
        text += namelist("PARM04", dict(usingSphericalPolarGrid=True, rSphere=lunar.RADIUS,
                                        xgOrigin=float(g.xg[0]), ygOrigin=float(g.yg[0]),
                                        delX=np.diff(g.xg), delY=np.diff(g.yg), delR=self.delr))
        files = dict(bathyFile="bathy.bin")
        for key, name in (("theta", "hydrogThetaFile"), ("salt", "hydrogSaltFile")):
            if key in self.start:
                fields[f"start_{key}.bin"] = self.start[key]
                files[name] = f"start_{key}.bin"
        if self.geothermal is not None:
            fields["geothermal.bin"] = np.where(g.wet, self.geothermal, 0.0)
            files["geothermalFile"] = "geothermal.bin"
        text += namelist("PARM05", files)
        ggl90 = dict(GGL90alpha=30.0, mxlMaxFlag=2, GGL90writeState=False)
        waves = ""
        if self.idemix is not None:     # the switch in the first group, the inputs in IDEMIX's own second group
            fields["idemix_tidal.bin"] = np.where(g.wet, self.idemix["tidal"], 0.0)
            fields["idemix_wind.bin"] = np.where(g.wet, self.idemix["wind"], 0.0)
            ggl90.update(useIDEMIX=True)
            waves = namelist("GGL90_PARM02", dict(IDEMIX_tidal_file="idemix_tidal.bin",
                                                  IDEMIX_wind_file="idemix_wind.bin"))
        texts = {"data.template": text,
                 "data.pkg": namelist("PACKAGES", dict(useEXF=True, useCAL=False, useDiagnostics=True, useGGL90=True)),
                 "eedata": " &EEPARMS\n nTx=1,\n nTy=1,\n &\n",
                 "data.ggl90": namelist("GGL90_PARM01", ggl90) + waves}
        nml1 = dict(useExfCheckRange=False, exf_iprec=64, useExfYearlyFields=False,
                    readStressOnAgrid=True)             # the stress is given at the cells' centres
        nml2 = {}
        if self.forcing:
            records = {len(v) for v in self.forcing.values()}
            if len(records) != 1 or set(self.forcing) - set(FORCING):
                raise ValueError("the forcing fields must share one number of records and be among " + ", ".join(FORCING))
            nml1["repeatPeriod"] = float(records.pop() * self.forcing_period_s)
            for name, values in self.forcing.items():
                fields[f"{name}.bin"] = values
                nml2.update({f"{name}file": f"{name}.bin", f"{name}period": float(self.forcing_period_s),
                             f"{name}StartTime": 0.5 * self.forcing_period_s})   # a record is its interval's mean
            if "climsst" in self.forcing:
                nml2["climsstTauRelax"] = float(self.sst_relax_s)
        if self.runoff is not None:
            fields["runoff.bin"] = self.runoff
            nml2.update(runofffile="runoff.bin", runoffperiod=0.0)
        nml3 = {}
        if self.tide_m is not None:     # the height at each record's instant, on the grid; the model takes -g times it
            fields["tidepot.bin"] = self.tide_m
            nml2.update(tidePotFile="tidepot.bin", tidePotPeriod=float(self.tide_period_s), tidePotStartTime=0.0,
                        tidePotRepCycle=float(len(self.tide_m) * self.tide_period_s))
            nml3["exf_inscal_tidePot"] = -lunar.GRAVITY
        if self.dye:
            mask = np.zeros((nr,) + g.depth.shape)
            mask[0] = g.wet
            fields["rbcs_mask.bin"], fields["rbcs_zero.bin"] = mask, np.zeros_like(mask)
            texts["data.pkg"] = namelist("PACKAGES", dict(useEXF=True, useCAL=False, useDiagnostics=True,
                                                          useGGL90=True, usePTRACERS=True, useRBCS=True))
            texts["data.ptracers"] = namelist("PTRACERS_PARM01", {
                "PTRACERS_numInUse": 1, "PTRACERS_Iter0": 0, "PTRACERS_advScheme(1)": 33,
                "PTRACERS_diffKh(1)": 0.0, "PTRACERS_diffKr(1)": self.diff_r,
                f"PTRACERS_ref(1:{nr},1)": f"{nr}*1.", "PTRACERS_names(1)": "old",
                **({"PTRACERS_initialFile(1)": "start_dye.bin"} if "dye" in self.start else {})})
            if "dye" in self.start:
                fields["start_dye.bin"] = self.start["dye"]
            texts["data.rbcs"] = (namelist("RBCS_PARM01", {"relaxMaskFile(3)": "rbcs_mask.bin"})
                                  + namelist("RBCS_PARM02", {"useRBCpTrNum(1)": True,
                                                             "tauRelaxPTR(1)": float(self.dye_surface_s),
                                                             "relaxPtracerFile(1)": "rbcs_zero.bin"}))
        elif self.tracers:
            wet3d = (np.cumsum(self.delr)[:, None, None] - self.delr[:, None, None] < g.depth[None]) & g.wet[None]
            ptr, rb1, rb2 = {"PTRACERS_numInUse": len(self.tracers), "PTRACERS_Iter0": int(self.tracers_from)}, {}, {}
            for i, tr in enumerate(self.tracers, start=1):
                start = np.broadcast_to(np.asarray(tr.start, dtype=float), wet3d.shape)
                fields[f"start_tr{i:02d}.bin"] = np.where(wet3d, start, 0.0)
                fields[f"rbcs_mask_tr{i:02d}.bin"], fields[f"rbcs_target_tr{i:02d}.bin"] = tr.relaxation(wet3d)
                ptr.update({f"PTRACERS_advScheme({i})": 33, f"PTRACERS_diffKh({i})": 0.0,
                            f"PTRACERS_diffKr({i})": self.diff_r, f"PTRACERS_names({i})": tr.name,
                            f"PTRACERS_initialFile({i})": f"start_tr{i:02d}.bin"})
                rb1[f"relaxMaskFile({2 + i})"] = f"rbcs_mask_tr{i:02d}.bin"
                rb2.update({f"useRBCpTrNum({i})": True, f"tauRelaxPTR({i})": float(tr.top_s),
                            f"relaxPtracerFile({i})": f"rbcs_target_tr{i:02d}.bin"})
            texts["data.pkg"] = namelist("PACKAGES", dict(useEXF=True, useCAL=False, useDiagnostics=True,
                                                          useGGL90=True, usePTRACERS=True, useRBCS=True))
            texts["data.ptracers"] = namelist("PTRACERS_PARM01", ptr)
            texts["data.rbcs"] = namelist("RBCS_PARM01", rb1) + namelist("RBCS_PARM02", rb2)
        texts["data.exf"] = (namelist("EXF_NML_01", nml1) + namelist("EXF_NML_02", nml2)
                             + namelist("EXF_NML_03", nml3) + namelist("EXF_NML_04", {}))
        more = {} if self.tide_m is None else dict(tide=dict(records=len(self.tide_m), period_s=self.tide_period_s))
        if self.tracers:
            more["tracers"] = [tr.name for tr in self.tracers]
            more["tracers_from_iteration"] = int(self.tracers_from)
        if self.idemix is not None:
            area = g.area * g.wet
            more["idemix_mw"] = {k: float((np.where(g.wet, v, 0.0) * area).sum() / 1e6) for k, v in self.idemix.items()}
        if self.geothermal is not None:
            more["geothermal_w_m2"] = float((np.where(g.wet, self.geothermal, 0.0) * g.area).sum() / g.area[g.wet].sum())
        texts["experiment.json"] = json.dumps(dict(
            name=self.name, nx=g.nx, ny=g.ny, nr=nr, wet_columns=int(g.wet.sum()), deltat=self.deltat,
            delr=[float(x) for x in self.delr], theta0=self.theta0, salt0=self.salt0, rho_const=self.rho_const,
            drag=self.drag, viscosity_smag=self.viscosity_smag, visc_h=self.visc_h, visc_r=self.visc_r,
            visc_grid=self.visc_grid, start=sorted(self.start),
            diff_r=self.diff_r, tiles=list(self.tiles), forcing=sorted(self.forcing),
            forcing_period_s=self.forcing_period_s, sst_relax_s=self.sst_relax_s, rivers=self.runoff is not None,
            closure="ggl90", equation_of_state="JMD95Z", dye=self.dye, dye_surface_s=self.dye_surface_s,
            gravity=lunar.GRAVITY, radius=lunar.RADIUS,
            rotation_period_s=lunar.ROTATION_PERIOD_S, **more), indent=1) + "\n"
        return texts, fields

    def write_diagnostics(self, active=True):
        """Means over each output interval: the surface's level and mixed layer, the top layer's water and flow, and,
        when asked, every level's."""
        dye = self.tracer_names()
        lists = [(self.surface_fields, "surface", self.output_s, None), (self.top_fields, "top", self.output_s, 1)]
        if self.deep_output_s > 0:
            lists.append((self.deep_fields + dye, "deep", self.deep_output_s, None))
        lists += list(self.extra_lists)
        text = " &DIAGNOSTICS_LIST\n"
        for n, (names, stem, interval, level) in enumerate(lists, start=1):
            text += f"  fields(1:{len(names)},{n}) = " + ",".join(f"'{name:<8s}'" for name in names) + ",\n"
            if level:
                text += f"  levels(1,{n}) = {float(level):.1f},\n"
            text += f"  fileName({n}) = '{stem}',\n  frequency({n}) = {interval if active else 1.0e12:.1f},\n"
        (self.folder / "data.diagnostics").write_text(text + " &\n &DIAG_STATIS_PARMS\n &\n")

    def restart_from(self, parent: Path, iteration: int, build: str = None):
        """Begin this run from another run's pickup at `iteration` (same grid and levels): copy its state and its
        turbulence into this run's folder, with the internal waves' energy at zero when this run has IDEMIX, and stamp
        the folder with this run's setup on `build` (its own binary's by default). Call after write(fresh=True); the
        tracers start from their files."""
        import shutil
        if self.tracers and self.tracers_from != iteration:
            raise ValueError("the tracers must start at the restart's iteration")
        parent = Path(parent)
        stem = f"{iteration:010d}"
        names = [f"pickup.{stem}.data", f"pickup.{stem}.meta"]
        if self.dye:                            # the dye continues from the parent's own
            names += [f"pickup_ptracers.{stem}.data", f"pickup_ptracers.{stem}.meta"]
        for name in names:
            shutil.copyfile(parent / name, self.folder / name)
        tke = (parent / f"pickup_ggl90.{stem}.data").read_bytes()
        meta = (parent / f"pickup_ggl90.{stem}.meta").read_text()
        if self.idemix is not None:            # the closure then reads a second record, the waves' energy
            tke += bytes(len(tke))
            meta, n = re.subn(r"nrecords = \[\s*1 \]", "nrecords = [          2 ]", meta)
            if n != 1:
                raise RuntimeError(f"{parent} holds a turbulence pickup of another layout")
        (self.folder / f"pickup_ggl90.{stem}.data").write_bytes(tke)
        (self.folder / f"pickup_ggl90.{stem}.meta").write_text(meta)
        (self.folder / "setup.sha256").write_text(self.setup_digest(build or self.binary().parent.name) + "\n")
        (self.folder / "restart.json").write_text(json.dumps(dict(parent=str(parent), iteration=iteration), indent=1) + "\n")
