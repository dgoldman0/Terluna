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

The old water. A dye starts at 1 everywhere and is drawn to 0 in the top layer within a day (pkg/ptracers with
pkg/rbcs), and is otherwise carried and mixed like salt. What is left of it at a depth is the share of that water
that has not been at the surface since the run began, so its fall measures how fast the deeps are renewed.
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


def exf_options(source: Path) -> str:
    """The forcing package's options with its bulk formulas compiled out: stress, heat and fresh water are read."""
    text = (source / "pkg" / "exf" / "EXF_OPTIONS.h").read_text()
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
    dye: bool = True                    # carry the old-water dye (TRAC01)
    dye_surface_s: float = 86400.0      # time scale drawing the top layer's dye to 0

    def __post_init__(self):
        self.nonlinear = self.nonlinear_surface = True

    def code(self) -> dict:
        source = builds.checkout("mitgcm")
        nr = len(self.delr)
        return {"SIZE.h": size_h(self.grid.nx, self.grid.ny, nr, *self.tiles), "packages.conf": PACKAGES_CONF,
                "CPP_OPTIONS.h": model_options(source), "EXF_OPTIONS.h": exf_options(source),
                # lists written at different intervals each hold their own levels, the top layer's a whole column
                "DIAGNOSTICS_SIZE.h": diagnostics_size(source, len(self.top_fields) + len(self.deep_fields) + 2)}

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
            viscAhGridMax=0.5, useFullLeith=False,
            viscAr=self.visc_r, diffKrT=self.diff_r, diffKrS=self.diff_r,
            implicitDiffusion=True, implicitViscosity=True,
            bottomDragQuadratic=self.drag, no_slip_sides=False, no_slip_bottom=True,
            tempStepping=True, saltStepping=True, tempAdvScheme=33, saltAdvScheme=33, staggerTimeStep=True,
            hFacMin=0.2, hFacMinDr=float(self.delr[0]),
            readBinaryPrec=64, writeBinaryPrec=32, useSingleCpuIO=True, globalFiles=True)
        text = namelist("PARM01", parm01)
        text += namelist("PARM02", dict(cg2dTargetResidual=self.cg2d_residual, cg2dMaxIters=2000))
        text += "PARM03_PLACEHOLDER"
        text += namelist("PARM04", dict(usingSphericalPolarGrid=True, rSphere=lunar.RADIUS,
                                        xgOrigin=float(g.xg[0]), ygOrigin=float(g.yg[0]),
                                        delX=np.diff(g.xg), delY=np.diff(g.yg), delR=self.delr))
        text += namelist("PARM05", dict(bathyFile="bathy.bin"))
        texts = {"data.template": text,
                 "data.pkg": namelist("PACKAGES", dict(useEXF=True, useCAL=False, useDiagnostics=True, useGGL90=True)),
                 "eedata": " &EEPARMS\n nTx=1,\n nTy=1,\n &\n",
                 "data.ggl90": namelist("GGL90_PARM01", dict(GGL90alpha=30.0, mxlMaxFlag=2, GGL90writeState=False))}
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
        if self.dye:
            mask = np.zeros((nr,) + g.depth.shape)
            mask[0] = g.wet
            fields["rbcs_mask.bin"], fields["rbcs_zero.bin"] = mask, np.zeros_like(mask)
            texts["data.pkg"] = namelist("PACKAGES", dict(useEXF=True, useCAL=False, useDiagnostics=True,
                                                          useGGL90=True, usePTRACERS=True, useRBCS=True))
            texts["data.ptracers"] = namelist("PTRACERS_PARM01", {
                "PTRACERS_numInUse": 1, "PTRACERS_Iter0": 0, "PTRACERS_advScheme(1)": 33,
                "PTRACERS_diffKh(1)": 0.0, "PTRACERS_diffKr(1)": self.diff_r,
                f"PTRACERS_ref(1:{nr},1)": f"{nr}*1.", "PTRACERS_names(1)": "old"})
            texts["data.rbcs"] = (namelist("RBCS_PARM01", {"relaxMaskFile(3)": "rbcs_mask.bin"})
                                  + namelist("RBCS_PARM02", {"useRBCpTrNum(1)": True,
                                                             "tauRelaxPTR(1)": float(self.dye_surface_s),
                                                             "relaxPtracerFile(1)": "rbcs_zero.bin"}))
        texts["data.exf"] = (namelist("EXF_NML_01", nml1) + namelist("EXF_NML_02", nml2)
                             + namelist("EXF_NML_03", {}) + namelist("EXF_NML_04", {}))
        texts["experiment.json"] = json.dumps(dict(
            name=self.name, nx=g.nx, ny=g.ny, nr=nr, wet_columns=int(g.wet.sum()), deltat=self.deltat,
            delr=[float(x) for x in self.delr], theta0=self.theta0, salt0=self.salt0, rho_const=self.rho_const,
            drag=self.drag, viscosity_smag=self.viscosity_smag, visc_h=self.visc_h, visc_r=self.visc_r,
            diff_r=self.diff_r, tiles=list(self.tiles), forcing=sorted(self.forcing),
            forcing_period_s=self.forcing_period_s, sst_relax_s=self.sst_relax_s, rivers=self.runoff is not None,
            closure="ggl90", equation_of_state="JMD95Z", dye=self.dye, dye_surface_s=self.dye_surface_s,
            gravity=lunar.GRAVITY, radius=lunar.RADIUS,
            rotation_period_s=lunar.ROTATION_PERIOD_S), indent=1) + "\n"
        return texts, fields

    def write_diagnostics(self, active=True):
        """Means over each output interval: the surface's level and mixed layer, the top layer's water and flow, and,
        when asked, every level's."""
        dye = ("TRAC01",) if self.dye else ()
        lists = [(self.surface_fields, "surface", self.output_s, None), (self.top_fields, "top", self.output_s, 1)]
        if self.deep_output_s > 0:
            lists.append((self.deep_fields + dye, "deep", self.deep_output_s, None))
        text = " &DIAGNOSTICS_LIST\n"
        for n, (names, stem, interval, level) in enumerate(lists, start=1):
            text += f"  fields(1:{len(names)},{n}) = " + ",".join(f"'{name:<8s}'" for name in names) + ",\n"
            if level:
                text += f"  levels(1,{n}) = {float(level):.1f},\n"
            text += f"  fileName({n}) = '{stem}',\n  frequency({n}) = {interval if active else 1.0e12:.1f},\n"
        (self.folder / "data.diagnostics").write_text(text + " &\n &DIAG_STATIS_PARMS\n &\n")
