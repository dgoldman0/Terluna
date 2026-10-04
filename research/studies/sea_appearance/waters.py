"""The colour of the Open Moon's waters and how deep their daylight reaches.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.waters

The biosphere's design guesses for the seas' contents (biosphere/living_water/waters.json) go through the
illumination domain's seawater optics (illumination/water_column/model.py) to a remote-sensing reflectance in each
of the solved sky's spectral channels. Multiplied by the daylight reaching the water, the solved spherical
sky's direct beam and diffuse light for the Sun 45 degrees up, it gives the light leaving the water: its
luminance and its chromaticity, under the Moon's sky and under the Earth control's. This is the water's own
colour, as seen looking down past the reflected sky; the reflected sky and glitter come with the results by
regime. The same optics give the depth at which the daylight's photosynthetic photons fall to 1%.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from illumination.water_column import model as optics
from illumination.water_column.fetch_inputs import INPUTS
from shared.constants import PLANCK, SPEED_OF_LIGHT, AVOGADRO
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / "results/water_colours.json"
WATERS = ROOT / "biosphere/living_water/waters.json"
SPHERICAL = ROOT / "research/runs/optical_comfort/spherical"
WORLDS = {"moon": SPHERICAL / "moon_1.2atm_standard.npz", "earth": SPHERICAL / "earth_control_standard.npz"}
SUN_ELEVATION = 45.0
BANDS = (412.0, 443.0, 490.0, 510.0, 555.0, 670.0)
DISPLAY = np.arange(400.0, 700.1, 2.0)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def daylight(path, elevation):
    """Downwelling daylight on the water in each channel, per unit source, with the channels' XYZ weights."""
    with np.load(path, allow_pickle=False) as z:
        a, sa = np.degrees(z["a"]), np.degrees(z["sa"])
        diffuse = np.array([np.interp(elevation, a, z["moments"][0, :, c, 5]) for c in range(z["moments"].shape[2])])
        direct = np.array([np.interp(elevation, sa, z["beam"][0, :, c, 1]) for c in range(z["beam"].shape[2])])
        return dict(wavelength=z["wavelength"], band=z["band"], xyz=z["xyz"], energy=z["energy"],
                    irradiance=direct + diffuse)


def euphotic_depth(light, kd):
    """Depth (m) where the photosynthetic photons, 400-700 nm, fall to 1% of their value under the surface."""
    inside = (light["wavelength"] >= 400) & (light["wavelength"] <= 700)
    photons = (light["irradiance"] * light["energy"] * light["wavelength"] * 1e-9 / (PLANCK * SPEED_OF_LIGHT))[inside]
    depth = np.linspace(0, 400, 40001)
    fraction = (photons[None, :] * np.exp(-np.outer(depth, kd[inside]))).sum(1) / photons.sum()
    surface_umol = float(photons.sum() / AVOGADRO * 1e6)
    return float(np.interp(-0.01, -fraction, depth)), surface_umol


def chromaticity(xyz):
    total = xyz.sum()
    return [round(float(xyz[0] / total), 4), round(float(xyz[1] / total), 4)]


def display_reflectance(water, chlorophyll, soil):
    """The reflectance on a regular 2 nm grid; the sky's channels are not ordered inside its absorption bands."""
    iop = optics.water(DISPLAY, chlorophyll, water["dissolved_organic_440_per_m"]["best"], water["fines_g_m3"]["best"], soil)
    return np.round(optics.remote_sensing_reflectance(iop["a"], iop["bb"]), 7).tolist()


def evaluate(name, water, chlorophyll, light, soil=None):
    w = light["wavelength"]
    iop = optics.water(w, chlorophyll, water["dissolved_organic_440_per_m"]["best"], water["fines_g_m3"]["best"],
                       soil or water["soil"])
    rrs = optics.remote_sensing_reflectance(iop["a"], iop["bb"])
    leaving = rrs * light["irradiance"]
    xyz = leaving @ light["xyz"]
    kd = optics.diffuse_attenuation(iop["a"], iop["bb"], 90 - SUN_ELEVATION)
    depth, _ = euphotic_depth(light, kd)
    return dict(xyz_cd_m2=[round(float(v), 4) for v in xyz], luminance_cd_m2=round(float(xyz[1]), 3),
                chromaticity_xy=chromaticity(xyz), depth_one_percent_photons_m=round(depth, 2),
                kd_490_per_m=round(float(np.interp(490, w, kd)), 4),
                rrs_bands_per_sr=dict(zip([f"{b:g}" for b in BANDS], np.round(np.interp(BANDS, w, rrs), 6).tolist())),
                rrs_per_sr=np.round(rrs, 7).tolist())


def build():
    record = json.loads(WATERS.read_text())
    lights = {world: daylight(path, SUN_ELEVATION) for world, path in WORLDS.items()}
    ratio = record["light_cycle"]["chlorophyll_dusk_over_dawn"]["best"]
    cases = []
    for name, water in record["waters"].items():
        chl = water["chlorophyll_mg_m3"]["best"]
        # The best guess is a mean over the cycle; dawn and dusk sit either side of it, their ratio as guessed.
        for phase, scale in (("mean", 1.0), ("dawn", 2 / (1 + ratio)), ("dusk", 2 * ratio / (1 + ratio))):
            cases.append(dict(water=name, phase=phase, soil=water["soil"], chlorophyll_mg_m3=round(chl * scale, 4),
                              rrs_display_per_sr=display_reflectance(water, chl * scale, water["soil"]),
                              **{world: evaluate(name, water, chl * scale, light) for world, light in lights.items()}))
    coast = record["waters"]["productive_coast"]
    for soil in record["soils"]:
        if soil != coast["soil"]:
            cases.append(dict(water="productive_coast", phase="mean", soil=soil,
                              chlorophyll_mg_m3=coast["chlorophyll_mg_m3"]["best"],
                              rrs_display_per_sr=display_reflectance(coast, coast["chlorophyll_mg_m3"]["best"], soil),
                              **{world: evaluate("productive_coast", coast, coast["chlorophyll_mg_m3"]["best"], light, soil)
                                 for world, light in lights.items()}))
    files = [Path(__file__), ROOT / "illumination/water_column/model.py"]
    surface = {world: round(euphotic_depth(light, np.zeros_like(light["wavelength"]))[1], 1)
               for world, light in lights.items()}
    return dict(
        schema="terluna.research.sea-water-colours/1",
        evidence=("The colour of the light leaving the Open Moon's waters under clear daylight, from design guesses "
                  "of their contents and Earth's seawater optics; a design-guess result, labelled as such."),
        reading_rule=("Remote-sensing reflectance per steradian in each channel of each world's solved sky "
                      "(channel_wavelength_nm); water-leaving XYZ in cd/m2 for the Sun 45 degrees up under "
                      "each world's clear sky, the water's own light without the reflected sky; chromaticity is CIE "
                      "1931 xy. rrs_display_per_sr repeats the reflectance on a regular 2 nm grid. Dawn and dusk scale the chlorophyll about its best guess by the light cycle's ratio. "
                      "Depths are where photosynthetic photons fall to 1% of their value just under the surface."),
        producer=dict(lane="research", runner=str(Path(__file__).relative_to(ROOT)),
                      files={str(p.relative_to(ROOT)): sha256(p) for p in files}, constants=constants_used(files)),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (WATERS, *WORLDS.values())} | {
            f"illumination/water_column/inputs/{f['name']}": f["sha256"] for f in INPUTS.manifest()["files"]},
        sun_elevation_deg=SUN_ELEVATION,
        channel_wavelength_nm={world: np.round(light["wavelength"], 3).tolist() for world, light in lights.items()},
        display_wavelength_nm=DISPLAY.tolist(),
        surface_photosynthetic_photons_umol_m2_s=surface,
        daylight_illuminance_lux={world: round(float(light["irradiance"] @ light["xyz"][:, 1]), 1)
                                  for world, light in lights.items()},
        cases=cases)


def main():
    product = build()
    RESULT.write_text(json.dumps(product, indent=1, allow_nan=False) + "\n")
    print("daylight", product["daylight_illuminance_lux"], "photons", product["surface_photosynthetic_photons_umol_m2_s"])
    for c in product["cases"]:
        m, e = c["moon"], c["earth"]
        print(f"{c['water']:17s} {c['phase']:5s} {c['soil']} chl {c['chlorophyll_mg_m3']:6.3f} | Moon L {m['luminance_cd_m2']:7.1f}"
              f" xy {m['chromaticity_xy']} 1% {m['depth_one_percent_photons_m']:6.1f} m | Earth L {e['luminance_cd_m2']:7.1f}"
              f" xy {e['chromaticity_xy']}")


if __name__ == "__main__":
    main()
