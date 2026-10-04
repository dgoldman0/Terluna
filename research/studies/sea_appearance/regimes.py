"""How the seas look by regime: the radiance and colour of sky and sea in every direction at four coasts.

    NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 <optical comfort environment>/bin/python \\
        -m research.studies.sea_appearance.regimes

(The sky's radiance in any direction comes from the solver's numba readout, as in the optical comfort study.)

At the four coasts chosen for the renderings, six moments of the shore month span the regimes: the Sun at its
highest; the Sun 4 degrees up in the afternoon, when its glitter runs to the horizon; 4, 15 and 35 degrees
below the horizon through the long evening; and the darkest hour of the night. For an eye 2 m above the water
the runner computes, at each moment, the radiance of sky and sea in every direction, in CIE XYZ:

- the sky: the solved spherical sky read from its cached scattering solution (illumination/sky/
  solved_transport.py) for the Sun and, weighted channel by channel by its spectrum, for the Earth
  (illumination/earthlight), with the calendar's 0.001-lux placeholder for stars and airglow;
- the sea: the sky reflected by rough water and the glitter of the Sun's and the Earth's disks
  (illumination/water_surface/reflection.py) with that hour's slopes from sea_slopes.json, and the water's own
  light (illumination/water_column) for the biosphere's productive coast.

The sea fills every direction below the horizon: the coastline is left to the renderings. The spread of
reflected radiance over the facets, relative to its mean, measures how strongly the waves show.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from illumination.water_column import model as water_optics
from illumination.water_surface import reflection
from research.studies.sea_appearance import lighting
from shared.constants import AU, EARTH_RADIUS, MOON_RADIUS, SUN_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / "results/regimes.json"
PANORAMAS = ROOT / "research/runs/sea_appearance/regimes.npz"
CALENDAR = HERE / "results/lighting_calendar.json"
SLOPES = HERE / "results/sea_slopes.json"
SLOPE_SERIES = HERE / "results/sea_slopes.npz"
WATERS = ROOT / "biosphere/living_water/waters.json"
SOLUTION = lighting.SPHERICAL / "moon_1.2atm_standard.npz"
COASTS = ("W Procellarum", "Smythii headland", "S Nubium", "Ingenii coast")
MOMENTS = (("noon", None), ("low Sun", 4.0), ("after sunset", -4.0), ("evening", -15.0),
           ("deep evening", -35.0), ("darkest", None))
WATER = "productive_coast"
EYE_HEIGHT_M = 2.0
AZIMUTHS = np.arange(0.5, 360.0, 1.0)
ELEVATIONS = np.arange(29.75, -30.0, -0.5)
FIELD_ELEVATIONS = np.r_[np.arange(0.0, 10.0, 0.25), np.arange(10.0, 30.0, 1.0), np.arange(30.0, 90.1, 2.5)]
FIELD_AZIMUTHS = np.arange(0.0, 180.1, 2.0)
D65 = (0.3127, 0.3290)
ORDER = 24


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(elevation_deg, azimuth_deg):
    """Unit vectors (east, north, up) for elevations and azimuths clockwise from north."""
    e, a = np.radians(elevation_deg), np.radians(azimuth_deg)
    return np.stack(np.broadcast_arrays(np.cos(e) * np.sin(a), np.cos(e) * np.cos(a), np.sin(e)), -1)


def moment(series, kind, target, total):
    """Hour index of a moment: the highest Sun, the darkest hour, or the Sun sinking through an elevation."""
    sun = np.asarray(series["sun_elevation_deg"])
    if kind == "noon":
        return int(np.argmax(sun))
    if kind == "darkest":
        return int(np.argmin(total))
    noon, n = int(np.argmax(sun)), len(sun)
    for step in range(1, n):
        i, j = (noon + step - 1) % n, (noon + step) % n
        if sun[i] >= target > sun[j]:
            return j
    raise ValueError(f"The Sun never sinks through {target} degrees")


class Fields:
    """The solved sky's radiance by direction for a source at any elevation, in XYZ per unit of its spectrum."""

    def __init__(self, path=SOLUTION):
        from illumination.sky.solved_transport import evaluate
        with np.load(path, allow_pickle=False) as z:
            self.solution = {k: z[k] for k in z.files if k != "meta"}
            self.meta = json.loads(str(z["meta"]))
        self.evaluate = evaluate
        self.xyz = self.solution["xyz"]

    def field(self, elevation_deg, weights=None):
        radiance = self.evaluate(self.solution, self.meta, np.array([elevation_deg]), FIELD_ELEVATIONS, FIELD_AZIMUTHS)[0]
        if weights is not None:
            radiance = radiance * weights
        return radiance @ self.xyz                                         # (elevations, azimuths, 3)


def lookup(field, elevation, relative_azimuth):
    """Bilinear readout of a field at elevations (deg, >= 0) and azimuths from its source (deg, 0-180)."""
    e = np.clip(elevation, FIELD_ELEVATIONS[0], FIELD_ELEVATIONS[-1])
    a = np.clip(relative_azimuth, 0.0, 180.0)
    i = np.clip(np.searchsorted(FIELD_ELEVATIONS, e) - 1, 0, len(FIELD_ELEVATIONS) - 2)
    j = np.clip(np.searchsorted(FIELD_AZIMUTHS, a) - 1, 0, len(FIELD_AZIMUTHS) - 2)
    u = ((e - FIELD_ELEVATIONS[i]) / (FIELD_ELEVATIONS[i + 1] - FIELD_ELEVATIONS[i]))[..., None]
    v = ((a - FIELD_AZIMUTHS[j]) / (FIELD_AZIMUTHS[j + 1] - FIELD_AZIMUTHS[j]))[..., None]
    return ((1 - u) * ((1 - v) * field[i, j] + v * field[i, j + 1]) +
            u * ((1 - v) * field[i + 1, j] + v * field[i + 1, j + 1]))


def sky_function(sources, floor_xyz):
    """XYZ sky radiance for direction vectors, summed over (field, source azimuth) pairs."""
    def sky(d):
        elevation = np.degrees(np.arcsin(np.clip(d[:, 2], -1, 1)))
        azimuth = np.degrees(np.arctan2(d[:, 0], d[:, 1])) % 360
        total = np.tile(floor_xyz, (len(d), 1))
        for field, source_azimuth in sources:
            relative = np.abs((azimuth - source_azimuth + 180) % 360 - 180)
            total = total + lookup(field, elevation, relative)
        return total
    return sky


def normal_beam(sky, elevation_deg):
    """Direct irradiance normal to the beam per channel and unit source, for a source at an elevation."""
    if elevation_deg <= -0.3:
        return np.zeros(sky.direct_channels.shape[1])
    horizontal = np.array([np.interp(elevation_deg, sky.beam_angles, sky.direct_channels[:, c])
                           for c in range(sky.direct_channels.shape[1])])
    return horizontal / np.sin(np.radians(max(elevation_deg, 0.5)))


def chromaticity(xyz):
    total = float(np.sum(xyz))
    return [round(float(xyz[0]) / total, 4), round(float(xyz[1]) / total, 4)] if total > 0 else None


def scene(coast, kind, target, calendar, slopes, series, fields, sky, earth_light, water, n_water, hours, jd):
    site = next(s for s in calendar["sites"] if s["short"] == coast)
    lon, lat = site["lon_lat"]
    s = site["series"]
    total = np.array(s["ground_lux_from_sun"]) + np.array(s["ground_lux_from_earth"]) + lighting.FLOOR_LUX
    index = moment(s, kind, target, total)
    g = lighting.geometry(jd[index:index + 1], lon, lat)
    sun_el, sun_az = float(g["sun_elevation"][0]), float(g["sun_azimuth"][0])
    earth_el, earth_az = float(g["earth_elevation"][0]), float(g["earth_azimuth"][0])
    phase = float(g["earth_phase_angle"][0])
    weights = earth_light.weights(g["earth_phase_angle"], g["earth_distance_m"], g["earth_sunlight_factor"])[0]
    # Slopes at that hour: the short waves hourly, the resolved waves from the nearest spectrum within a day, or the
    # nearside shores' median where no wave run reaches.
    key = coast.replace(" ", "_").lower()
    short = series[f"{key}/short_covariance"][index].reshape(2, 2).astype(float)
    resolved_hours = series[f"{key}/resolved_hours"]
    if len(resolved_hours) and np.min(np.abs(resolved_hours - hours[index])) <= 24:
        k = int(np.argmin(np.abs(resolved_hours - hours[index])))
        resolved = series[f"{key}/resolved_covariance"][k].reshape(2, 2).astype(float)
        resolved_source = f"spectrum {abs(resolved_hours[k] - hours[index]):.0f} h from the moment"
    else:
        stack = [series[f"{c.replace(' ', '_').lower()}/resolved_covariance"].reshape(-1, 2, 2) for c in
                 ("W Procellarum", "S Imbrium", "S Nubium", "Nectaris")]
        resolved = np.median(np.concatenate(stack), axis=0).astype(float)
        resolved_source = "median of the nearside shores' spectra (no wave run here)"
    covariance = short + resolved
    # The sky: the Sun, the Earth (while it lights the air) and the placeholder for stars and airglow.
    floor = lighting.FLOOR_LUX / np.pi * np.array([D65[0] / D65[1], 1.0, (1 - sum(D65)) / D65[1]])
    sources = [(fields.field(sun_el), sun_az)]
    if earth_el > -30:
        sources.append((fields.field(earth_el, weights), earth_az))
    sky_xyz = sky_function(sources, floor)
    # Direct beams for the glitter, the water's own light from the downwelling light of both sources.
    sun_normal = normal_beam(sky, sun_el) @ fields.xyz
    earth_normal = (normal_beam(sky, earth_el) * weights) @ fields.xyz
    downwelling = sky.ground_channels(sun_el)[0] + weights * sky.ground_channels(earth_el)[0]
    w = water["waters"][WATER]
    iop = water_optics.water(sky.channel_wavelength, w["chlorophyll_mg_m3"]["best"],
                             w["dissolved_organic_440_per_m"]["best"], w["fines_g_m3"]["best"], w["soil"])
    rrs = water_optics.remote_sensing_reflectance(iop["a"], iop["bb"])
    leaving_nadir = (rrs * downwelling) @ fields.xyz
    # The panorama: sky above, sea below the dip of the horizon.
    dip = np.degrees(np.sqrt(2 * EYE_HEIGHT_M / MOON_RADIUS))
    ee, aa = np.meshgrid(ELEVATIONS, AZIMUTHS, indexing="ij")
    image = np.zeros(ee.shape + (3,))
    above = ee > -dip
    image[above] = sky_xyz(unit(ee[above], aa[above]))
    below = ~above
    views = unit(-ee[below], aa[below] + 180.0)                 # from the water toward the eye
    reflected, square = reflection.sky_reflection_batch(views, covariance, n_water, sky_xyz, order=ORDER)
    sun_dir, earth_dir = unit(sun_el, sun_az), unit(earth_el, earth_az)
    sun_radius = np.arcsin(SUN_RADIUS / (AU * float(g["earth_sunlight_factor"][0]) ** -0.5))
    earth_radius = np.arcsin(EARTH_RADIUS / float(g["earth_distance_m"][0]))
    glitter_sun = reflection.disk_glint(sun_dir, sun_radius, views, covariance, n_water)[:, None] * sun_normal
    glitter_earth = reflection.disk_glint(earth_dir, earth_radius, views, covariance, n_water)[:, None] * earth_normal
    transmit = (1 - reflection.fresnel(views[:, 2], n_water)) / (1 - reflection.fresnel(1.0, n_water))
    leaving = transmit[:, None] * leaving_nadir
    image[below] = reflected + glitter_sun + glitter_earth + leaving
    contrast = np.zeros(ee.shape)
    contrast[below] = np.sqrt(np.maximum(square[:, 1] - reflected[:, 1] ** 2, 0)) / np.maximum(reflected[:, 1], 1e-30)
    # Readouts: toward the brighter source and away from it.
    bright_az = sun_az if sun_el > -40 or earth_el <= 0 else earth_az

    def at(elevation, azimuth):
        r = int(np.argmin(np.abs(ELEVATIONS - elevation)))
        c = int(np.argmin(np.abs((AZIMUTHS - azimuth + 180) % 360 - 180)))
        return image[r, c], contrast[r, c]

    readout = {}
    for name, azimuth in (("toward", bright_az), ("away", (bright_az + 180) % 360)):
        for label, elevation in (("sky_30", 29.75), ("sky_5", 4.75), ("sky_1", 0.75), ("sea_1", -0.75),
                                 ("sea_5", -4.75), ("sea_20", -19.75)):
            xyz, cv = at(elevation, azimuth)
            readout[f"{name}/{label}"] = dict(luminance_cd_m2=round(float(xyz[1]), 6), chromaticity_xy=chromaticity(xyz),
                                              **({"wave_contrast": round(float(cv), 4)} if label.startswith("sea") else {}))
    glitter = np.zeros(ee.shape)
    glitter[below] = glitter_sun[:, 1] + glitter_earth[:, 1]
    peak = np.unravel_index(int(np.argmax(glitter)), glitter.shape)
    bright = (glitter > reflected_map(reflected, below, ee.shape)) & below
    record = dict(
        coast=coast, moment=kind, hour=float(hours[index]), date_tt=float(jd[index]),
        sun=dict(elevation_deg=round(sun_el, 2), azimuth_deg=round(sun_az, 1)),
        earth=dict(elevation_deg=round(earth_el, 2), azimuth_deg=round(earth_az, 1), phase_angle_deg=round(phase, 1),
                   lit_fraction=round(float(g["earth_lit_fraction"][0]), 3)),
        light_on_ground_lux=round(float(total[index]), 4),
        slopes=dict(mean_square=round(float(np.trace(covariance)), 5), short_waves=round(float(np.trace(short)), 5),
                    resolved_waves=round(float(np.trace(resolved)), 5), resolved_source=resolved_source),
        water=WATER, water_leaving_nadir_cd_m2=round(float(leaving_nadir[1]), 6),
        glitter=dict(peak_cd_m2=round(float(glitter[peak]), 3), peak_elevation_deg=float(ELEVATIONS[peak[0]]),
                     peak_azimuth_deg=float(AZIMUTHS[peak[1]]),
                     share_of_sea_outshining_the_reflected_sky=round(float(bright.sum() / below.sum()), 5)),
        readout=readout)
    return record, image.astype(np.float32), contrast.astype(np.float32)


def reflected_map(reflected, below, shape):
    out = np.zeros(shape)
    out[below] = reflected[:, 1]
    return out


def build():
    calendar = json.loads(CALENDAR.read_text())
    slopes = json.loads(SLOPES.read_text())
    water = json.loads(WATERS.read_text())
    month = json.loads(lighting.SHORE_MONTH.read_text())
    hours = np.array(calendar["time_hours"], float)
    jd = month["tide_month"]["middle_jd_tt"] + (hours - 1.5 * lighting.CYCLE_HOURS) / 24
    sky = lighting.Sky()
    sky.channel_wavelength = np.load(SOLUTION, allow_pickle=False)["wavelength"]
    earth_light = lighting.Earthlight(sky)
    fields = Fields()
    n_water = float(reflection.refractive_index(550.0))
    with np.load(SLOPE_SERIES) as z:
        series = {k: z[k] for k in z.files}
    records, images, contrasts = [], [], []
    for coast in COASTS:
        for kind, target in MOMENTS:
            record, image, contrast = scene(coast, kind, target, calendar, slopes, series, fields, sky, earth_light,
                                            water, n_water, hours, jd)
            records.append(record); images.append(image); contrasts.append(contrast)
            print(coast, kind, "sun", record["sun"], "earth", record["earth"]["elevation_deg"],
                  "ground", record["light_on_ground_lux"], "glitter peak", record["glitter"]["peak_cd_m2"], flush=True)
    PANORAMAS.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PANORAMAS, xyz=np.stack(images), wave_contrast=np.stack(contrasts),
                        elevation_deg=ELEVATIONS, azimuth_deg=AZIMUTHS)
    files = [Path(__file__), ROOT / "illumination/water_surface/reflection.py",
             ROOT / "illumination/water_column/model.py", ROOT / "illumination/earthlight/model.py",
             ROOT / "research/studies/sea_appearance/lighting.py", ROOT / "illumination/sky/solved_transport.py"]
    own = files[:5]
    return dict(
        schema="terluna.research.sea-regimes/1",
        evidence=("The clear-sky radiance and colour of sky and sea in every direction for an eye 2 m above the water, "
                  "at six moments of the shore month at four coasts: the solved sky for the Sun and the Earth, rough-"
                  "water reflection with that hour's slopes, glitter of both disks, and the design-guess water's own "
                  "light. A calculation from models, without clouds, foam, polarization or the coastline."),
        reading_rule=("Panoramas (research/runs/sea_appearance/regimes.npz, not committed) hold CIE XYZ radiance in "
                      "cd/m2 on 1 by 0.5 degree cells, azimuth clockwise from north, elevation from the eye's "
                      "horizontal; cells below the horizon's dip are sea. Readouts give luminance and chromaticity "
                      "toward the brighter source (the Sun until 40 degrees down, then the Earth where it is up) and "
                      "away from it, at 30, 5 and 1 degrees above and 1, 5 and 20 below. wave_contrast is the spread "
                      "of reflected sky radiance over the facets divided by its mean."),
        producer=dict(lane="research", runner=str(Path(__file__).relative_to(ROOT)),
                      files={str(p.relative_to(ROOT)): sha256(p) for p in files}, constants=constants_used(own)),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (CALENDAR, SLOPES, SLOPE_SERIES, WATERS, SOLUTION)},
        panoramas=dict(path=str(PANORAMAS.relative_to(ROOT)), sha256=sha256(PANORAMAS), shape=list(np.stack(images).shape)),
        eye_height_m=EYE_HEIGHT_M, water_index=round(n_water, 5), facet_quadrature_order=ORDER,
        scenes=records)


def main():
    product = build()
    RESULT.write_text(json.dumps(product, indent=1, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
