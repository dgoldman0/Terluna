"""The solver's atmospheres as three-wavelength parameters for engine skies.

Real-time skies such as Unreal Engine's SkyAtmosphere describe an atmosphere with
coefficients at three wavelengths; its Earth defaults are Bruneton's molecular and
ozone values at 680, 550 and 440 nm, which the spectral solver also uses. This
module derives the same parameters for each of the solver's atmospheres from
atmospheres.py, so that an engine sky is configured from the illumination domain
rather than by hand. It also reports how far such a sky's direct sunlight departs
from the spectral calculation, which grows toward the horizon.

Run: python illumination/sky/engine_atmosphere.py [out.json]
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np

try:
    from .atmospheres import SW, SF, EARTH, MOON, MOON_ZERO, Atmosphere, rayleigh, ozone_absorption
    from .colour_matching import colour_weights
except ImportError:  # run as a script from this folder
    from atmospheres import SW, SF, EARTH, MOON, MOON_ZERO, Atmosphere, rayleigh, ozone_absorption
    from colour_matching import colour_weights

SCHEMA = 'terluna.illumination.engine-atmosphere/1'
ENGINE_WAVELENGTHS_NM = (680., 550., 440.)
LUMINANCE = np.array([0.2126, 0.7152, 0.0722])
HERE = Path(__file__).resolve().parent
SPECTRUM = np.arange(380., 800.01, 5.)
CHECK_ELEVATIONS_DEG = (90, 45, 20, 10, 5, 2, 0)
FIT_ELEVATIONS_DEG = (90, 60, 45, 30, 20, 15, 10, 7, 5, 3, 2)


def sun_rgb(atm: Atmosphere):
    """Top-of-atmosphere sunlight as photopic-weighted linear sRGB (lux per channel)."""
    return (np.interp(SPECTRUM, SW, SF) * atm.visible_filter) @ colour_weights(SPECTRUM)


def spectral_direct(atm: Atmosphere, elevation_deg: float):
    """Direct sunlight at the ground by the spectral calculation (linear sRGB lux)."""
    gas, ozone = paths(atm, elevation_deg)
    solar = np.interp(SPECTRUM, SW, SF) * atm.visible_filter
    tr = np.exp(-rayleigh(SPECTRUM, atm) * gas - ozone_absorption(SPECTRUM, atm) * ozone)
    return (solar * tr) @ colour_weights(SPECTRUM), gas, ozone


def fitted_rayleigh(atm: Atmosphere):
    """Per-channel molecular coefficients (1/m) that best reproduce the spectral direct
    sunlight over Sun elevations (least squares in log transmission), keeping the sampled
    ozone. Channels where the spectral sunlight leaves the sRGB gamut are not fitted there."""
    lam, sun = np.array(ENGINE_WAVELENGTHS_NM), sun_rgb(atm)
    a = ozone_absorption(lam, atm)
    k = np.zeros(3)
    for c in range(3):
        num = den = 0.
        for e in FIT_ELEVATIONS_DEG:
            rgb, gas, ozone = spectral_direct(atm, e)
            if rgb[c] <= 1e-6 * sun[c]:
                continue
            y = np.log(rgb[c] / sun[c]) + a[c] * ozone
            num += gas * y
            den += gas * gas
        k[c] = -num / den
    return k


def engine_parameters(atm: Atmosphere) -> dict:
    lam = np.array(ENGINE_WAVELENGTHS_NM)
    sun = sun_rgb(atm)
    return {
        'bottom_radius_km': atm.radius_m / 1000.,
        'atmosphere_height_km': atm.scale_height_m * atm.top_scale_heights / 1000.,
        # Sampled at the three wavelengths, as engine defaults are for Earth.
        'rayleigh_scattering_per_km': (rayleigh(lam, atm) * 1000.).tolist(),
        # Fitted to the spectral direct sunlight over Sun elevations (see direct_sun_check).
        'rayleigh_scattering_fitted_per_km': (fitted_rayleigh(atm) * 1000.).tolist(),
        'rayleigh_scale_height_km': atm.scale_height_m / 1000.,
        # The solver's atmosphere is molecular only: no aerosols.
        'mie_scattering_per_km': [0., 0., 0.],
        'mie_absorption_per_km': [0., 0., 0.],
        # Ozone: a tent from 10 to 40 km, scaled with the scale height, peak coefficient.
        'absorption_per_km': (ozone_absorption(lam, atm) * 1000.).tolist(),
        'absorption_tent': {'tip_altitude_km': 25. * atm.ozone_scale, 'width_km': 15. * atm.ozone_scale},
        'ground_albedo': [atm.ground_albedo] * 3,
        'sun_illuminance_lux': float(sun @ LUMINANCE),
        'sun_colour_linear_srgb': (sun / sun.max()).tolist(),
        'sun_angular_diameter_deg': 2. * atm.sun_radius_deg,
    }


def paths(atm: Atmosphere, elevation_deg: float, n: int = 40000):
    """Gas and ozone columns (m, relative to surface density) from the ground to the top."""
    R, H = atm.radius_m, atm.scale_height_m
    e = np.radians(elevation_deg)
    top = atm.top - R
    length = -R * np.sin(e) + np.sqrt((R * np.sin(e)) ** 2 + 2 * R * top + top * top)
    s = (np.arange(n) + .5) / n * length
    h = np.sqrt(R * R + s * s + 2 * R * s * np.sin(e)) - R
    q = h / atm.ozone_scale
    ozone = np.clip(np.minimum((q - 10000.) / 15000., (40000. - q) / 15000.), 0., None)
    ds = length / n
    return float(np.exp(-h / H).sum() * ds), float(ozone.sum() * ds)


def direct_sun_check(atm: Atmosphere) -> list:
    """Direct sunlight at the ground: spectral calculation vs three-wavelength models
    with sampled and with fitted molecular coefficients."""
    lam = np.array(ENGINE_WAVELENGTHS_NM)
    sun, fitted, a = sun_rgb(atm), fitted_rayleigh(atm), ozone_absorption(lam, atm)
    rows = []
    for e in CHECK_ELEVATIONS_DEG:
        spectral, gas, ozone = spectral_direct(atm, e)
        sampled = sun * np.exp(-rayleigh(lam, atm) * gas - a * ozone)
        fit = sun * np.exp(-fitted * gas - a * ozone)
        rows.append({
            'sun_elevation_deg': e,
            'spectral_rgb_lux': spectral.tolist(),
            'sampled_rgb_lux': sampled.tolist(),
            'fitted_rgb_lux': fit.tolist(),
            'sampled_luminance_ratio': float((sampled @ LUMINANCE) / (spectral @ LUMINANCE)),
            'fitted_luminance_ratio': float((fit @ LUMINANCE) / (spectral @ LUMINANCE)),
        })
    return rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def product() -> dict:
    return {
        'schema': SCHEMA,
        'producer': {
            'domain': 'illumination',
            'model': 'illumination/sky/atmospheres.py',
            'model_sha256': sha256(HERE / 'atmospheres.py'),
            'colour': 'illumination/sky/colour_matching.py',
            'colour_sha256': sha256(HERE / 'colour_matching.py'),
            'exporter': 'illumination/sky/engine_atmosphere.py',
        },
        'evidence': ('The spectral solver\'s atmospheres sampled at three wavelengths, as engine '
                     'skies sample Earth\'s. For configuring a display sky; the spectral atlas '
                     'remains the reference, and direct_sun_check shows the three-wavelength '
                     'error for sunlight, which grows toward the horizon.'),
        'units': 'lengths in km; coefficients in 1/km; illuminance in lux; colours linear sRGB',
        'wavelengths_nm': list(ENGINE_WAVELENGTHS_NM),
        'atmospheres': {a.name: engine_parameters(a) for a in (MOON, MOON_ZERO, EARTH)},
        'direct_sun_check': {a.name: direct_sun_check(a) for a in (MOON, EARTH)},
    }


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / 'products' / 'engine_atmosphere.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(product(), indent=1))
    print(out)
