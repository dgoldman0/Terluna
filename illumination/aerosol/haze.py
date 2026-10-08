"""Haze in the Open Moon's sky: the aerosol study's particles as optics, and what they do to visibility, daylight and the
long twilight.

    python -m illumination.aerosol.haze
    # -> illumination/aerosol/results/haze.json

The Open-Moon aerosol study (research/studies/open_moon_aerosol, atmospheric-electricity branch) estimates lognormal
particle modes near the ground for five kinds of region, clean, central and loaded, by day and by night, with each
mode's hygroscopicity. This module turns them into optics at 550 nm: the modes grow with the region's humidity by
kappa-Koehler theory (Petters and Kreidenweis 2007, without the Kelvin term, humidity capped at 95%), mix with water
by volume, and scatter by Mie theory (mie.py). Dry refractive indices follow OPAC (Hess, Koepke and Schult 1998) for
water-soluble, sea-salt and mineral particles; the biological modes take a weakly absorbing 1.53 + 0.003i, an
assumption. The haze fills the region's mixed layer, falling off above the ground with that depth as its scale
height. Visibility is the meteorological optical range, 3.912 over the extinction at the ground, molecules included.
The column's daylight comes from the domain's delta two-stream solver (atmosphere/radiative_convective/shortwave.py)
with the Open Moon's molecular optical depth at 550 nm; the twilight arrives as diffuse light above the haze, which
the haze layer passes with its diffuse transmittance.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from atmosphere.radiative_convective.shortwave import column_fluxes, layer_properties
from illumination.aerosol.mie import lognormal_mode
from shared.constants import MOON_GM, MOON_RADIUS, AVOGADRO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.illumination.haze/1'
AEROSOL = ROOT / 'research' / 'studies' / 'open_moon_aerosol' / 'results' / 'open_moon_aerosol.json'
CALENDAR = ROOT / 'illumination' / 'calendar' / 'results' / 'transfer.json'
OUT = HERE / 'results' / 'haze.json'
WAVELENGTH_UM = 0.55
INDEX = {'aitken': 1.53 + 0.006j, 'accumulation': 1.53 + 0.006j, 'spray': 1.50 + 1e-8j,
         'dust': 1.53 + 0.0055j, 'biological': 1.53 + 0.003j, 'fungal': 1.53 + 0.003j}
WATER = 1.333 + 1e-9j
RH_CAP = 0.95
SURFACE_PRESSURE_PA = 121590.0        # the design atmosphere, 1.2 atm (research/decisions.md)
SURFACE_K = 290.0
EARTH_RAYLEIGH_TAU_550 = 0.0973       # Earth's standard column at 550 nm (Bodhaine et al. 1999)
EARTH_COLUMN_KG_M2 = 101325.0 / 9.80665
EARTH_RAYLEIGH_PER_KM_550 = 0.01162   # at 1013.25 hPa and 288.15 K
ALBEDO = 0.15
DUSK_LUX = 2.98
DEG_PER_HOUR = 360.0 / (29.530589 * 24.0)
EVIDENCE = ('Mie optics of the aerosol study\'s ground-level modes (Earth analogues changed by lunar factors) with '
            'assumed refractive indices, a mixed-layer profile and kappa-Koehler growth; one wavelength; the haze '
            'aloft and fog are outside it. Daylight from a two-layer delta two-stream column; twilight as clear-sky '
            'diffuse light passed through the haze layer.')
READING_RULE = ('Coefficients per km at 550 nm. aod is the haze\'s optical depth through its mixed-layer profile; '
                'visibility_km includes the air\'s own scattering. daylight gives the share of sunlight reaching level '
                'ground as direct and diffuse at three Sun zenith angles, clear and hazy. twilight gives the '
                'equatorial hours after sunset until 2.98 lux, clear and hazy.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def mode_optics(mode, rh):
    rh = min(rh, RH_CAP)
    growth = (1.0 + mode['kappa'] * rh / (1.0 - rh)) ** (1.0 / 3.0)
    dry = INDEX.get(mode['name'], 1.53 + 0.003j)
    wet = (dry + (growth ** 3 - 1.0) * WATER) / growth ** 3
    return lognormal_mode(mode['number_cm3'], mode['d_nm'] * 1e-3 * growth, mode['sigma'], wet, WAVELENGTH_UM)


def haze(modes, rh):
    ext = sca = gsc = 0.0
    for m in modes:
        o = mode_optics(m, rh)
        ext += o['extinction_per_km']; sca += o['scattering_per_km']; gsc += o['scattering_per_km'] * o['asymmetry']
    return dict(extinction_per_km=ext, single_scattering_albedo=sca / ext if ext else 1.0, asymmetry=gsc / sca if sca else 0.0)


def molecular_tau():
    gravity = MOON_GM / MOON_RADIUS ** 2
    return EARTH_RAYLEIGH_TAU_550 * (SURFACE_PRESSURE_PA / gravity) / EARTH_COLUMN_KG_M2


def daylight(aod, ssa, g, mixed_km, tau_mol):
    """Shares of top-of-air sunlight reaching the ground, direct and diffuse, at three Sun zenith angles."""
    scale_km = 8.314462618 * SURFACE_K / (0.02897 * MOON_GM / MOON_RADIUS ** 2) / 1e3
    low = tau_mol * (1.0 - np.exp(-mixed_km / scale_km))
    out = {}
    for zenith in (0.0, 60.0, 80.0):
        mu0 = np.cos(np.radians(zenith))
        rows = {}
        for label, a in (('clear', 0.0), ('hazy', aod)):
            tau_sca = np.array([[tau_mol - low], [low + a * ssa]])
            tau_abs = np.array([[0.0], [a * (1.0 - ssa)]])
            gg = np.array([[0.0], [(a * ssa * g) / max(low + a * ssa, 1e-12)]])
            direct, diffuse, _ = column_fluxes(tau_abs, tau_sca, gg, ALBEDO, mu0)
            rows[label] = dict(direct=round(float(direct[-1, 0] / mu0), 4), diffuse=round(float(diffuse[-1, 0] / mu0), 4))
        out[f'{zenith:g}'] = rows
    return out


def twilight_hours(aod, ssa, g, calendar):
    e = np.array(calendar['diffuse_elevation_deg']); lux = np.array(calendar['solar']['diffuse_lux'])
    _, t_dif, _, _, _ = layer_properties(np.array([[aod]]), np.array([[ssa]]), np.array([[g]]), np.array([[0.5]]))
    passed = float(t_dif[0, 0])
    order = np.argsort(lux)
    clear = -np.interp(DUSK_LUX, lux[order], e[order])
    hazy = -np.interp(DUSK_LUX / passed, lux[order], e[order])
    return dict(diffuse_transmittance=round(passed, 4), clear_depression_deg=round(clear, 2), hazy_depression_deg=round(hazy, 2),
                clear_hours=round(clear / DEG_PER_HOUR, 1), hazy_hours=round(hazy / DEG_PER_HOUR, 1))


def main(argv=None) -> int:
    aero, calendar = json.loads(AEROSOL.read_text()), json.loads(CALENDAR.read_text())
    tau_mol = molecular_tau()
    beta_mol = EARTH_RAYLEIGH_PER_KM_550 * (SURFACE_PRESSURE_PA / 101325.0) * (288.15 / SURFACE_K)
    regions = {}
    for name, region in aero['regions'].items():
        regions[name] = dict(area_share=region['area_share'], cases={})
        for case, by_time in region['cases'].items():
            regions[name]['cases'][case] = {}
            for time in ('day', 'night'):
                w = region['weather'][time]
                h = haze(by_time[time]['modes'], w['relative_humidity'])
                aod = h['extinction_per_km'] * w['mixed_layer_km']
                row = dict(relative_humidity=round(w['relative_humidity'], 3), mixed_layer_km=round(w['mixed_layer_km'], 2),
                           extinction_per_km=round(h['extinction_per_km'], 5), single_scattering_albedo=round(h['single_scattering_albedo'], 4),
                           asymmetry=round(h['asymmetry'], 3), aod=round(aod, 4),
                           visibility_km=round(3.912 / (h['extinction_per_km'] + beta_mol), 1))
                if time == 'day':
                    row['daylight'] = daylight(aod, h['single_scattering_albedo'], h['asymmetry'], w['mixed_layer_km'], tau_mol)
                else:
                    row['twilight'] = twilight_hours(aod, h['single_scattering_albedo'], h['asymmetry'], calendar)
                regions[name]['cases'][case][time] = row
    product = dict(schema=SCHEMA, producer=dict(domain='illumination', files={'aerosol/haze.py': digest(__file__),
                                                                               'aerosol/mie.py': digest(HERE / 'mie.py')},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (AEROSOL, CALENDAR)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, wavelength_um=WAVELENGTH_UM,
                   refractive_index={k: [v.real, v.imag] for k, v in INDEX.items()},
                   molecular_tau_550=round(tau_mol, 3), molecular_extinction_per_km_ground=round(beta_mol, 4),
                   clear_visibility_km=round(3.912 / beta_mol, 0), regions=regions)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(f"molecular optical depth at 550 nm {tau_mol:.3f}; clear-air visibility {3.912 / beta_mol:.0f} km")
    for name, r in regions.items():
        for case, c in r['cases'].items():
            d, n = c['day'], c['night']
            print(f"{name:15s} {case:8s} day AOD {d['aod']:.3f} (ssa {d['single_scattering_albedo']:.3f}, g {d['asymmetry']:.2f}) "
                  f"vis {d['visibility_km']:6.1f} km | noon direct {d['daylight']['0']['clear']['direct']:.3f}->{d['daylight']['0']['hazy']['direct']:.3f} "
                  f"diffuse {d['daylight']['0']['clear']['diffuse']:.3f}->{d['daylight']['0']['hazy']['diffuse']:.3f} | "
                  f"night AOD {n['aod']:.3f}, dusk {n['twilight']['clear_hours']}->{n['twilight']['hazy_hours']} h")
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
