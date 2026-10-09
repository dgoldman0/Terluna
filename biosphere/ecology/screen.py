"""A first screen of what the Open Moon's committed products already say about its ecology.

    python -m biosphere.ecology.screen
    # -> biosphere/ecology/results/first_screen.json

Three results, each read from committed products with the literature values named in sources.json:
- **Ultraviolet.** The ground's light under the titania stack (illumination/surface_light, the design Moon against
  its Earth control, clear sky over ground of albedo 0.1) weighted by biological action spectra: skin erythema
  (the UV index), the band that makes vitamin D, sunlight disinfection, and the visual pigments of bees and birds
  (Govardovskii et al. 2000 templates).
- **Settling.** How fast pollen, spores and bacteria fall at lunar gravity in the design air, and how long the day's
  mixed layer holds them, against Earth.
- **Trace gases.** Under the titania stack the air holds no OH and no O(1D) (atmosphere/middle_atmosphere), so
  methane leaves only into soils and N2O has no photochemical sink. Where methane settles for given lake and
  wetland emissions, how fast N2O builds from Earth-like soils and seas, at Earth's rates per area.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import (AVOGADRO, GAS_CONSTANT, MOON_RADIUS, MOON_SURFACE_GRAVITY, PLANCK, SPEED_OF_LIGHT,
                              STANDARD_GRAVITY)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.biosphere.ecology-first-screen/1'
SURFACE = ROOT / 'illumination' / 'surface_light' / 'results' / 'surface_light.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
AIR = ROOT / 'research' / 'baselines' / 'feasibility' / 'reference.json'
OUT = HERE / 'results' / 'first_screen.json'

ALBEDO = 0.1
SUN_DEG = ('90.0', '30.0')
RECEPTORS_NM = {                           # visual pigment peaks (sources.json)
    'bee_uv': 344.0, 'bee_blue': 436.0, 'bee_green': 544.0,        # honeybee, Peitsch et al. 1992
    'bird_uvs': 370.0, 'bird_vs': 410.0,                           # the two kinds of bird short-wave cone
}
DISINFECTION_DECADE_NM = (30.0, 50.0)      # inactivation falling a decade per 30-50 nm above 300 nm
PREVITAMIN_D_END_NM = 330.0                # the CIE 174:2006 previtamin D3 action spectrum ends at 330 nm

# Settling. Particles are spheres with Cunningham slip; the day's mixed layer holds them until they settle out.
PARTICLES = {                              # diameter m, density kg/m3
    'pollen_25um': (25e-6, 1000.0), 'spore_10um': (10e-6, 1100.0), 'spore_5um': (5e-6, 1100.0),
    'bacterium_1um': (1e-6, 1100.0),
}
MIXED_LAYER_M = {'moon_wet_land_day': 5000.0, 'moon_cm1_noon': 12000.0, 'earth_day': 1500.0}
WIND_M_S = 3.0                             # a light wind for the distance the residence carries a particle
EARTH_AIR = (101325.0, 288.15)             # Pa, K
AIR_MOLAR_MASS = 0.02897                   # kg/mol, for Earth; the Moon's comes from the reference inventory

# Trace gases at Earth's rates per area (sources.json). Ranges are low-high.
CH4_LAKE_G_M2_YR = (2.0, 30.0)
CH4_WETLAND_G_M2_YR = (10.0, 50.0)
WETLAND_SHARE_OF_MOON = (0.013, 0.053)     # a scenario: 0.5-2 million km2 of wet margins, no product yet
CENTRAL = dict(lake_g_m2_yr=10.0, wetland_g_m2_yr=20.0, wetland_share=0.026)   # 1 million km2 of wet margins
EARTH_SOIL_CH4_SINK_TG_YR = (11.0, 30.0, 49.0)
EARTH_UPLAND_SOIL_KM2 = 1.3e8
EARTH_CH4_PPM = 1.85
N2O_SOIL_G_N_M2_YR = (0.038, 0.050)        # natural soils, 4.9-6.5 Tg N/yr over about 1.3e8 km2
N2O_SEA_G_N_M2_YR = (0.007, 0.012)         # oceans, 2.5-4.3 Tg N/yr over 3.6e8 km2
EARTH_N2O_LIFETIME_YR = 116.0
GAP_TRANSMISSION = (3e-5, 2.6e-3)          # the swarm's UV transmission by design level (protection_architecture)

EVIDENCE = ('A screen on committed products: clear-sky surface spectra, the atlas and drainage areas and the '
            'reference air inventory, with Earth rates per area and published action spectra. No ecosystem, soil '
            'or chemistry model is run; the trace-gas rates are Earth analogues and the wetland area a scenario.')
READING_RULE = ('ultraviolet.by_sun[sun_deg] gives Moon/Earth ratios of weighted photon or energy fluxes on a '
                'horizontal surface over ground of albedo 0.1; receptor catches use photon flux. settling gives '
                'fall speeds in m/s and residence = mixed layer / fall speed in days. trace_gases gives Tg a year, '
                'ppm (or ppb) by mole of the whole air, and years.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def govardovskii(wavelength_nm, peak_nm):
    """A1 visual pigment template (Govardovskii et al. 2000): alpha band plus beta band, peak normalised near 1."""
    x = peak_nm / wavelength_nm
    a = 0.8795 + 0.0459 * np.exp(-(peak_nm - 300.0) ** 2 / 11940.0)
    alpha = 1.0 / (np.exp(69.7 * (a - x)) + np.exp(28.0 * (0.922 - x)) + np.exp(-14.9 * (1.104 - x)) + 0.674)
    beta_peak, beta_width = 189.0 + 0.315 * peak_nm, -40.5 + 0.195 * peak_nm
    return alpha + 0.26 * np.exp(-((wavelength_nm - beta_peak) / beta_width) ** 2)


def erythema(wavelength_nm):
    """CIE erythema action spectrum (ISO 17166); 40 m2/W times the weighted irradiance is the UV index."""
    w = wavelength_nm
    return np.where(w <= 298, 1.0, np.where(w <= 328, 10 ** (0.094 * (298 - w)),
                                            np.where(w <= 400, 10 ** (0.015 * (140 - w)), 0.0)))


def ultraviolet(product):
    w = np.array(product['centres_nm'])
    photons = w * 1e-9 / (PLANCK * SPEED_OF_LIGHT) / AVOGADRO * 1e6          # umol s-1 per W, by bin

    def down(case, sun):
        c = product['cases'][case]
        s = c['spectra'][sun]
        r = np.array(product['columns'][c['column']]['reflectance_from_below'])
        return (np.array(s['direct_w_m2_nm']) + np.array(s['diffuse_black_w_m2_nm'])) / (1 - ALBEDO * r)

    band = lambda lo, hi: (w >= lo) & (w < hi)
    out = {}
    for sun in SUN_DEG:
        moon, earth = down('moon_1.2atm/design', sun), down('earth_control/unfiltered', sun)
        pm, pe = moon * photons, earth * photons
        ratio = lambda weights, m=pm, e=pe: float(np.sum(m * weights) / np.sum(e * weights))
        row = dict(
            uv_index=dict(moon=round(float(40 * np.sum(moon * erythema(w))), 3), earth=round(float(40 * np.sum(earth * erythema(w))), 2)),
            uv_b_photons_umol_m2_s=dict(moon=round(float(np.sum(pm[band(280, 315)])), 6), earth=round(float(np.sum(pe[band(280, 315)])), 3)),
            previtamin_d_band_photons_ratio=float(f'{ratio(band(280, PREVITAMIN_D_END_NM)):.2g}'),
            uv_a_photons_ratio=round(ratio(band(315, 400)), 3),
            uv_a_to_par_photons=dict(moon=round(float(np.sum(pm[band(315, 400)]) / np.sum(pm[band(400, 700)])), 4),
                                     earth=round(float(np.sum(pe[band(315, 400)]) / np.sum(pe[band(400, 700)])), 4)),
            disinfection_ratio={f'decade_per_{d:g}_nm': round(ratio(10 ** (-(w - 300) / d) * band(280, 450)), 4)
                                for d in DISINFECTION_DECADE_NM},
            receptor_catch_ratio={k: round(ratio(govardovskii(w, nm)), 3) for k, nm in RECEPTORS_NM.items()},
        )
        row['bee_uv_to_green_ratio'] = round(row['receptor_catch_ratio']['bee_uv'] / row['receptor_catch_ratio']['bee_green'], 3)
        out[sun] = row
    return dict(albedo=ALBEDO, receptors_nm=RECEPTORS_NM, by_sun=out)


def viscosity(temperature_k):
    """Air viscosity by Sutherland's law (Pa s); it does not depend on pressure."""
    return 1.716e-5 * (temperature_k / 273.15) ** 1.5 * (273.15 + 110.4) / (temperature_k + 110.4)


def fall_speed(diameter, density, gravity, pressure, temperature, molar_mass):
    mu = viscosity(temperature)
    rho_air = pressure * molar_mass / (GAS_CONSTANT * temperature)
    free_path = mu / pressure * np.sqrt(np.pi * GAS_CONSTANT * temperature / (2 * molar_mass))
    kn = 2 * free_path / diameter
    slip = 1 + kn * (1.257 + 0.4 * np.exp(-1.1 / kn))
    return (density - rho_air) * gravity * diameter ** 2 * slip / (18 * mu), rho_air


def settling(column, molar_mass):
    p, t = column['surface_pressure_pa'], column['surface_temperature_k']
    rows = {}
    for name, (d, rho) in PARTICLES.items():
        v_moon, rho_moon = fall_speed(d, rho, MOON_SURFACE_GRAVITY, p, t, molar_mass)
        v_earth, rho_earth = fall_speed(d, rho, STANDARD_GRAVITY, *EARTH_AIR, AIR_MOLAR_MASS)
        days = {k: round(h / (v_moon if k.startswith('moon') else v_earth) / 86400, 1) for k, h in MIXED_LAYER_M.items()}
        rows[name] = dict(fall_m_s=dict(moon=float(f'{v_moon:.3g}'), earth=float(f'{v_earth:.3g}')),
                          residence_days=days,
                          carried_km_in_moon_wet_land_day=round(WIND_M_S * days['moon_wet_land_day'] * 86.4, 0))
    plumed = float(np.sqrt(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY * rho_earth / rho_moon))
    return dict(air=dict(pressure_pa=p, temperature_k=t, density_kg_m3=round(float(rho_moon), 3),
                         viscosity_pa_s=float(f'{viscosity(t):.4g}')),
                mixed_layer_m=MIXED_LAYER_M, wind_m_s=WIND_M_S, particles=rows,
                plumed_seed_fall_ratio=round(plumed, 3), plumed_seed_reach_ratio=round(1 / plumed, 2))


def trace_gases(air, atlas, drainage):
    area = 4 * np.pi * MOON_RADIUS ** 2
    air_mol = air['mass_kg'] / air['mean_molar_mass']
    lakes = drainage['lakes']['area_share']
    land = atlas['main_land_share'] + atlas['island_share']
    soil = land - lakes                                    # dry land, the soils that take methane up
    tg = lambda g_m2, share: g_m2 * share * area / 1e12
    lake_tg = [tg(r, lakes) for r in CH4_LAKE_G_M2_YR]
    wet_tg = [tg(r, s) for r, s in zip(CH4_WETLAND_G_M2_YR, WETLAND_SHARE_OF_MOON)]
    sources = [lake_tg[0] + wet_tg[0], lake_tg[1] + wet_tg[1]]
    per_ppm_tg = air_mol * 1e-6 * 16.04 / 1e12            # Tg of CH4 in 1 ppm of the whole air
    sink_per_ppm = [s / EARTH_CH4_PPM / EARTH_UPLAND_SOIL_KM2 * soil * area / 1e6 for s in EARTH_SOIL_CH4_SINK_TG_YR]
    steady = [sources[0] / sink_per_ppm[2], sources[1] / sink_per_ppm[0]]
    lifetime = [per_ppm_tg / k for k in sink_per_ppm[::-1]]
    central_tg = tg(CENTRAL['lake_g_m2_yr'], lakes) + tg(CENTRAL['wetland_g_m2_yr'], CENTRAL['wetland_share'])
    central = dict(inputs=CENTRAL, sources_tg_yr=round(central_tg, 1), steady_ppm=round(central_tg / sink_per_ppm[1], 1),
                   lifetime_years=round(per_ppm_tg / sink_per_ppm[1], -1))
    methane = dict(lake_share=lakes, dry_land_share=round(soil, 4), wetland_share_scenario=WETLAND_SHARE_OF_MOON,
                   sources_tg_yr=dict(lakes=[round(x, 1) for x in lake_tg], wetlands=[round(x, 1) for x in wet_tg],
                                      total=[round(x, 1) for x in sources]),
                   tg_per_ppm=round(per_ppm_tg, 0), soil_sink_tg_yr_per_ppm=[round(x, 2) for x in sink_per_ppm],
                   lifetime_years=[round(x, -1) for x in lifetime],
                   steady_ppm=[round(x, 1) for x in steady],
                   rise_ppm_per_year_without_methanotrophs=[round(s / per_ppm_tg, 4) for s in sources],
                   central=central, earth_ppm=EARTH_CH4_PPM)
    sea = atlas['water_share']
    n2o_tg_n = [(soil * a + sea * b) * area / 1e12 for a, b in zip(N2O_SOIL_G_N_M2_YR, N2O_SEA_G_N_M2_YR)]
    per_ppb_tg_n = air_mol * 1e-9 * 28.014 / 1e12
    rise = [x / per_ppb_tg_n for x in n2o_tg_n]
    lifetime_gaps = [EARTH_N2O_LIFETIME_YR / f for f in GAP_TRANSMISSION[::-1]]
    nitrous = dict(sources_tg_n_yr=[round(x, 2) for x in n2o_tg_n], tg_n_per_ppb=round(per_ppb_tg_n, 2),
                   rise_ppb_per_year=[round(x, 2) for x in rise],
                   ppm_after_500_years=[round(x * 500 / 1000, 2) for x in rise],
                   ppm_after_10000_years=[round(x * 10000 / 1000, 1) for x in rise],
                   lifetime_years_from_gap_light=[float(f'{x:.2g}') for x in lifetime_gaps],
                   note='No light below 320 nm reaches the air, so N2O is neither photolysed nor attacked by O(1D); '
                        'the swarm\'s gaps pass 3e-5 to 2.6e-3 of sunlight, which scales Earth\'s 116-year lifetime.')
    return dict(air_mol=float(f'{air_mol:.4g}'), methane=methane, nitrous_oxide=nitrous)


def main(argv=None) -> int:
    surface, atlas, drainage, air = (json.loads(p.read_text()) for p in (SURFACE, ATLAS, DRAINAGE, AIR))
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={'biosphere/ecology/screen.py': digest(__file__)},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in (SURFACE, ATLAS, DRAINAGE, AIR)},
                      constants=constants_used(['biosphere/ecology/screen.py'])),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        ultraviolet=ultraviolet(surface),
        settling=settling(surface['columns']['moon_1.2atm'], air['inventory']['mean_molar_mass']),
        trace_gases=trace_gases(air['inventory'], atlas, drainage))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(json.dumps({k: product[k] for k in ('ultraviolet', 'settling', 'trace_gases')}, indent=1)[:6000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
