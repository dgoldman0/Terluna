#!/usr/bin/env python3
"""The design air's hydrogen escape with its methane: what diffusion through the homopause can carry up, for each
carrier and for methane broken or unbroken there, checked against the oxygen step's full transport for H2 and water.

    python -m atmosphere.loss_response.hydrogen        # seconds: results/hydrogen_escape.json

Above the homopause each gas spreads with its own scale height, so a light gas's escape is limited by how fast it
diffuses up through the heavier air (Hunten 1973): b f (m_a - m) g / (k T) molecules per square metre and second, with
b = n D the gas's binary diffusion parameter in air, f its share of the air below the homopause, m_a the air's mean
molecular mass, m the gas's and T the temperature there. Over a sphere g r^2 = G M, so the Moon loses
4 pi G M b f (m_a - m) / (k T) molecules a second wherever its homopause lies. At its exobase hydrogen's Jeans
parameter is below 2, so what reaches it leaves, and the diffusion sets the escape.

The carriers: H2 at the middle atmosphere's 0.53 ppm (Earth's present air); water at the share the cold trap lets past
the tropopause (the middle atmosphere's design cases, at grounds of 288-308 K, and two and a half times that with the
cold cloud at lunar gravity, which warms the GCM's coldest level 3.7 K); and the design air's 1.9 ppm of methane
(register, 13 October), which the middle atmosphere's chemistry leaves out. Unbroken at the homopause, methane
diffuses as methane, slowly and with 0.56 of the air's mass; broken below it, its four hydrogen atoms arrive as H2 and
H and diffuse as H2. Behind the titania stack the sky's Lyman-alpha glow, which no Sun-facing shield blocks, breaks
it above the 0.3 Pa base, where O2 lets the glow through. A steady tracer mixed upward at the eddy diffusion K and
broken at the glow's rate J falls off as exp(-z / sqrt(K / J)), so the share broken below the homopause is
1 - exp(-dz / sqrt(K / J)), dz the height from the base to methane's homopause, where K equals its molecular
diffusivity. K above the base is the middle atmosphere's mixing for the Moon at its cap.

Diffusion parameters: Banks and Kockarts (1973) for H2, as the oxygen step takes them; Massman (1998) for methane and
water in air. The temperature at the homopause is the 0.3 Pa base's in the three design treatments with the curved
light (curved_light.json): behind the film alone the thermal column above it stays within 0.4 K of it. The oxygen
step's transport with chemistry (oxygen.py) gives the H2 and water's escape behind the stack, against which the
limiting flux is checked. Cycle times: the escape against the air's hydrogen in H2 and methane, against methane's
inventory, and as water against the waters the Moon holds and their delivery.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

from shared.constants import AVOGADRO, BOLTZMANN, MOON_GM, MOON_RADIUS
from atmosphere.radiative_convective import thermodynamics as th

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'hydrogen_escape.json'
SCHEMA = 'terluna.atmosphere.hydrogen-escape/1'
AMU_KG = 1.66053907e-27
MASS_AMU = dict(H=1.008, H2=2.016, H2O=18.015, CH4=16.043)
# Binary diffusion parameters b = n D = A (T / T0)^s, cm^-1 s^-1: Banks and Kockarts (1973) for H2 in N2, as the
# oxygen step has it (A T^s with T0 = 1 K); Massman (1998) for CH4 and H2O in air, D at 273.15 K and 101.325 kPa
# (0.1952 and 0.2178 cm^2/s) times the Loschmidt number, with his exponent 1.81 on D (0.81 on b).
LOSCHMIDT_CM3 = 2.6867811e19
DIFFUSION = dict(H2=(2.73e17, 1.0, 0.741), CH4=(0.1952 * LOSCHMIDT_CM3, 273.15, 0.81),
                 H2O=(0.2178 * LOSCHMIDT_CM3, 273.15, 0.81))
H2_PPM = 0.53                       # the middle atmosphere's H2 (equilibrium.Case.h2_ppm)
METHANE_PPM = 1.9                   # the design air's methane (register, 13 October)
COLD_CLOUD_WATER = 2.5              # the cold cloud's factor on the water past the GCM's coldest level (A28R_301k_cold)
CH4_LYMAN_ALPHA_CM2 = 1.85e-17      # methane at 121.6 nm (Mount, Warden and Moos 1977)
EDDY_CAP_CM2_S = 1.0e7 * 36.0       # the middle atmosphere's mixing cap times the Moon's 36 (chemistry.Mixing)
BASE_PA = 0.3
AIR_KG = 2.8e18                     # the 1.2 atm atmosphere (middle_atmosphere/README.md, losses)
WATERS_KM3 = 13.8e6                 # the common-level waters and the rain-fed lakes (studies/sea_chemistry)
DELIVERY_KM3_YR = dict(first_thousand_years=(1500.0, 4100.0), to_10000_years=1700.0, to_100000_years=240.0,
                       by_a_million_years=10.0)   # the top-up the crust's uptake takes (geography/README.md)
YEAR_S = 3.15576e7
CURVED = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results' / 'curved_light.json'
CASES = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results' / 'cases'
OXYGEN = HERE / 'results' / 'oxygen_escape.json'
TREATMENTS = dict(lte='moon_1.2atm_titania_stack', all_heats='moon_1.2atm_titania_stack_nonlte',
                  collisional='moon_1.2atm_titania_stack_nonlte_nir_collisional')
GROUNDS = dict(ts288='moon_1.2atm_titania_stack', ts298='moon_1.2atm_titania_stack_ts298',
               ts308='moon_1.2atm_titania_stack_ts308')


def air_mass_amu():
    return th.earthlike_air(121590.0, 400.0).molar_mass * 1e3


def diffusion_parameter(gas, t):
    """b = n D (m^-1 s^-1) of a gas in air at temperature t."""
    a, t0, s = DIFFUSION[gas]
    return 100.0 * a * (t / t0) ** s


def limiting_rate(gas, share, t, m_air_amu=None):
    """Molecules a second the whole Moon loses of a gas at a share of the air below the homopause (Hunten 1973)."""
    m_air = (m_air_amu or air_mass_amu()) * AMU_KG
    return 4.0 * math.pi * MOON_GM * diffusion_parameter(gas, t) * share * (m_air - MASS_AMU[gas] * AMU_KG) / (
        BOLTZMANN * t)


def limiting_flux_per_m2(gas, share, t, radius_m, m_air_amu=None):
    """The same per square metre at a radius: b f (m_a - m) g / (k T)."""
    m_air = (m_air_amu or air_mass_amu()) * AMU_KG
    g = MOON_GM / radius_m ** 2
    return diffusion_parameter(gas, t) * share * (m_air - MASS_AMU[gas] * AMU_KG) * g / (BOLTZMANN * t)


def hydrogen_kg_s(rates):
    """kg/s of hydrogen in molecules a second of each gas."""
    atoms = dict(H2=2, H2O=2, CH4=4)
    return sum(atoms[g] * r for g, r in rates.items()) * MASS_AMU['H'] * AMU_KG


def glow_rate():
    """Methane's photolysis rate (s^-1) by the sky's Lyman-alpha glow where O2 lets it through: twice the flux onto
    the column, as the oxygen step's glow rates take the actinic flux of isotropic light from above."""
    from atmosphere.loss_response import model as lr
    intensity = lr.GLOW_RAYLEIGH * 1e10 / (4.0 * math.pi) * 1e-4      # photons cm^-2 s^-1 sr^-1
    return CH4_LYMAN_ALPHA_CM2 * 2.0 * math.pi * intensity


def broken_share(t, base_radius_m, m_air_amu=None):
    """The share of methane broken by the glow between the base and its homopause, and the heights behind it."""
    m_air = (m_air_amu or air_mass_amu()) * AMU_KG
    g = MOON_GM / base_radius_m ** 2
    scale_m = BOLTZMANN * t / (m_air * g)
    n_base = BASE_PA / (BOLTZMANN * t) * 1e-6                          # cm^-3
    d_base = diffusion_parameter('CH4', t) * 1e-2 / n_base            # cm^2/s
    n_homopause = diffusion_parameter('CH4', t) * 1e-2 / EDDY_CAP_CM2_S
    dz_m = scale_m * math.log(n_base / n_homopause)
    j = glow_rate()
    decay_m = math.sqrt(EDDY_CAP_CM2_S * 1e-4 / j)
    return dict(scale_height_km=scale_m / 1e3, methane_diffusivity_at_base_cm2_s=d_base,
                homopause_pa=n_homopause * 1e6 * BOLTZMANN * t, height_to_homopause_km=dz_m / 1e3,
                glow_rate_s=j, glow_lifetime_days=1.0 / j / 86400.0, decay_length_km=decay_m / 1e3,
                share=1.0 - math.exp(-dz_m / decay_m))


def escape_rates(t, water_share, broken):
    """Hydrogen escape (kg/s) by carrier with methane broken in the given share below the homopause."""
    m_air = air_mass_amu()
    rates = dict(H2=limiting_rate('H2', H2_PPM * 1e-6, t, m_air), H2O=limiting_rate('H2O', water_share, t, m_air))
    methane = dict(unbroken=dict(CH4=limiting_rate('CH4', METHANE_PPM * 1e-6, t, m_air)),
                   broken=dict(H2=limiting_rate('H2', 2.0 * METHANE_PPM * 1e-6, t, m_air)))
    kg = dict(h2=hydrogen_kg_s({'H2': rates['H2']}), water=hydrogen_kg_s({'H2O': rates['H2O']}),
              methane_unbroken=hydrogen_kg_s(methane['unbroken']), methane_broken=hydrogen_kg_s(methane['broken']))
    kg['methane'] = (1.0 - broken) * kg['methane_unbroken'] + broken * kg['methane_broken']
    kg['without_methane'] = kg['h2'] + kg['water']
    kg['total'] = kg['without_methane'] + kg['methane']
    kg['total_unbroken'] = kg['without_methane'] + kg['methane_unbroken']
    kg['total_broken'] = kg['without_methane'] + kg['methane_broken']
    return kg


def cycle_times(kg_s):
    """Years for an escape of kg_s of hydrogen to take the air's hydrogen in H2 and methane, methane's own, and, as
    water, the Moon's waters; and the water lost a year against the delivery."""
    m_air = air_mass_amu()
    h_free = AIR_KG * (H2_PPM * 1e-6 * 2 * MASS_AMU['H'] + METHANE_PPM * 1e-6 * 4 * MASS_AMU['H']) / m_air
    methane_kg = AIR_KG * METHANE_PPM * 1e-6 * MASS_AMU['CH4'] / m_air
    water_kg_s = kg_s * MASS_AMU['H2O'] / (2 * MASS_AMU['H'])
    water_km3_yr = water_kg_s * YEAR_S / 1e12
    return dict(air_hydrogen_kg=h_free, air_hydrogen_years=h_free / kg_s / YEAR_S,
                methane_kg=methane_kg, methane_hydrogen_years=methane_kg * 4 * MASS_AMU['H'] / MASS_AMU['CH4'] / kg_s / YEAR_S,
                water_equivalent_kg_s=water_kg_s, water_km3_per_year=water_km3_yr,
                waters_years=WATERS_KM3 / water_km3_yr,
                share_of_delivery={k: (water_km3_yr / v[1], water_km3_yr / v[0]) if isinstance(v, tuple)
                                   else water_km3_yr / v for k, v in DELIVERY_KM3_YR.items()})


def analyse():
    from atmosphere.loss_response import model as lr
    curved = json.loads(CURVED.read_text())['exobase']
    temps = {k: float(curved[c]['curved']['base_temperature_k']) for k, c in TREATMENTS.items()}
    water = {k: float(json.loads((CASES / f'{c}.json').read_text())['stratospheric_h2o']) for k, c in GROUNDS.items()}
    base_r = lr.base_radius_R('titania_stack', 'lte') * MOON_RADIUS
    oxygen = json.loads(OXYGEN.read_text())['summary']['moon_mixing']['titania_stack']['hydrogen_escape_kg_s']
    out = dict(schema=SCHEMA, air_molar_mass_amu=air_mass_amu(), homopause_temperatures_k=temps,
               water_past_the_tropopause=water, base_radius_R=base_r / MOON_RADIUS, treatments={})
    for name, t in temps.items():
        share = broken_share(t, base_r)
        rows = {}
        for ground, x in water.items():
            for label, w in (('', x), ('_cold_cloud', x * COLD_CLOUD_WATER)):
                kg = escape_rates(t, w, share['share'])
                rows[f'{ground}{label}'] = dict(water_share=w, **kg, cycles=cycle_times(kg['total']))
        out['treatments'][name] = dict(temperature_k=t, methane=share, grounds=rows)
    lte = out['treatments']['lte']['grounds']['ts288']
    out['check'] = dict(oxygen_step_kg_s=list(oxygen), limiting_kg_s=lte['without_methane'],
                        ratio=[x / lte['without_methane'] for x in oxygen])
    return out


def main(argv=None) -> int:
    out = analyse()
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    out['evidence'] = ('The diffusion-limited escape of hydrogen in each form through the homopause (Hunten 1973), with '
                       'measured binary diffusion and a steady tracer for the glow\'s breaking of methane above the '
                       'base; screening estimates from the middle atmosphere\'s design cases and the oxygen step, '
                       'without methane\'s chemistry below the base or what its hydrogen does to the air\'s H2.')
    out['reading_rule'] = ('treatments: for each design treatment\'s homopause temperature, methane\'s share broken '
                           'below its homopause and, for each ground\'s water past the tropopause (and the cold '
                           'cloud\'s 2.5 times it), the hydrogen escape in kg/s by carrier: total takes methane at '
                           'that share, total_unbroken and total_broken its bounds; cycles gives years for the total '
                           'to take the air\'s hydrogen in H2 and methane, methane\'s hydrogen and the Moon\'s waters, '
                           'and the water a year against the delivery. check sets the limiting flux of H2 and water '
                           'beside the oxygen step\'s transport.')
    out['producer'] = dict(domain='atmosphere', files={'loss_response/hydrogen.py': digest(__file__)},
                           inputs={str(p.relative_to(ROOT)): dict(path=str(p.relative_to(ROOT)), sha256=digest(p))
                                   for p in [CURVED, OXYGEN] + [CASES / f'{c}.json' for c in GROUNDS.values()]})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f"limiting flux of H2 and water {out['check']['limiting_kg_s']:.5f} kg/s against the oxygen step's "
          f"{out['check']['oxygen_step_kg_s']}")
    for name, tr in out['treatments'].items():
        m = tr['methane']
        print(f"{name} ({tr['temperature_k']:.1f} K): methane broken {m['share']:.2f} below its homopause "
              f"({m['height_to_homopause_km']:.0f} km up, decay {m['decay_length_km']:.0f} km, glow lifetime "
              f"{m['glow_lifetime_days']:.0f} days)")
        for g, r in tr['grounds'].items():
            c = r['cycles']
            print(f"   {g}: water {r['water_share']:.2e}; H2 {r['h2']:.5f}, water {r['water']:.6f}, methane "
                  f"{r['methane']:.5f} ({r['methane_unbroken']:.5f}-{r['methane_broken']:.5f}); total {r['total']:.4f} "
                  f"kg/s; air's hydrogen {c['air_hydrogen_years']:.2e} yr, waters {c['waters_years']:.2e} yr")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
