"""The array's computing: what the ring fleet's light can supply against what computing tasks demand.

    python -m provisioning.computing
    # -> provisioning/results/computing.json

Supply. Computing draws electricity from collectors flying with the fleet on platforms of their own, in light the
Moon never needs: the night half intercepts about 850 PW (joint synthesis). Per terawatt of computing the shield
study gives collectors, radiators at 330 K and their masses (array_heat.json); this screen adds collector and processor
mass from the literature, the energy per operation from today's accelerators to the irreversible floor, and where
the heat lands. Scenarios run from 1 to 10,000 TW.

Demand. Each task's computing comes from a published estimate or from arithmetic on a published measurement, with
the assumptions stated in its basis: the fleet's control, the Moon's weather and climate twin from kilometre to metre
scale, cell-level biology, radio correlation, cosmology, lattice QCD, molecular simulation, AI training and services
for 20-30 billion people, digital minds, ancestor simulations and the limits of computation.

Comparison. For each task, its power at today's efficiency and at the bounds, the share of the night half's light
it takes, the heat it puts on the Moon from orbit and if it ran on the Moon, and where it can run by latency.

Every literature value cites its source (Author year), each listed in computing_sources.json.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path

from shared.constants import (AU, BOLTZMANN, EARTH_GM, EARTH_MOON_DISTANCE, EARTH_RADIUS, JULIAN_DAY, JULIAN_YEAR_DAYS,
                              MOON_GM, MOON_RADIUS, PLANCK, SOLAR_CONSTANT, SPEED_OF_LIGHT, STEFAN_BOLTZMANN, SUN_GM)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.provisioning.computing/1'
HEAT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'array_heat.json'
LEDGERS = ROOT / 'research' / 'studies' / 'joint_synthesis' / 'results' / 'ledgers.json'
NIGHT = ROOT / 'illumination' / 'fleet_light' / 'results' / 'night_light.json'
SCREEN = HERE / 'results' / 'first_screen.json'
SOURCES = HERE / 'computing_sources.json'
OUT = HERE / 'results' / 'computing.json'
INPUTS = (HEAT, LEDGERS, NIGHT, SCREEN, SOURCES)

YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY
DAY = JULIAN_DAY
LN2 = math.log(2.0)
HBAR = PLANCK / (2 * math.pi)
SUN_LUMINOSITY = 4 * math.pi * AU ** 2 * SOLAR_CONSTANT          # W, from the solar constant at 1 AU
MOON_AREA = 4 * math.pi * MOON_RADIUS ** 2
EARTH_AREA = 4 * math.pi * EARTH_RADIUS ** 2
LUMENS_PER_PW = 0.94e17                                          # sunlight, the fleet-light product's reading rule

# Supply scenarios, TW of electricity used for computing.
SCENARIOS_TW = (1.0, 10.0, 100.0, 1000.0, 10000.0)

# Energy per operation.
H100_FP16_DENSE = 989.4e12        # FLOP/s, dense 16-bit, 1,979 with sparsity (NVIDIA, H100 specifications)
H100_FP32 = 67e12                 # FLOP/s (NVIDIA, H100 specifications)
H100_FP64 = (34e12, 67e12)        # FLOP/s, vector and tensor core (NVIDIA, H100 specifications)
H100_W = 700.0                    # W, maximum thermal design power (NVIDIA, H100 specifications)
P100_FP32, P100_W = 9.3e12, 250.0     # Tesla P100 PCIe (NVIDIA 2016)
A100_FP32, A100_W = 19.5e12, 400.0    # A100 SXM4 (NVIDIA 2020)
CMOS_LIMIT_FP4_PER_J = 4.7e15     # geometric mean of the maximum CMOS efficiency, 4-bit operations (Ho et al. 2023)
CMOS_LIMIT_SPREAD_OOM = 0.7       # its standard deviation in orders of magnitude (Ho et al. 2023)
FP16_OVER_FP4 = 16.0              # ideal energy grows with the square of precision, (16/4)^2 (Ho et al. 2023)
SWITCHES_PER_FP4 = (30.0, 3000.0) # transistor switches an optimal device needs per 4-bit operation (Ho et al. 2023)
EFFICIENCY_DOUBLING_YEARS = (1.57, 2.7)   # 1946-2000 and after 2000, as quoted by Ho et al. 2023
RADIATOR_TEMPERATURES_K = (150.0, 250.0, 330.0, 500.0)

# Mass beyond the radiators that array_heat.json gives.
COLLECTOR_W_PER_KG = (1000.0, 100.0)  # electric W per kg of array system: the shield study's exploratory 1,000 W/kg and
                                      # its low 100 W/kg case (solar_shield_array/electromagnetic.md; 300 W/kg central)
STARCLOUD_KG_PER_KW = 100e3 / 40e3    # about 40 MW of computing in a launch container of about 100 t (Feilden et al. 2024)
STARCLOUD_SHIELD_KG_PER_KW = 1.0      # radiation shielding per kW of computing (Feilden et al. 2024)
DGX_H100_KG, DGX_H100_KW = 130.45, 10.2   # server weight and maximum system power (NVIDIA 2024, DGX H100)
HARDWARE_LIFE_YEARS = 5.0             # useful mean life of a GPU (Ho et al. 2023)

# Today's world, for scale.
DATA_CENTRES_TWH_2024 = 415.0         # the world's data centres in 2024 (IEA 2025)
AI_H100_EQUIVALENTS_2026 = 15e6       # installed AI accelerators, January 2026, drawing over 10 GW (Epoch AI 2026)
NVIDIA_STOCK_FLOP_S_2024 = 4e21       # dense 16-bit, Q4 2024, about 4 million H100-equivalents (Epoch AI 2025a)
TOP500_LEAD_FLOP_S = 2.198e18         # LineShine, HPL at 64 bits, June 2026 (TOP500 2026)
EL_CAPITAN = (1.809e18, 29.685e6)     # HPL FLOP/s and W, June 2026 (TOP500 2026)
GOOGLE_TOKENS_PER_MONTH = 1.3e15      # all Google surfaces, third quarter of 2025 (Alphabet 2025)
WORLD_POPULATION_2024 = 8.2e9         # (UN 2024)
METROPOLIS_PEOPLE = 100e6             # the summit metropolis (decisions register, 27 September 2026)

# The array's control.
HPSC_FLOP_S = 1e12                    # a radiation-hardened spaceflight processor, 1 TFLOP/s in BF16 (Microchip 2024)
TILE_FILTER_STATES = 20               # orbit, attitude, store and trim states per tile (assumed for scale)
FILTER_HZ = 1.0                       # updates per second (assumed)
FILTER_ENSEMBLE = 1000                # members of a fleet-wide ensemble (assumed)

# The Moon's weather and climate twin, scaled from a measured km-scale run.
COSMO_DX_M = 930.0                    # near-global COSMO grid spacing (Schulthess et al. 2019; Fuhrer et al. 2018)
COSMO_LEVELS = 60                     # surface to 25 km (Schulthess et al. 2019)
COSMO_DT_S = 6.0                      # split-explicit step (Schulthess et al. 2019)
COSMO_EARTH_SHARE = 0.984             # share of Earth's surface covered (Schulthess et al. 2019)
COSMO_MWH_PER_SY = 596.0              # energy per simulated year on Piz Daint's P100 nodes (Schulthess et al. 2019)
COUPLING_FACTOR = 1.2                 # land surface, ocean, waves and sea ice coupled, Table 2 (Schulthess et al. 2019)
GOAL_LEVELS = (180, 360)              # the 1 km goal's 180 levels to 100 km (Schulthess et al. 2019), doubled at the
                                      # high end for the Moon's taller air (assumed)
REFINED_LAYER_M = 1000.0              # the lowest kilometre refined to the grid spacing (assumed)
GRIDS_M = (1000.0, 100.0, 10.0, 1.0)
FORECAST_SERVICE = 50 * 20            # an ensemble of 50 members run 20 times faster than real time (assumed)

# Cells and soil.
MINIMAL_CELL_GPU_HOURS = 15000.0 / 50 # A100 GPU-hours per 105-minute cell cycle, 50 replicates (Thornburg et al. 2026)
MINIMAL_CELL_CYCLE_S = 105 * 60.0     # JCVI-syn3A doubling time (Thornburg et al. 2026)
SOIL_CELLS_PER_G = 1e10               # up to ten billion microbial cells per gram of soil (Torsvik & Ovreas 2002)
EARTH_PROKARYOTES = (4e30, 6e30)      # (Whitman et al. 1998)
MINIMAL_CELL_ATOMS = 6e9              # more than six billion atoms, 561 million coarse-grained beads (Stevens et al. 2023)
MD_FLOP_PER_ENTITY_STEP = 1e3         # the roadmap's level-10 cost per molecule and step (Sandberg & Bostrom 2008)
MD_STEPS_PER_S = 1e15                 # femtosecond steps (Sandberg & Bostrom 2008)

# Science.
FARVIEW_DIPOLES = 1e5                 # (Polidan et al. 2024; NASA, FarView)
FARVIEW_BAND_HZ = 40e6 - 5e6          # 5-40 MHz (NASA, FarView)
LARGE_ARRAY_ELEMENTS = 1e7            # a hundred times FarView's dipoles (for scale)
FLOP_PER_CMAC = 8                     # a complex multiply-accumulate in real operations
FFT_FLOP_PER_POINT = 5.0              # a complex FFT takes about 5 N log2 N real operations; N log N (Tegmark & Zaldarriaga 2009)
SKA_SDP_FLOP_S = 2 * 135e15           # two science data processors of about 135 PFLOP/s (SKAO, Science Processing Centres)
FRONTIER_E = (513.1e15, 7 * DAY, 4e12)    # FLOP/s peak, just over one week, particles (Frontiere et al. 2025)
LQCD_EXAFLOP_HOURS = (1.5, 12000.0)   # proposed ensembles, sustained exaflop-hours (Boyle et al. 2022)
FOLDING_AT_HOME_FLOP_S = 1.01e18      # peak, 2020 (Zimmerman et al. 2021)
ALPHAFOLD_CHIPS = 128 / 2             # 128 TPU v3 cores (Jumper et al. 2021), two per chip (Google Cloud, TPU v3)
ALPHAFOLD_DAYS = 7 + 4                # about a week, then about four more days (Jumper et al. 2021)
TPU_V3_FLOP_S = 123e12                # per chip, bf16 (Google Cloud, TPU v3)

# AI.
FLOP_PER_TOKEN_PER_PARAMETER = 2      # forward pass (Kaplan et al. 2020)
TOKENS_PER_S_PER_PERSON = (100.0, 1e4)    # served to each person (assumed)
MODEL_PARAMETERS = 1e12               # active parameters (assumed)
POPULATION = (20e9, 30e9)             # the author's working scale (decisions register, 9 October 2026)
BRAIN_EQUIVALENTS_PER_PERSON = (1.0, 100.0)   # AI services per person, in human brains (assumed)
LARGEST_RUN_FLOP_2025 = 4.6e26        # Grok 3 (Epoch AI 2025b)
GROK4_GWH = 310.0                     # Grok 4's training energy (Epoch AI 2025c)
FRONTIER_GROWTH_PER_YEAR = (4.0, 5.0) # training compute of frontier models (Sevilla & Roldan 2024)
RUN_2030_FLOP = 2e29                  # feasible by 2030 (Sevilla et al. 2024)
LATENCY_WALL_FLOP = (3e30, 1e32)      # (Sevilla et al. 2024)
TRAINING_RUN_DAYS = 100.0             # (assumed)

# Minds and limits.
BRAIN_FLOP_S = (1e13, 1e15, 1e17, 1e21)   # mechanistic range, the central figure, the <10% bound (Carlsmith 2020)
WBE_LEVELS = {'4 spiking neural network': 1e18, '5 electrophysiology': 1e22, '6 metabolome': 1e25,
              '7 proteome': 1e26, '8 states of protein complexes': 1e27, '9 distribution of complexes': 1e30,
              '10 stochastic behaviour of single molecules': 1e43}     # Table 9 (Sandberg & Bostrom 2008)
DIGITAL_MIND_FLOP_S = (1e15, 1e18)    # Carlsmith's central figure to the roadmap's level 4
EMS = 1e12                            # "trillions" of ems (Hanson 2016)
EM_SPEED = 1000.0                     # a typical em's speed against a human's (Hanson 2016)
ANCESTOR_OPS = (1e33, 1e36)           # (Bostrom 2003)
ANCESTOR_TERMS = (100e9, 50.0, 30e6, (1e14, 1e17))   # people, years each, s a year, operations/s a brain (Bostrom 2003)
EVOLUTION_ANCHOR_FLOP = 1e41          # (Cotra 2020, as summarised by Ho 2022)
MATRIOSHKA = (3e42, 3e26)             # multi-layer: operations/s and W (Bradbury 1999)
PLANETARY_COMPUTER_OPS = 1e42         # (Bostrom 2003, after Bradbury 1999)
LLOYD_LAPTOP = (5.4258e50, 2.13e31)   # operations/s and bits for 1 kg in 1 litre (Lloyd 2000)
ORDINARY_MATTER_OPS_PER_KG = 1e40     # (Lloyd 2000)
UNIVERSE = (1e120, 1e90)              # operations and bits (Lloyd 2002)
AESTIVATION_GAIN = 1e30               # (Sandberg, Armstrong & Cirkovic 2017)
SCHNEIER_T_K = 3.2                    # (Schneier 1996)
CMB_K = 2.7

# Latency and links.
CONVERSATION_GAP_S = 0.208            # mean response offset across ten languages (Stivers et al. 2009)
TBIRD_BPS = 200e9                     # from low Earth orbit (Schieler et al. 2023)
LLCD_BPS = 622e6                      # from lunar orbit (Boroson & Robinson 2014)
VIDEO_BPS = 10e6                      # one stream per person (assumed)
BYTES_PER_TOKEN = 4.0                 # (assumed)
EARTH_POPULATION_PEAK = 10.3e9        # (UN 2024)

EVIDENCE = ('A screen on committed products: the shield study\'s array heat (collectors, radiators, placement, fleet, '
            'film plant), the joint synthesis\'s night-half interception, the fleet-light model\'s night-side '
            'reflection and the provisioning first screen\'s heat response, with the literature values in '
            'computing_sources.json. Closed-form arithmetic; no platform, device, network or workload is modelled. '
            'Task demands are published estimates or arithmetic on published measurements, each with the assumptions '
            'named in its basis; the per-person AI demand, the weather twin\'s levels and service factor, and the '
            'control filter sizes are assumptions for scale.')
READING_RULE = ('Power in W unless a key names TW or PW; computing in FLOP/s at 16 bits unless stated, or total FLOP '
                'for campaigns. "today" is the H100\'s 1.4e12 FLOP/J at 16 bits; tasks measured on older GPUs carry '
                'their measured energy scaled to the H100 by 32-bit FLOP per joule, and their FLOP/s are '
                'H100-equivalent. "cmos_limit" divides today\'s energy by the ratio of Ho et al.\'s 16-bit limit to '
                'the H100; "landauer_330K" by the ratio of the irreversible floor at 330 K. Light shares are of the '
                'night half\'s interception at 20-30% conversion. Heat on the Moon in W/m2 over its whole surface, '
                'warming in K by the first screen\'s 0.63-1.0 K per W/m2. Ranges pair the low ends and the high ends.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def sig(x, n=3):
    """x rounded to n significant figures, recursively through lists, tuples and dicts."""
    if isinstance(x, dict):
        return {k: sig(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [sig(v, n) for v in x]
    if isinstance(x, bool) or x is None or isinstance(x, str):
        return x
    x = float(x)
    if x == 0 or not math.isfinite(x):
        return x
    return float(f'{x:.{n}g}')


def moon_share(distance_m):
    """Share of heat radiated evenly from a point at distance_m that strikes the Moon: its solid angle over 4 pi."""
    s = MOON_RADIUS / distance_m
    return (1 - math.sqrt(1 - s * s)) / 2


def efficiency(today):
    """FLOP per joule at 16 bits: today's accelerator, the CMOS limit and the irreversible floor."""
    cmos = CMOS_LIMIT_FP4_PER_J / FP16_OVER_FP4
    switches = [s * FP16_OVER_FP4 for s in SWITCHES_PER_FP4]
    def floor(t):
        bit = BOLTZMANN * t * LN2
        return [1 / (switches[1] * bit), 1 / (switches[0] * bit)]
    gen = dict(p100=P100_FP32 / P100_W, a100=A100_FP32 / A100_W, h100=H100_FP32 / H100_W)
    return dict(
        today_flop_per_j=today, h100_check=H100_FP16_DENSE / H100_W,
        cmos_limit_flop_per_j=cmos,
        cmos_limit_range_flop_per_j=[cmos / 10 ** CMOS_LIMIT_SPREAD_OOM, cmos * 10 ** CMOS_LIMIT_SPREAD_OOM],
        cmos_over_today=cmos / today,
        years_to_cmos_at_past_doubling=[math.log2(cmos / today) * d for d in EFFICIENCY_DOUBLING_YEARS],
        switches_per_fp16=switches,
        landauer_bit_j={f'{t:g}K': BOLTZMANN * t * LN2 for t in (300.0, 330.0, CMB_K)},
        landauer_floor_flop_per_j={f'{t:g}K': floor(t) for t in RADIATOR_TEMPERATURES_K},
        floor_330K_over_today=[f / today for f in floor(330.0)],
        today_j_per_flop_in_kT_ln2_at_330K=1 / today / (BOLTZMANN * 330.0 * LN2),
        fp32_flop_per_j=gen, h100_over_p100_fp32=gen['h100'] / gen['p100'], h100_over_a100_fp32=gen['h100'] / gen['a100'],
        fp64_energy_over_fp16=[H100_FP16_DENSE / H100_FP64[1], H100_FP16_DENSE / H100_FP64[0]],
        # The CMOS limit scales as the square of precision; the H100's own units do not.
        cmos_gain_by_precision=dict(fp16=cmos / today, fp32=cmos / 4 / gen['h100'],
                                    fp64=[cmos / 16 / (H100_FP64[1] / H100_W), cmos / 16 / (H100_FP64[0] / H100_W)]),
        reversible='no floor: logically reversible, adiabatic logic can dissipate below kT ln 2 per operation, trading '
                   'speed for energy (Bennett 1973; Frank et al. 2020)')


def radiators(design):
    """Radiator area, mass and the Landauer-limited rate per square metre at each temperature."""
    eps, kg = design['radiator_emissivity'], design['radiator_kg_m2']
    out = {}
    for t in RADIATOR_TEMPERATURES_K:
        flux = 2 * eps * STEFAN_BOLTZMANN * t ** 4                 # two faces to deep space
        area = 1e12 / flux                                        # m2 per TW
        bit = BOLTZMANN * t * LN2
        out[f'{t:g}K'] = dict(net_w_m2=flux, km2_per_tw=area / 1e6, mt_per_tw=[area * k / 1e9 for k in kg],
                              landauer_bits_per_j=1 / bit, landauer_bits_per_s_per_m2=flux / bit)
    ref = out['330K']['landauer_bits_per_s_per_m2']
    for v in out.values():
        v['per_m2_over_330K'] = v['landauer_bits_per_s_per_m2'] / ref
    return out


def supply(ah, night_w, night_area_km2, turned_w, eff, ctx):
    comp, fleet, design, placement = ah['computing'], ah['fleet'], ah['design'], ah['placement']
    collectors = sorted({c['conversion']: c['heat_per_electric'] for c in ah['collectors']}.items())
    eta = [collectors[1][0], collectors[0][0]]                    # 0.3 then 0.2: low ends pair with low ends
    heat_per_electric = [collectors[1][1], collectors[0][1]]      # 1.93 then 3.4 W of collector heat per W
    collector_km2 = sorted(comp['collector_km2_per_TW'])
    radiator_km2 = comp['radiator_panel_km2_per_TW']
    radiator_mt = comp['radiator_Mt_per_TW']
    collector_mt = [1e12 / w / 1e9 for w in COLLECTOR_W_PER_KG]
    hardware_kg_per_kw = [STARCLOUD_KG_PER_KW + STARCLOUD_SHIELD_KG_PER_KW, DGX_H100_KG / DGX_H100_KW]
    hardware_mt = list(hardware_kg_per_kw)                        # a kilogram per kW is a megatonne per TW
    mass_mt = [r + c + h for r, c, h in zip(radiator_mt, collector_mt, hardware_mt)]
    fleet_gt = sum(fleet['optical_mass_Gt']) / 2
    fleet_km2 = sum(fleet['area_km2']) / 2
    upkeep = fleet['film_upkeep_Gt_per_year']
    orbit_w_m2 = placement['moon_W_m2_per_TW_released_in_orbit']
    moon_w_m2 = placement['moon_W_m2_per_TW_used_on_moon']
    k_per_w_m2 = ctx['k_per_w_m2']
    released = [1 + h for h in heat_per_electric]                 # TW of heat released in orbit per TW of computing
    orbit_k_per_tw = [released[0] * orbit_w_m2 * k_per_w_m2[0], released[1] * orbit_w_m2 * k_per_w_m2[1]]
    intercepted = ah['tiles']['intercepted_PW']
    per_tw = dict(
        light_tw=[1 / e for e in eta], collector_km2=collector_km2, radiator_km2=radiator_km2, radiator_mt=radiator_mt,
        collector_mt=collector_mt, hardware_kg_per_kw=hardware_kg_per_kw, hardware_mt=hardware_mt, mass_mt=mass_mt,
        hardware_renewal_mt_per_year=[h / HARDWARE_LIFE_YEARS for h in hardware_mt],
        heat_released_in_orbit_tw=released, orbit_moon_w_m2=[r * orbit_w_m2 for r in released],
        orbit_moon_k=orbit_k_per_tw, on_moon_w_m2=moon_w_m2, on_moon_k=[moon_w_m2 * k for k in k_per_w_m2],
        platform_side_km=[math.sqrt(c + radiator_km2) for c in collector_km2],
        flop_s_per_kg_today=[1e12 * eff['today_flop_per_j'] / (m * 1e9) for m in mass_mt[::-1]])
    scen = {}
    for tw in SCENARIOS_TW:
        light_pw = [tw / e / 1e3 for e in eta]
        area = [tw * (c + radiator_km2) for c in collector_km2]
        mass_gt = [tw * m / 1e3 for m in mass_mt]
        scen[f'{tw:g}'] = dict(
            flop_s=dict(today=tw * 1e12 * eff['today_flop_per_j'], cmos_limit=tw * 1e12 * eff['cmos_limit_flop_per_j'],
                        landauer_330K=[tw * 1e12 * f for f in eff['landauer_floor_flop_per_j']['330K']]),
            light_pw=light_pw, share_of_night_half_light=[p * 1e15 / night_w for p in light_pw],
            share_of_fleet_interception=[light_pw[0] / intercepted[1], light_pw[1] / intercepted[0]],
            collector_and_radiator_km2=area, share_of_fleet_area=[a / fleet_km2 for a in area],
            share_of_night_half_tile_area=[a / night_area_km2 for a in area],
            radiator_mt=[tw * r for r in radiator_mt], collector_mt=[tw * c for c in collector_mt],
            hardware_mt=[tw * h for h in hardware_mt], mass_gt=mass_gt, share_of_fleet_mass=[m / fleet_gt for m in mass_gt],
            years_of_film_plant_throughput=[mass_gt[0] / upkeep[1], mass_gt[1] / upkeep[0]],
            hardware_renewal_gt_per_year=[tw * h / HARDWARE_LIFE_YEARS / 1e3 for h in hardware_mt],
            orbit_moon_w_m2=[tw * x for x in per_tw['orbit_moon_w_m2']], orbit_moon_k=[tw * x for x in orbit_k_per_tw],
            orbit_share_of_glow_low_end=[tw * x / ctx['glow_w_m2'][0] for x in per_tw['orbit_moon_w_m2']],
            if_on_moon_k=[tw * x for x in per_tw['on_moon_k']],
            times_data_centres_2024=tw / ctx['data_centres_tw'],
            times_ai_compute_2026=tw * 1e12 * eff['today_flop_per_j'] / ctx['ai_flop_s_2026'],
            brain_equivalents_today=tw * 1e12 * eff['today_flop_per_j'] / BRAIN_FLOP_S[1],
            brain_equivalents_cmos_limit=tw * 1e12 * eff['cmos_limit_flop_per_j'] / BRAIN_FLOP_S[1])
    # Heat from orbital computing falls with the inverse square of its distance from the Moon.
    hill = EARTH_MOON_DISTANCE * (MOON_GM / (3 * EARTH_GM)) ** (1 / 3)
    sel1 = AU * (EARTH_GM / (3 * SUN_GM)) ** (1 / 3)
    distances = {'ring_radius': design['orbit_radius_km'] * 1e3, 'earth_moon_l1_l2': hill, 'at_100000_km': 1e8,
                 'sun_earth_l1_l2': sel1}
    by_distance = {k: dict(distance_km=d / 1e3, moon_w_m2_per_tw_heat=1e12 * moon_share(d) / MOON_AREA)
                   for k, d in distances.items()}
    # A flat emitter with one of its two faces turned to the Moon sends it half its heat times the view factor (R/d)^2.
    ring = distances['ring_radius']
    face_on = dict(share=(MOON_RADIUS / ring) ** 2 / 2, over_isotropic=(MOON_RADIUS / ring) ** 2 / 2 / moon_share(ring))
    k_bounds = {f'{k:g}K': [k / x for x in orbit_k_per_tw[::-1]] for k in (0.1, 1.0)}
    return dict(
        night_half_intercepts_pw=night_w / 1e15, night_half_tile_km2=night_area_km2, night_half_turns_pw=turned_w / 1e15,
        turned_light_supports_tw=[turned_w / 1e12 * e for e in eta[::-1]],
        fleet_intercepts_pw=intercepted, fleet_area_km2=fleet_km2, fleet_mass_gt=fleet_gt, film_plant_gt_per_year=upkeep,
        film_plant_tw=ah['film_plant']['power_TW'], conversion=eta, per_tw=per_tw, scenarios=scen,
        orbital_tw_for_warming=k_bounds, heat_by_distance=by_distance, face_turned_to_moon=face_on,
        fleet_own_loads_tw=dict(trim_bistable=ah['trim']['bistable']['fleet_TW'],
                                trim_electrochromic=ah['trim']['electrochromic']['fleet_TW'],
                                store=[min(r['fleet_GW'][0] for r in ah['store']['rotors']) / 1e3,
                                       max(r['fleet_GW'][1] for r in ah['store']['rotors']) / 1e3]))


def weather(eff):
    """The Moon's weather and climate twin, from COSMO's measured energy per cell and step."""
    cosmo_cells = EARTH_AREA * COSMO_EARTH_SHARE / COSMO_DX_M ** 2 * COSMO_LEVELS
    j_per_cell_step = COSMO_MWH_PER_SY * 3.6e9 / (cosmo_cells * YEAR / COSMO_DT_S)
    j_today = j_per_cell_step / eff['h100_over_p100_fp32']
    rows = {}
    for dx in GRIDS_M:
        refined = REFINED_LAYER_M / dx if dx < 1000 else 0.0
        levels = [g + refined for g in GOAL_LEVELS]
        dt = COSMO_DT_S * dx / COSMO_DX_M                          # the step shrinks with the grid
        cells = [MOON_AREA / dx ** 2 * n for n in levels]
        real_time = [c / dt * j_today * COUPLING_FACTOR for c in cells]
        rows[f'{dx:g}m'] = dict(levels=levels, step_s=dt, cells=cells, real_time_w=real_time,
                                forecast_service_w=[w * FORECAST_SERVICE for w in real_time])
    return dict(cosmo_cells=cosmo_cells, j_per_cell_step_p100=j_per_cell_step, j_per_cell_step_today=j_today,
                coupling=COUPLING_FACTOR, forecast_service_factor=FORECAST_SERVICE, grids=rows)


def demand(ctx, eff, wx):
    """Each task's computing: FLOP/s for continuous services, FLOP for campaigns, with its power today."""
    today = eff['today_flop_per_j']
    rate = lambda f: dict(kind='rate', flop_s=list(f), power_today_w=[x / today for x in f])
    measured = lambda w: dict(kind='rate', flop_s=[x * today for x in w], power_today_w=list(w), h100_equivalent=True)
    campaign = lambda f: dict(kind='campaign', flop=list(f), energy_today_j=[x / today for x in f])
    tiles = ctx['tiles']
    filter_flop_s = [n * 2 * TILE_FILTER_STATES ** 3 * FILTER_HZ * FILTER_ENSEMBLE for n in tiles]
    cell_w = MINIMAL_CELL_GPU_HOURS * 3600 * A100_W / MINIMAL_CELL_CYCLE_S / eff['h100_over_a100_fp32']
    fx = lambda n: n * (n - 1) / 2 * FLOP_PER_CMAC * FARVIEW_BAND_HZ
    fft = lambda n: FFT_FLOP_PER_POINT * n * math.log2(n) * FARVIEW_BAND_HZ
    g = wx['grids']
    people_brains = [POPULATION[0] * BRAIN_EQUIVALENTS_PER_PERSON[0] * BRAIN_FLOP_S[1],
                     POPULATION[1] * BRAIN_EQUIVALENTS_PER_PERSON[1] * BRAIN_FLOP_S[1]]
    tokens = [POPULATION[0] * TOKENS_PER_S_PER_PERSON[0] * FLOP_PER_TOKEN_PER_PARAMETER * MODEL_PARAMETERS,
              POPULATION[1] * TOKENS_PER_S_PER_PERSON[1] * FLOP_PER_TOKEN_PER_PARAMETER * MODEL_PARAMETERS]
    training = [f / (TRAINING_RUN_DAYS * DAY) for f in LATENCY_WALL_FLOP]
    product = ANCESTOR_TERMS[0] * ANCESTOR_TERMS[1] * ANCESTOR_TERMS[2]
    t = {}
    # The array itself.
    t['fleet_control_tile_processors'] = dict(rate([n * HPSC_FLOP_S for n in tiles]), where='tiles', latency_s=1e-3,
        basis='one 1 TFLOP/s spaceflight processor on each of 28.1-28.2 million tiles (Microchip 2024)')
    t['fleet_twin_ensemble'] = dict(rate(filter_flop_s), where='array', latency_s=1.0,
        basis='every tile\'s 20-state filter at 1 Hz in a 1,000-member ensemble (assumed sizes, 2n^3 FLOP per update)')
    # The Moon's management.
    for dx, label in (('1000m', '1 km'), ('100m', '100 m'), ('10m', '10 m'), ('1m', '1 m')):
        t[f'weather_{dx}_forecast'] = dict(measured(g[dx]['forecast_service_w']), where='array', latency_s=60.0,
            basis=f'the Moon\'s coupled weather and climate twin at {label}, 50 members at 20 times real time, scaled '
                  f'from COSMO\'s measured 596 MWh per simulated year (Schulthess et al. 2019)')
    t['weather_1m_real_time'] = dict(measured(g['1m']['real_time_w']), where='array', latency_s=1.0,
        basis='the same at 1 m, one member in real time')
    t['cells_one_gram_of_soil'] = dict(measured([SOIL_CELLS_PER_G * cell_w] * 2), where='anywhere', latency_s=None,
        basis='a cell-by-cell twin of 1 g of soil, 1e10 cells (Torsvik & Ovreas 2002), each a minimal cell in 4D at '
              'about 300 A100 GPU-hours per 105-minute cycle (Thornburg et al. 2026)')
    t['cells_earth_prokaryotes'] = dict(measured([n * cell_w for n in EARTH_PROKARYOTES]), where='beyond', latency_s=None,
        basis='the same for Earth\'s 4-6e30 prokaryotes (Whitman et al. 1998)')
    # Science.
    t['radio_farview_fx'] = dict(rate([fx(FARVIEW_DIPOLES)] * 2), where='near the telescope', latency_s=None,
        basis='full correlation of 100,000 dipoles of one polarization over 5-40 MHz, 8 FLOP per baseline and sample '
              '(Polidan et al. 2024; NASA, FarView); two polarizations take four times as much')
    t['radio_1e7_fx'] = dict(rate([fx(LARGE_ARRAY_ELEMENTS)] * 2), where='near the telescope', latency_s=None,
        basis='the same for ten million elements')
    t['radio_1e7_fft'] = dict(rate([fft(LARGE_ARRAY_ELEMENTS)] * 2), where='near the telescope', latency_s=None,
        basis='ten million elements on a regular grid, correlated by FFT at 5 N log2 N per sample (Tegmark & '
              'Zaldarriaga 2009)')
    t['ska_science_data'] = dict(rate([SKA_SDP_FLOP_S] * 2), where='anywhere', latency_s=None,
        basis='the SKA\'s two science data processors at about 135 PFLOP/s each (SKAO, Science Processing Centres)')
    t['folding_at_home_2020'] = dict(rate([FOLDING_AT_HOME_FLOP_S] * 2), where='anywhere', latency_s=None,
        basis='Folding@home at its 2020 peak, which simulated 0.1 s of the SARS-CoV-2 proteome (Zimmerman et al. 2021)')
    t['alphafold2_training'] = dict(campaign([ALPHAFOLD_CHIPS * TPU_V3_FLOP_S * ALPHAFOLD_DAYS * DAY] * 2),
        where='anywhere', latency_s=None, basis='128 TPU v3 cores for about 11 days at peak (Jumper et al. 2021; '
                                                'Google Cloud, TPU v3)')
    t['cosmology_frontier_e'] = dict(campaign([FRONTIER_E[0] * FRONTIER_E[1]] * 2), where='anywhere', latency_s=None,
        basis='Frontier-E, 4 trillion particles at 513.1 PFLOP/s peak for about a week (Frontiere et al. 2025)')
    t['lattice_qcd_ensembles'] = dict(campaign([x * 1e18 * 3600 for x in LQCD_EXAFLOP_HOURS]), where='anywhere',
        latency_s=None, basis='the proposed ensembles, 1.5-12,000 sustained exaflop-hours (Boyle et al. 2022)')
    t['all_atom_minimal_cell_cycle'] = dict(
        campaign([MINIMAL_CELL_ATOMS * MD_FLOP_PER_ENTITY_STEP * MD_STEPS_PER_S * MINIMAL_CELL_CYCLE_S] * 2),
        where='anywhere', latency_s=None, basis='6e9 atoms (Stevens et al. 2023) for one 105-minute cycle at 1e3 FLOP '
        'per atom and femtosecond step, the roadmap\'s cost per molecule (Sandberg & Bostrom 2008)')
    # AI.
    t['ai_world_2026'] = dict(rate([ctx['ai_flop_s_2026']] * 2), where='Earth', latency_s=None,
        basis='15 million H100-equivalents at 989 TFLOP/s, drawing over 10 GW (Epoch AI 2026)')
    t['ai_largest_run_2025'] = dict(campaign([LARGEST_RUN_FLOP_2025] * 2), where='anywhere', latency_s=None,
        basis='Grok 3 at 4.6e26 FLOP (Epoch AI 2025b); Grok 4 took about 310 GWh (Epoch AI 2025c)')
    t['ai_training_latency_wall'] = dict(rate(training), where='anywhere, compact', latency_s=None,
        basis='one run of 3e30-1e32 FLOP, the latency wall (Sevilla et al. 2024), over 100 days (assumed)')
    t['ai_services_20_30_billion'] = dict(rate(people_brains), where='array for the Moon', latency_s=CONVERSATION_GAP_S,
        basis='1-100 brain-equivalents of 1e15 FLOP/s (Carlsmith 2020) for each of 20-30 billion people (assumed use)')
    t['ai_services_token_check'] = dict(rate(tokens), where='array for the Moon', latency_s=CONVERSATION_GAP_S,
        basis='100-10,000 tokens/s per person from 1e12-parameter models at 2N FLOP per token (Kaplan et al. 2020)')
    # Minds and the far out.
    t['digital_population_like_ours'] = dict(
        rate([POPULATION[0] * DIGITAL_MIND_FLOP_S[0], POPULATION[1] * DIGITAL_MIND_FLOP_S[1]]),
        where='array or anywhere', latency_s=CONVERSATION_GAP_S,
        basis='20-30 billion digital minds at 1e15 (Carlsmith 2020) to 1e18 FLOP/s (level 4, Sandberg & Bostrom 2008)')
    t['em_economy'] = dict(rate([EMS * EM_SPEED * x for x in DIGITAL_MIND_FLOP_S]), where='anywhere', latency_s=None,
        basis='a trillion ems at 1,000 times human speed (Hanson 2016), 1e15-1e18 FLOP/s each at human speed')
    for name, f in WBE_LEVELS.items():
        t['wbe_level_' + name.split()[0]] = dict(rate([f, f]), where='anywhere', latency_s=None,
            basis=f'one human brain emulated at level {name} (Sandberg & Bostrom 2008)')
    t['ancestor_simulation'] = dict(campaign(ANCESTOR_OPS), where='anywhere', latency_s=None,
        product_of_terms=[product * x for x in ANCESTOR_TERMS[3]],
        basis='the mental history of humankind, 1e33-1e36 operations (Bostrom 2003); the product of his figures '
              '(100 billion people, 50 years each, 30 million s a year, 1e14-1e17 operations/s) is 1.5e34-1.5e37')
    t['evolution_anchor'] = dict(campaign([EVOLUTION_ANCHOR_FLOP] * 2), where='anywhere', latency_s=None,
        basis='the computation of all nervous systems through evolution, about 1e41 FLOP (Cotra 2020, via Ho 2022)')
    t['planetary_computer'] = dict(rate([PLANETARY_COMPUTER_OPS] * 2), where='beyond', latency_s=None,
        basis='a computer of planetary mass, 1e42 operations/s (Bostrom 2003, after Bradbury 1999)')
    t['matrioshka_brain'] = dict(rate([MATRIOSHKA[0]] * 2), where='beyond', latency_s=None,
        basis='a multi-layer Matrioshka brain on the Sun\'s 3e26 W (Bradbury 1999)')
    return t


def smallest(tw):
    for s in SCENARIOS_TW:
        if tw <= s:
            return s
    return None


def compare(tasks, eff, sup):
    """Each task against the supply: power at the bounds, light, heat and the smallest scenario that holds it."""
    cmos = eff['cmos_over_today']
    floor = eff['floor_330K_over_today']                          # (low, high) gains over today
    eta = sup['conversion']
    night = sup['night_half_intercepts_pw'] * 1e15
    per = sup['per_tw']
    rows = {}
    for key, t in tasks.items():
        r = dict(where=t['where'], latency_s=t['latency_s'], basis=t['basis'])
        if t['kind'] == 'rate':
            p = t['power_today_w']
            r.update(kind='rate', flop_s=t['flop_s'], power_today_tw=[x / 1e12 for x in p],
                     power_cmos_limit_tw=[x / cmos / 1e12 for x in p],
                     power_landauer_330K_tw=[p[0] / floor[1] / 1e12, p[1] / floor[0] / 1e12],
                     share_of_night_half_light_today=[p[0] / eta[0] / night, p[1] / eta[1] / night],
                     orbit_moon_w_m2_today=[p[0] / 1e12 * per['orbit_moon_w_m2'][0], p[1] / 1e12 * per['orbit_moon_w_m2'][1]],
                     orbit_moon_k_today=[p[0] / 1e12 * per['orbit_moon_k'][0], p[1] / 1e12 * per['orbit_moon_k'][1]],
                     if_on_moon_k_today=[p[0] / 1e12 * per['on_moon_k'][0], p[1] / 1e12 * per['on_moon_k'][1]],
                     smallest_scenario_tw=dict(today=smallest(p[1] / 1e12), cmos_limit=smallest(p[1] / cmos / 1e12)),
                     times_sun=[x / SUN_LUMINOSITY for x in p])
            if t.get('h100_equivalent'):
                r['flop_s_is_h100_equivalent'] = True
        else:
            e = t['energy_today_j']
            r.update(kind='campaign', flop=t['flop'], energy_today_j=e, energy_today_twh=[x / 3.6e15 for x in e],
                     days_at_1_tw_today=[x / 1e12 / DAY for x in e], days_at_1000_tw_today=[x / 1e15 / DAY for x in e],
                     years_at_1000_tw_cmos_limit=[x / cmos / 1e15 / YEAR for x in e],
                     years_at_10000_tw_landauer_330K=[e[0] / floor[1] / 1e16 / YEAR, e[1] / floor[0] / 1e16 / YEAR])
        if 'product_of_terms' in t:
            r['product_of_terms'] = t['product_of_terms']
        rows[key] = r
    return rows


def limits(night_w):
    """Bounds on computing per kilogram and per joule, and the far end of the scale."""
    sun_year = SUN_LUMINOSITY * YEAR
    bit_changes = sun_year / (BOLTZMANN * SCHNEIER_T_K)           # Schneier counts kT per bit change
    erase_256 = 2.0 ** 256 * BOLTZMANN * CMB_K * LN2
    return dict(
        bremermann_bits_per_s_per_kg=SPEED_OF_LIGHT ** 2 / PLANCK,
        margolus_levitin_ops_per_s_per_kg=2 * SPEED_OF_LIGHT ** 2 / (math.pi * HBAR),
        lloyd_laptop_ops_per_s=LLOYD_LAPTOP[0], lloyd_laptop_bits=LLOYD_LAPTOP[1],
        ordinary_matter_ops_per_s_per_kg=ORDINARY_MATTER_OPS_PER_KG, universe_ops=UNIVERSE[0], universe_bits=UNIVERSE[1],
        aestivation_gain=AESTIVATION_GAIN, sun_luminosity_w=SUN_LUMINOSITY, night_half_share_of_sun=night_w / SUN_LUMINOSITY,
        matrioshka_ops_per_j=MATRIOSHKA[0] / MATRIOSHKA[1], matrioshka_power_over_night_half=MATRIOSHKA[1] / night_w,
        matrioshka_power_over_sun=MATRIOSHKA[1] / SUN_LUMINOSITY,
        sun_at_landauer_330K_bits_per_s=SUN_LUMINOSITY / (BOLTZMANN * 330.0 * LN2),
        schneier=dict(sun_year_j=sun_year, bit_changes_at_3_2K=bit_changes, counter_bits=math.log2(bit_changes)),
        erase_2_256_at_2_7K_j=erase_256, erase_2_256_sun_years=erase_256 / sun_year)


def latency(ring_m):
    hill = EARTH_MOON_DISTANCE * (MOON_GM / (3 * EARTH_GM)) ** (1 / 3)
    sel1 = AU * (EARTH_GM / (3 * SUN_GM)) ** (1 / 3)
    c = SPEED_OF_LIGHT
    return dict(moon_array_round_trip_s=2 * ring_m / c, earth_moon_one_way_s=EARTH_MOON_DISTANCE / c,
                earth_moon_round_trip_s=2 * EARTH_MOON_DISTANCE / c, moon_to_earth_moon_l1_l2_one_way_s=hill / c,
                moon_to_sun_earth_l1_l2_one_way_s=[(sel1 - EARTH_MOON_DISTANCE) / c, (sel1 + EARTH_MOON_DISTANCE) / c],
                conversation_gap_s=CONVERSATION_GAP_S)


def links():
    text = [EARTH_POPULATION_PEAK * n * BYTES_PER_TOKEN * 8 for n in TOKENS_PER_S_PER_PERSON]
    video = EARTH_POPULATION_PEAK * VIDEO_BPS
    return dict(tbird_bps=TBIRD_BPS, llcd_bps=LLCD_BPS, text_to_earth_bps=text,
                text_tbird_links=[x / TBIRD_BPS for x in text], video_to_earth_bps=video,
                video_tbird_links=video / TBIRD_BPS)


def moon_side(ctx, sup):
    """Latency-bound computing kept on the Moon, at today's data-centre power per person and a hundred times it."""
    per_person = ctx['data_centres_tw'] * 1e12 / WORLD_POPULATION_2024
    cases = {}
    for people, label in ((METROPOLIS_PEOPLE, 'metropolis'), (1e9, 'per_billion_people')):
        for factor in (1, 100):
            tw = per_person * people * factor / 1e12
            cases[f'{label}_x{factor}'] = dict(tw=tw, k=[tw * k for k in sup['per_tw']['on_moon_k']])
    return dict(data_centre_w_per_person_2024=per_person, cases=cases)


def per_person(eff):
    """A brain-equivalent of computing in watts, and what 0.1-1 kW a person (the population note's range) buys."""
    today, cmos = BRAIN_FLOP_S[1] / eff['today_flop_per_j'], BRAIN_FLOP_S[1] / eff['cmos_limit_flop_per_j']
    return dict(w_per_brain_equivalent_today=today, w_per_brain_equivalent_cmos_limit=cmos,
                ai_services_kw_per_person_today=[n * today / 1e3 for n in BRAIN_EQUIVALENTS_PER_PERSON],
                brain_equivalents_in_0_1_to_1_kw_today=[100 / today, 1000 / today],
                brain_equivalents_in_0_1_to_1_kw_cmos_limit=[100 / cmos, 1000 / cmos])


def context(ah, ledgers, screen):
    hu = screen['heat']['human_use']
    return dict(
        tiles=ah['fleet']['tiles'], night_w=ledgers['energy']['fleet_night_half_intercepts_w'],
        k_per_w_m2=[k / hu['w_m2_per_tw'] for k in hu['k_per_tw']], glow_w_m2=screen['heat']['fleet_infrared_w_m2'],
        data_centres_tw=DATA_CENTRES_TWH_2024 / (YEAR / 3600),
        ai_flop_s_2026=AI_H100_EQUIVALENTS_2026 * H100_FP16_DENSE)


def main(argv=None) -> int:
    ah, ledgers, night, screen = (json.loads(p.read_text()) for p in (HEAT, LEDGERS, NIGHT, SCREEN))
    ctx = context(ah, ledgers, screen)
    eff = efficiency(ah['design']['accelerator_flop_per_J'])
    turned = (night['specular_reflected_lm'] + night['diffuse_reflected_lm']) / LUMENS_PER_PW * 1e15
    sup = supply(ah, ctx['night_w'], night['night_half_tile_area_km2'], turned, eff, ctx)
    wx = weather(eff)
    files = ['provisioning/computing.py']
    month = YEAR / 12
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='provisioning', files={f: digest(ROOT / f) for f in files},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in INPUTS}, constants=constants_used(files)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        today=sig(dict(data_centres_tw_2024=ctx['data_centres_tw'], ai_flop_s_2026=ctx['ai_flop_s_2026'],
                       nvidia_flop_s_2024=NVIDIA_STOCK_FLOP_S_2024, top500_lead_flop_s=TOP500_LEAD_FLOP_S,
                       el_capitan_flop_s_w=list(EL_CAPITAN), google_tokens_per_s=GOOGLE_TOKENS_PER_MONTH / month,
                       tokens_per_s_per_person=GOOGLE_TOKENS_PER_MONTH / month / WORLD_POPULATION_2024,
                       ai_flop_s_per_person=ctx['ai_flop_s_2026'] / WORLD_POPULATION_2024,
                       largest_run_2025_flop=LARGEST_RUN_FLOP_2025, grok4_gwh=GROK4_GWH,
                       frontier_growth_per_year=list(FRONTIER_GROWTH_PER_YEAR), run_2030_flop=RUN_2030_FLOP), 4),
        efficiency=sig(eff, 4), radiators=sig(radiators(ah['design']), 4), supply=sig(sup, 4), weather=sig(wx, 4),
        comparison=sig(compare(demand(ctx, eff, wx), eff, sup), 3), moon_side=sig(moon_side(ctx, sup), 3),
        per_person=sig(per_person(eff), 3),
        limits=sig(limits(ctx['night_w']), 4),
        latency=sig(latency(ah['design']['orbit_radius_km'] * 1e3), 4), links=sig(links(), 3))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for k, r in product['comparison'].items():
        if r['kind'] == 'rate':
            print(f"{k:34s} today {r['power_today_tw']} TW  cmos {r['power_cmos_limit_tw']} TW  "
                  f"light {r['share_of_night_half_light_today']}")
        else:
            print(f"{k:34s} {r['flop']} FLOP  days at 1000 TW {r['days_at_1000_tw_today']}  "
                  f"years at 1000 TW cmos {r['years_at_1000_tw_cmos_limit']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
