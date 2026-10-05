"""The Open Moon's aerosol near the ground, region by region, by day and by night, and what it does to the air's
conductivity and to clouds.

Each region's particles are built from a measured Earth analogue and the lunar factors the literature notes on the data
drive give (aerosol_marine.md, aerosol_biogenic.md, aerosol_dust_fog.md, aerosol_budget.md; each value's section is
named where it is set), every one as a clean, central and loaded case:

- Over clean forest and clean seas the number of particles comes from particles formed high up and brought down by
  storms, not from what the surface emits; surface emissions set their mass, growth and κ (biogenic Sect. 7.0; marine
  Sect. 7). The analogues are the clean wet-season Amazon and the remote ocean by day.
- The Open Moon's low ionization takes away much of the part of that formation that ions drive: 27-40 % of cloud
  nuclei on Earth, with ion-induced rates 7-19 times lower at lunar ionization (budget Sects. 4.3, 7.7).
- Its forests make about half the secondary organic aerosol of an Amazon day (cooler air; biogenic Sect. 7.2), while
  particles live 5-10 times longer and mix through a layer 3-6 times deeper (budget Sects. 7.1-7.2): the accumulation
  mode's mass is about the Amazon's, carried on fewer particles.
- The plains are clean: the seas' sulfur reaches them at a few parts per trillion of sulfur dioxide, against the
  0.7 ppb under which savanna plains form particles daily (biogenic Sect. 3.4; derived in sulfur_dioxide_ppt below).
- Dust comes from bare ground in storm outflows by day, never from ground wet with dew or fog, and the saline lake beds'
  magnesium and calcium chloride bitterns stay deliquesced (dust Sect. 8).
- At night particle formation stops; the air near the ground loses particles to the ground and to fog, while humid
  nights add fungal nanoparticle bursts (biogenic Sect. 7.3; marine night; budget Sect. 7.6).
- In fog the small ions go to the droplets (dust Sect. 8.8).
"""
from __future__ import annotations

import math

from research.studies.open_moon_aerosol.particles import Mode

CASES = ('clean', 'central', 'loaded')       # from the fewest particles to the most
CM3 = 1.0e6

# --- shared lunar factors -------------------------------------------------------------------------------------------
# share of the particle number ion-induced nucleation supplies on Earth (Gordon et al. 2017 via budget Sect. 4.2), and
# the fraction of it left at lunar ionization (rates 7-19 times lower, budget Sect. 4.3)
ION_SHARE = dict(clean=0.40, central=0.33, loaded=0.27)
ION_LEFT = dict(clean=1.0 / 19.0, central=1.0 / 10.0, loaded=1.0 / 7.0)
# the lunar forest's secondary organic aerosol source against an Amazon day (biogenic Sect. 7.2: 33 / 91 / 380 against
# about 180 ug m-2 h-1), and how particle number follows a precursor supply (an exponent: 1 if number follows the
# supply, 0 if not; the formation rate rises about as its square, Kirkby et al. 2016, but sinks shrink with it; choice)
SOA_RATIO = dict(clean=0.18, central=0.5, loaded=2.1)
NUMBER_EXPONENT = dict(clean=1.0, central=0.5, loaded=0.0)
# accumulation-mode mass against the analogue: SOA source times lifetime over mixing depth, (30 d / 5 km) against
# (4 d / 1.5 km), budget Sects. 7.1-7.2 (lifetime x dilution 0.4-1.7)
LIFE_DILUTION = dict(clean=0.4, central=1.0, loaded=1.7)


def ionization_factor(case: str) -> float:
    return 1.0 - ION_SHARE[case] * (1.0 - ION_LEFT[case])


def sulfur_dioxide_ppt(dms_umol_m2_day: float, sea_share: float, land_share: float, mixed_layer_m: float,
                       lifetime_days: float = 1.5, yield_so2: float = 0.7, air_mol_m3: float = 49.0) -> float:
    """Sulfur dioxide (ppt) over land if all the seas' dimethyl sulfide came inland and turned to it at the yield,
    living lifetime_days in the land's mixed layer (derived here)."""
    flux = dms_umol_m2_day * 1e-6 * yield_so2 * sea_share / land_share          # mol per m2 of land per day
    return flux * lifetime_days / mixed_layer_m / air_mol_m3 * 1e12


# --- the regions' particles, day and night ----------------------------------------------------------------------------
# marine day (marine Sect. 7: Heintzenberg et al. 2000, mostly formed aloft; κ 0.4 and 0.5) and sea spray (Gong 2003
# scaled to lunar winds: 0.7 / 2 / 5 cm-3; submicron κ 0.35 with organics); the night's end (Aitken gone, accumulation
# 0.03 / 0.1 / 0.3 of the day's)
MARINE_DAY = dict(clean=(100, 45e-9, 60, 165e-9, 0.7), central=(250, 45e-9, 200, 165e-9, 2.0),
                  loaded=(600, 45e-9, 250, 165e-9, 5.0))
MARINE_NIGHT_END = dict(clean=0.03, central=0.1, loaded=0.3)
# forest day (biogenic Sect. 7.3: Pöhlker et al. 2018 pristine, 2016 clean wet season, Rizzo et al. 2018's wet-season
# 90th percentiles); coarse biological particles 0.3 cm-3 at 1.6 um, κ 0.02
FOREST_DAY = dict(clean=(160, 69e-9, 0.12, 86, 157e-9, 0.18), central=(246, 70e-9, 0.13, 145, 170e-9, 0.21),
                  loaded=(540, 67e-9, 0.14, 480, 172e-9, 0.22))
# humid-night fungal nanoparticle bursts (Lawler et al. 2020, biogenic Sect. 7.3): 455 cm-3 at 30 nm, the share of
# night hours they last (23 % at Oklahoma; the lunar forests' warm decomposing nights a choice above it), κ
FUNGAL_SHARE = dict(wet_land=dict(clean=0.0, central=0.23, loaded=0.5), lakes=dict(clean=0.0, central=0.0, loaded=0.0),
                    fog_desert=dict(clean=0.0, central=0.1, loaded=0.3),
                    polar_dry_land=dict(clean=0.0, central=0.05, loaded=0.2))
FUNGAL_KAPPA = dict(clean=0.1, central=0.2, loaded=0.47)
# dust (dust Sect. 8): number near an active Earth source during emission (11 / 50 / 120 cm-3 above 0.3 um), the lunar
# flux factor at the same friction velocity (1.0 / 1.45 / 7.6), and the share of daytime hours the fog desert's air
# holds that much (outflows over ground left bare by crusts and grazing, with what lingers after them; choice), κ of
# developed-soil dust with some salt (0.003 / 0.03 / 0.3)
DUST_EVENT = dict(clean=11, central=50, loaded=120)
DUST_FLUX = dict(clean=1.0, central=1.45, loaded=7.6)
DUST_DUTY = dict(fog_desert=dict(clean=0.005, central=0.03, loaded=0.15),
                 polar_dry_land=dict(clean=0.0, central=0.0, loaded=0.02),
                 wet_land=dict(clean=0.0, central=0.0, loaded=0.0))
DUST_KAPPA = dict(clean=0.003, central=0.03, loaded=0.3)
# polar dry land: the Arctic's summer (Tunved et al. 2013, budget Sect. 3.3): 50-250 cm-3; no rain over it, calm air
POLAR_DAY = dict(clean=(50, 50e-9, 50, 200e-9), central=(100, 50e-9, 120, 200e-9), loaded=(250, 50e-9, 300, 200e-9))
# the night near the ground (budget Sect. 7.6): no new particles; the surface layer loses particles to the ground at
# 0.3 / 0.05 / 0.01 cm/s (clean / central / loaded; highest over forest), and to fog the activated ones, those above
# ~0.35 um at fog's 0.04 % supersaturation (budget Sect. 7.4; dust Sect. 8.8)
DEPOSITION_M_S = dict(clean=0.003, central=0.0005, loaded=0.0001)
NIGHT_HOURS = 238.0                                                   # dark hours at the equator (ecology register)


def night_mean(ratio_end: float) -> float:
    """The night's mean of a quantity decaying exponentially to ratio_end of its dusk value by the night's end."""
    if ratio_end >= 1.0:
        return 1.0
    k = -math.log(max(ratio_end, 1e-6))
    return (1.0 - math.exp(-k)) / k


def deposition_end(depth_m: float, case: str) -> float:
    return math.exp(-DEPOSITION_M_S[case] * NIGHT_HOURS * 3600.0 / depth_m)


def forest_modes(case: str, phase: str, region: str, night_depth_m: float) -> list:
    na, da, ka, nc, dc, kc = FOREST_DAY[case]
    f_number = ionization_factor(case) * SOA_RATIO[case] ** NUMBER_EXPONENT[case]
    f_mass = SOA_RATIO[case] * LIFE_DILUTION[case]
    if region not in ('wet_land', 'lakes'):           # plains and polar land grow their particles on less vapour
        f_mass *= dict(clean=0.1, central=0.3, loaded=1.0)[case]
    grow = (f_mass / f_number) ** (1.0 / 3.0)
    modes = [Mode(na * f_number * CM3, da, 1.6, ka, 'aitken'),
             Mode(nc * f_number * CM3, dc * grow, 1.5, kc, 'accumulation'),
             Mode(0.3 * CM3, 1.6e-6, 1.6, 0.02, 'biological')]
    if phase == 'night':
        end = deposition_end(night_depth_m, case)
        modes = [m.scaled(night_mean(end)) for m in modes]
        modes.append(Mode(455 * CM3 * FUNGAL_SHARE[region][case], 30e-9, 1.5, FUNGAL_KAPPA[case], 'fungal'))
    return modes


def marine_modes(case: str, phase: str) -> list:
    na, da, nc, dc, spray = MARINE_DAY[case]
    f_number = ionization_factor(case)
    modes = [Mode(na * f_number * CM3, da, 1.5, 0.4, 'aitken'), Mode(nc * f_number * CM3, dc, 1.5, 0.5, 'accumulation'),
             Mode(spray * CM3, 0.25e-6, 2.0, 0.35, 'spray')]
    if phase == 'night':
        modes = [modes[0].scaled(night_mean(0.01)), modes[1].scaled(night_mean(MARINE_NIGHT_END[case])),
                 modes[2].scaled(night_mean(MARINE_NIGHT_END[case]))]
    return modes


def dust_modes(case: str, phase: str, region: str, night_depth_m: float) -> list:
    n = DUST_EVENT[case] * DUST_FLUX[case] * DUST_DUTY[region][case]
    if n <= 0.0:
        return []
    mode = Mode(n * CM3, 0.7e-6, 1.6, DUST_KAPPA[case], 'dust')
    if phase == 'night':                                  # none raised at night; fog takes nearly all of it (dust 8.8)
        mode = mode.scaled(night_mean(deposition_end(night_depth_m, case) * 0.3))
    return [mode]


def polar_modes(case: str, phase: str, night_depth_m: float) -> list:
    na, da, nc, dc = POLAR_DAY[case]
    f_number = ionization_factor(case)
    modes = [Mode(na * f_number * CM3, da, 1.6, 0.3, 'aitken'), Mode(nc * f_number * CM3, dc, 1.5, 0.3, 'accumulation')]
    if phase == 'night':
        end = deposition_end(night_depth_m, case)
        modes = [modes[0].scaled(night_mean(0.01)), modes[1].scaled(night_mean(end))]
        modes.append(Mode(455 * CM3 * FUNGAL_SHARE['polar_dry_land'][case], 30e-9, 1.5, FUNGAL_KAPPA[case], 'fungal'))
    return modes


def modes(region: str, case: str, phase: str, night_depth_m: float) -> list:
    """The particles of a region near the ground, by case and phase."""
    if region in ('seas',):
        return marine_modes(case, phase)
    if region in ('lakes', 'wet_land'):                  # the lakes lie in the forest belt; no forest-floor bursts
        return forest_modes(case, phase, region, night_depth_m)
    if region == 'fog_desert':
        return forest_modes(case, phase, region, night_depth_m) + dust_modes(case, phase, region, night_depth_m)
    if region == 'polar_dry_land':
        return polar_modes(case, phase, night_depth_m) + dust_modes(case, phase, region, night_depth_m)
    raise KeyError(region)


# fog: droplets per cm3 (dust Sect. 8.8: Hammer et al. 2014, 10-150; Weston et al. 2022 in the Namib, 40-80) at 8 um
# radius, about 0.13 g/m3 at the central number; the small ions go to the droplets and the conductivity follows from
# the same ion balance (the dust note's 0.013-0.3 of clear air's as a check)
FOG_DROPLETS = dict(clean=10, central=60, loaded=150)
FOG_RADIUS_M = 8.0e-6
FOG_RATIO_NOTE = dict(clean=0.3, central=0.05, loaded=0.013)
