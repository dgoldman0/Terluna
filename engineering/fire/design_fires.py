"""Design fires: EN 1991-1-2 Annex E's occupancies and growth classes, and what is measured of how fires and materials
change at lunar gravity and in air with less oxygen.

A design fire grows as Q = 1 MW (t / t_alpha)^2 until it reaches its peak heat release per unit of floor burning
(EN 1991-1-2 E.4). Annex E gives t_alpha and that peak for common occupancies, and their fire load densities.

Two effects on the fire itself are measured, on Earth or in short low-gravity tests:

- Less oxygen. Peatross and Beyler (1997) found pool fires' burning rate falls linearly with the oxygen around them,
  m/m0 = 10 X_O2 - 1.1, reaching zero near 11% (at 1 atm and Earth gravity). Normalised to Earth's air, the Open
  Moon's 17.5% gives 0.65.
- Lunar gravity. Drop-tower centrifuge tests (Ferkul and Olson 2011) found three spacecraft materials burned in
  oxygen 2.2 to 5.9 points lower at lunar gravity than their limits in NASA's 1 g upward test, and later tests found
  acrylic most flammable near lunar gravity (Olson et al. 2024). In partial-gravity flights, upward flame spread over
  very thin fuels was roughly proportional to gravity (Feier et al. 2002). No fire of room size has burned at
  partial gravity.

So a fire at lunar gravity is likely slower to grow and smaller per unit area than its Earth design fire, while more
materials sustain one. `lunar_estimate` gives the slower bracket: one growth class slower (t_alpha doubled) and the
peak scaled by the burning-rate factor. Designs keep the Earth fire as the other bracket until full-scale tests at
0.16 g exist.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

GROWTH_TIME_S = dict(slow=600.0, medium=300.0, fast=150.0, ultrafast=75.0)   # time to reach 1 MW (EN 1991-1-2 E.4)
EARTH_OXYGEN = 0.2095


@dataclass(frozen=True)
class DesignFire:
    growth_time_s: float          # t_alpha: the time to reach 1 MW
    peak_kw_m2: float             # heat release per m2 of floor burning, fuel-controlled
    fire_load_mj_m2: float        # 80% fractile fire load density

    @property
    def alpha_kw_s2(self) -> float:
        return 1000.0 / self.growth_time_s ** 2

    def hrr_kw(self, t_s: float, cap_kw: float = math.inf) -> float:
        return min(self.alpha_kw_s2 * t_s * t_s, cap_kw)


@dataclass(frozen=True)
class Occupancy:
    name: str
    fire_load_mean_mj_m2: float
    fire_load_80_mj_m2: float
    growth: str
    peak_kw_m2: float

    def earth_fire(self) -> DesignFire:
        return DesignFire(GROWTH_TIME_S[self.growth], self.peak_kw_m2, self.fire_load_80_mj_m2)


# EN 1991-1-2:2002 Annex E, Tables E.4 (fire load densities, MJ/m2: mean and 80% fractile) and E.5 (growth and RHRf).
EN1991_ANNEX_E = {o.name: o for o in (
    Occupancy('dwelling', 780, 948, 'medium', 250),
    Occupancy('hospital (room)', 230, 280, 'medium', 250),
    Occupancy('hotel (room)', 310, 377, 'medium', 250),
    Occupancy('library', 1500, 1824, 'fast', 500),
    Occupancy('office', 420, 511, 'medium', 250),
    Occupancy('classroom of a school', 285, 347, 'medium', 250),
    Occupancy('shopping centre', 600, 730, 'fast', 250),
    Occupancy('theatre (cinema)', 300, 365, 'fast', 500),
    Occupancy('transport (public space)', 100, 122, 'slow', 250),
)}


def burning_rate_factor(oxygen_mole_fraction: float) -> float:
    """Burning rate relative to Earth's air, from Peatross and Beyler's m/m0 = 10 X_O2 - 1.1 (zero below 11%)."""
    return max(0.0, 10.0 * oxygen_mole_fraction - 1.1) / (10.0 * EARTH_OXYGEN - 1.1)


def lunar_estimate(fire: DesignFire, oxygen_mole_fraction: float) -> DesignFire:
    """The slower bracket for a fire at lunar gravity in air with this oxygen: one growth class slower and the peak
    heat release per unit area scaled by the burning-rate factor. The fire load is unchanged."""
    return DesignFire(2.0 * fire.growth_time_s, fire.peak_kw_m2 * burning_rate_factor(oxygen_mole_fraction),
                      fire.fire_load_mj_m2)


# Upward limiting oxygen index (ULOI) and maximum oxygen concentration for upward flame spread (MOC), in mole % O2,
# at Earth and lunar gravity: Ferkul and Olson (2011), drop-tower centrifuge, as tabled by Miller et al. in their
# decadal-survey paper on spacecraft materials fire safety (Table I). The low-gravity limits use a shorter extinction
# criterion than the 1 g test, forced by the drop time.
LUNAR_OXYGEN_LIMITS = {
    'Mylar G film (70 kPa)': dict(uloi_1g=21.2, moc_1g=20.0, uloi_lunar=15.6, moc_lunar=14.1),
    'Ultem 1000 (70 kPa)': dict(uloi_1g=23.5, moc_1g=23.0, uloi_lunar=21.0, moc_lunar=19.9),
    'Nomex HT90-40 fabric (101 kPa)': dict(uloi_1g=23.5, moc_1g=22.1, uloi_lunar=21.0, moc_lunar=19.9),
}


def lunar_oxygen_shifts() -> dict:
    """How many points of oxygen each material's limits fall at lunar gravity."""
    return {k: dict(uloi=round(v['uloi_1g'] - v['uloi_lunar'], 2), moc=round(v['moc_1g'] - v['moc_lunar'], 2))
            for k, v in LUNAR_OXYGEN_LIMITS.items()}


def materials_test_oxygen_percent(design_oxygen_mole_fraction: float, margin_points: float, step: float = 0.5) -> float:
    """The oxygen (mole %) at which a material must stop burning in the 1 g upward test so that it will not burn at
    lunar gravity in the design air: the design fraction plus the margin, rounded up to the step."""
    return math.ceil((100.0 * design_oxygen_mole_fraction + margin_points) / step) * step
