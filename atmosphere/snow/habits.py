"""Snow morphology and mechanical-load diagnostics, separate from the legacy model.

The habit anchors describe laboratory growth conditions; they are not probabilities
of habits in a lunar cloud. Aggregate mass/area families are sensitivity closures,
not a crystal-shape generator. In particular their extension to decimetre snow is
untested. Source identity and reading extent are in habits_sources.json.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from atmosphere.snow import flakes as sf


@dataclass(frozen=True)
class AggregateGeometry:
    id: str
    mass_coefficient: float
    mass_exponent: float
    area_ratio: float
    capacitance_ratio: float
    evidence: str
    source_ids: tuple[str, ...]

    def mass(self, diameter_m):
        """kg; D is maximum dimension, except the recorded legacy mismatch."""
        return self.mass_coefficient * np.asarray(diameter_m) ** self.mass_exponent

    def size(self, mass_kg):
        """Maximum dimension in m, with the family's stated extrapolation."""
        return (np.asarray(mass_kg) / self.mass_coefficient) ** (1.0 / self.mass_exponent)

    def area(self, diameter_m):
        """Projected area in m²; area_ratio always uses a circumscribing circle."""
        return self.area_ratio * math.pi / 4.0 * np.asarray(diameter_m) ** 2

    def to_dict(self):
        return asdict(self)


def aggregate_area_geometry(habit: str, monomers: int):
    """Seifert et al. (2026), mean horizontal phi and ellipse-area q, eqs 9–17.

    These regressions were fitted to synthetic *unrimed* aggregates with 2–2048
    monomers. Return q*phi to express area relative to a Dmax circle. This product
    of the separate fitted means is a representative geometry, not the mean of
    the jointly distributed product. No extrapolation in monomer number allowed.
    """
    if not isinstance(monomers, (int, np.integer)) or not 2 <= monomers <= 2048:
        raise ValueError('monomers must be an integer from 2 to 2048')
    n = monomers
    if habit == 'needles':
        phi = 0.742 - 0.227 * (2.0 / n) ** 0.95
        q = 0.327 + 0.386 * (2.0 / n) ** 0.75
    elif habit == 'plates':
        phi = 0.736 + 0.141 * (2.0 / n) ** 0.65
        q = (0.900 - 0.305 * (2.0 / n) ** 1.5 if n < 32
             else 0.378 + 0.217 * (32.0 / n) ** 0.45)
    elif habit == 'dendrites':
        phi = 0.732 + 0.172 * (2.0 / n) ** 0.6
        q = (0.404 - 0.0252 * (2.0 / n) ** 2 if n < 16
             else 0.178 + 0.266 * (16.0 / n) ** 0.1)
    else:
        raise ValueError('habit must be needles, plates or dendrites')
    return dict(habit=habit, monomers=int(n), horizontal_aspect_ratio=phi,
                ellipse_area_ratio=q, circle_area_ratio=phi * q,
                source_id='seifert_2026')


_VALIDITY = ('Earth cloud aggregate relation; not validated for 0.15–0.60 m snow. '
             'No riming, crystal habit fractions or geometric scatter represented.')
GEOMETRIES = {
    'legacy_bf95': AggregateGeometry(
        'legacy_bf95', sf.MASS_COEFFICIENT, sf.MASS_EXPONENT, sf.AREA_RATIO, sf.CAPACITANCE,
        'Historical closure: Brown–Francis mean-image dimension used as maximum dimension. '
        'Retained for controlled comparison. ' + _VALIDITY,
        ('hogan_2012', 'westbrook_2008')),
    'bf95_dmax': AggregateGeometry(
        'bf95_dmax', 0.0121, 1.9, sf.AREA_RATIO, sf.CAPACITANCE,
        'Hogan et al. eq. 4 corrects the Brown–Francis dimension convention; '
        'area ratio 0.35 retained to isolate this correction. ' + _VALIDITY,
        ('hogan_2012', 'westbrook_2008')),
}
for _habit in ('needles', 'plates', 'dendrites'):
    _key = f'bf95_dmax_{_habit}_area'
    GEOMETRIES[_key] = AggregateGeometry(
        _key, 0.0121, 1.9, aggregate_area_geometry(_habit, 1024)['circle_area_ratio'], sf.CAPACITANCE,
        'Fixed representative area sensitivity: q(1024)*phi(1024) from synthetic '
        f'{_habit} aggregates. Holding this ratio at all sizes and combining it '
        'with Brown–Francis is a modelling scenario, not a habit-specific mass fit. ' + _VALIDITY,
        ('hogan_2012', 'seifert_2026', 'westbrook_2008'))


# Temperature anchors from Libbrecht 2012, introduction and Fig. 1. They are
# examples of favored forms, not boundaries or an exhaustive habit catalogue.
HABIT_ANCHORS = (
    (-15.0, 'thin plates; branching and dendrites at sufficiently high supersaturation'),
    (-7.0, 'thin-walled hollow columns'),
    (-5.0, 'slender columns and needles'),
    (-2.0, 'thin plate-like forms'),
)


def habit_ladder(column):
    """Map known growth-temperature examples into this hypothetical air column.

    Supercooled liquid saturation is the existing snow model's imposed boundary
    condition. The lapse gives heights, but does not establish cloud occupancy,
    droplet abundance, nucleation or the fraction of snow in any crystal habit.
    """
    rows = []
    for temperature_c, habit in HABIT_ANCHORS:
        t = sf.T0 + temperature_c
        z = column.height_of(t)
        s = column.state(z)
        rows.append(dict(temperature_c=temperature_c, height_m=z, pressure_pa=s['p'],
                         water_saturation_ice_supersaturation=sf.liquid_saturation_pressure(t)
                         / sf.ice_saturation_pressure(t) - 1.0,
                         growth_form_context=habit, source_id='libbrecht_2012'))
    return rows


def humidity_regime(temperature_c: float, relative_humidity_water: float):
    """Thermodynamic growth context, not a crystal-habit prediction.

    The same subzero air can be unsaturated over water yet supersaturated over
    ice. Matching a temperature alone therefore cannot identify a snow habit.
    """
    if not -40.0 <= temperature_c <= 0.0 or not math.isfinite(relative_humidity_water):
        raise ValueError('temperature must be -40..0 C and relative humidity finite')
    if relative_humidity_water < 0:
        raise ValueError('relative humidity must be nonnegative')
    t = sf.T0 + temperature_c
    ratio = (relative_humidity_water * sf.liquid_saturation_pressure(t)
             / sf.ice_saturation_pressure(t))
    if ratio < 1.0:
        regime = 'ice sublimation'
    elif relative_humidity_water < 1.0:
        regime = 'vapour growth possible; subsaturated over liquid water'
    else:
        regime = 'vapour growth possible; liquid saturation reached or exceeded'
    return dict(supersaturation_ice=ratio - 1.0, regime=regime)


def differential_load(diameter_m: float, air_density_kg_m3: float, velocity_difference_m_s: float,
                      area_ratio: float, drag_coefficient: float = 1.0):
    """Reference nonuniform aerodynamic load; does not decide fracture.

    q=0.5*rho*du² and F=Cd*q*A. du is an imposed velocity difference across a
    flake, not its shared translation with the wind. Treating that pressure as
    acting across the entire projected area is a diagnostic loading convention.
    Real bending and tensile loads depend on drag distribution, porosity, local
    slip and connectivity. Cd=1 is an explicit reference, not a measured value.
    """
    if not all(math.isfinite(x) for x in (diameter_m, air_density_kg_m3, velocity_difference_m_s,
                                         area_ratio, drag_coefficient)):
        raise ValueError('all load inputs must be finite')
    if diameter_m <= 0 or air_density_kg_m3 <= 0 or velocity_difference_m_s < 0:
        raise ValueError('diameter and density must be positive; velocity difference nonnegative')
    if not 0 < area_ratio <= 1 or drag_coefficient <= 0:
        raise ValueError('area ratio must be 0..1 and drag coefficient positive')
    pressure = 0.5 * air_density_kg_m3 * velocity_difference_m_s ** 2
    area = area_ratio * math.pi / 4.0 * diameter_m ** 2
    force = drag_coefficient * pressure * area
    return dict(dynamic_pressure_pa=pressure, projected_area_m2=area,
                reference_differential_force_n=force,
                evidence='Reference load only; no measured aggregate failure strength assumed.')


def critical_velocity_difference(force_capacity_n: float, diameter_m: float,
                                 air_density_kg_m3: float, area_ratio: float,
                                 drag_coefficient: float = 1.0):
    """du at which the reference load equals an explicitly supplied force capacity.

    Capacity belongs to the aggregate's load-bearing contact network. This module
    supplies no preferred capacity and imports no snowpack or solid-ice strength.
    """
    if not math.isfinite(force_capacity_n) or force_capacity_n < 0:
        raise ValueError('force capacity must be finite and nonnegative')
    reference = differential_load(diameter_m, air_density_kg_m3, 1.0, area_ratio, drag_coefficient)
    return math.sqrt(force_capacity_n / reference['reference_differential_force_n'])
