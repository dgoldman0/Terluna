"""The rain each audition hears, read from the atmosphere domain's arrivals product.

A scene is the rain arriving on one surface in one target climate, on the Moon or on Earth, at a rain rate on the
ground: the drops of a bounded cloud-base spectrum that survive the fall (atmosphere/rain/results/
rain_arrivals.json), each class with its size, its speed at impact and its number per square metre per second; for
grass and leaves a share of the same water arriving as drips released at rest from leaves; and, as a sensitivity, a
share arriving as drops of 8-12 mm. The air at 2 m, which sets the sound's speed and its absorption, comes from the
same product, and the sound speed from the thunder model (atmosphere/electricity/thunder.py).

The Moon's rain falls at the GCM's rates for its target (atmosphere/rain/results/rain_rates.json): a median rainy
6.5 hours on the land and a heavy one. The Earth texture stands for Earth rain of 5 mm/h, the rain work's assumption
for its recordings, whose own rain went unmeasured. With independent drops the arriving spectrum keeps its shape at
any rate, and the arrivals product gives the rain leaving the cloud base that a rate at the ground takes. Rain falls
vertically onto a horizontal surface in still air, so the impact speed is the speed normal to the surface.
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from atmosphere.rain.distributions import MM_H_PER_M_S  # noqa: E402
from atmosphere.rain.drops import WATER_DENSITY, WATER_MOLAR_MASS, surface_tension  # noqa: E402
from atmosphere.electricity import thunder  # noqa: E402

ARRIVALS = ROOT / 'atmosphere' / 'rain' / 'results' / 'rain_arrivals.json'
ARRIVALS_SCHEMA = 'terluna.atmosphere.rain-arrivals/1'
RATES = ROOT / 'atmosphere' / 'rain' / 'results' / 'rain_rates.json'
RATES_SCHEMA = 'terluna.atmosphere.rain-rates/1'
TARGETS = ('hot_world', 'resort', 'lock')
WORLDS = ('moon', 'earth')
EARTH_REFERENCE_MM_H = 5.0          # the rain the Earth texture stands for
# The grass and leaves' assumed water route, as in the rain work: a quarter of the rain reaches the surface as 3 mm
# drips released at rest 1 m above it. Interception, storage and release wait on biosphere/canopy.
DEFAULT_DRIPS = dict(fraction=0.25, diameter_mm=3.0, height_m=1.0)


def load_arrivals(path=ARRIVALS) -> dict:
    product = json.loads(Path(path).read_text())
    if product.get('schema') != ARRIVALS_SCHEMA:
        raise ValueError(f'Expected {ARRIVALS_SCHEMA}, found {product.get("schema")}')
    return product


def load_rates(path=RATES) -> dict:
    product = json.loads(Path(path).read_text())
    if product.get('schema') != RATES_SCHEMA:
        raise ValueError(f'Expected {RATES_SCHEMA}, found {product.get("schema")}')
    return product


def gcm_rates(rates, target) -> dict:
    """A target's rain rates on the land (mm/h) from the GCM's 6.5-hour means: the median wet mean, a median rainy
    6.5 hours, and the 99th percentile of the wet means, a heavy one."""
    wet = rates['targets'][target]['six_and_a_half_hour_means']['wet']
    return dict(median=wet['p50'], heavy=wet['p99'])


@dataclass(frozen=True)
class Population:
    """Impacts on a horizontal surface: diameters (m), number flux (per m2 per s) and speed normal to it (m/s)."""
    label: str
    diameter_m: np.ndarray
    number_flux_m2_s: np.ndarray
    normal_speed_m_s: np.ndarray
    liquid_density_kg_m3: float = WATER_DENSITY

    def __post_init__(self):
        d, j, u = (np.asarray(x, dtype=float) for x in (self.diameter_m, self.number_flux_m2_s, self.normal_speed_m_s))
        if d.ndim != 1 or not d.size or j.shape != d.shape or u.shape != d.shape:
            raise ValueError('Diameter, flux and speed must be matching non-empty vectors')
        if not (np.all(np.isfinite(d)) and np.all(np.isfinite(j)) and np.all(np.isfinite(u))):
            raise ValueError('Non-finite impact input')
        if np.any(d <= 0) or np.any(j < 0) or np.any((j > 0) & (u <= 0)):
            raise ValueError('Diameters must be positive, fluxes non-negative and arriving drops moving')
        for name, value in (('diameter_m', d), ('number_flux_m2_s', j), ('normal_speed_m_s', u)):
            object.__setattr__(self, name, value)

    @property
    def rain_mm_h(self) -> float:
        return float(np.dot(self.number_flux_m2_s, math.pi / 6 * self.diameter_m**3) * MM_H_PER_M_S)

    @property
    def arrivals_m2_s(self) -> float:
        return float(self.number_flux_m2_s.sum())

    def describe(self) -> dict:
        j, d, u = self.number_flux_m2_s, self.diameter_m, self.normal_speed_m_s
        mass = self.liquid_density_kg_m3 * math.pi / 6 * d**3
        return dict(label=self.label, rain_mm_h=self.rain_mm_h, arrivals_m2_s=self.arrivals_m2_s,
                    mean_diameter_mm_by_number=float(np.dot(j, d) / j.sum() * 1e3) if j.sum() else None,
                    mean_speed_m_s_by_number=float(np.dot(j, u) / j.sum()) if j.sum() else None,
                    momentum_flux_pa=float(np.dot(j, mass * u)), kinetic_power_w_m2=float(np.dot(j, 0.5 * mass * u**2)))


def _side(product, target, world):
    if target not in TARGETS or world not in WORLDS:
        raise ValueError(f'target must be one of {TARGETS} and world one of {WORLDS}')
    return product['targets'][target][world]


def direct_rain(product, target, world, rain_mm_h, spectrum=None) -> Population:
    """The drops of a bounded cloud-base spectrum arriving at rain_mm_h on the ground (the central one by default)."""
    side = _side(product, target, world)
    key = spectrum or (product['central_spectrum']['dm_mm'], product['central_spectrum']['shape_mu'])
    row = next(s for s in side['spectra'] if (s['dm_mm'], s['shape_mu']) == tuple(key))
    if row['ground_flux_per_mm_h'] is None:
        raise ValueError(f'The {key} spectrum delivers no water to the ground in {target} on {world}')
    flux = np.asarray(row['ground_flux_per_mm_h']) * float(rain_mm_h)
    keep = flux > 0
    size = np.array([np.nan if x is None else x for x in side['classes']['surface_diameter_mm']]) * 1e-3
    speed = np.array([np.nan if x is None else x for x in side['classes']['impact_speed_m_s']])
    return Population(f'{world} direct rain, spectrum {key[0]:g} mm / shape {key[1]:g}', size[keep], flux[keep],
                      speed[keep])


def drips(product, target, world, rain_mm_h, diameter_mm=3.0, height_m=1.0) -> Population:
    """Drips of one size released at rest from one height, delivering rain_mm_h at the ground."""
    side = _side(product, target, world)
    row = next((r for r in side['drips'] if r['release_diameter_mm'] == diameter_mm
                and r['release_height_m'] == height_m), None)
    if row is None:
        raise ValueError(f'The product has no {diameter_mm} mm drip from {height_m} m')
    d = row['surface_diameter_mm'] * 1e-3
    flux = float(rain_mm_h) / MM_H_PER_M_S / (math.pi / 6 * d**3)
    return Population(f'{world} drips, {diameter_mm:g} mm from {height_m:g} m', np.array([d]), np.array([flux]),
                      np.array([row['impact_speed_m_s']]))


def tail(product, target, world, rain_mm_h) -> Population:
    """Drops of 8-12 mm leaving the cloud base, their water spread evenly over the sizes, delivering rain_mm_h."""
    rows = [r for r in _side(product, target, world)['tail'] if r['status'] == 'surface']
    if not rows:
        raise ValueError(f'No 8-12 mm drop reaches the ground in {target} on {world}')
    d = np.array([r['surface_diameter_mm'] for r in rows]) * 1e-3
    flux = float(rain_mm_h) / len(rows) / MM_H_PER_M_S / (math.pi / 6 * d**3)
    return Population(f'{world} drops of 8-12 mm', d, flux, np.array([r['impact_speed_m_s'] for r in rows]))


def air(product, target, world) -> dict:
    """The air at 2 m for the sound: temperature, pressure, humidity, oxygen and nitrogen, its sound speed, and
    the surface tension of water at its temperature."""
    a = _side(product, target, world)['air']
    x_v = a['vapour_mole_fraction']
    mixing_ratio = x_v * WATER_MOLAR_MASS / ((1.0 - x_v) * a['dry_molar_mass_kg_mol'])
    # The thunder model's oxygen and nitrogen are fractions of the dry air, as the thunder study passes them.
    return dict(temperature_k=a['temperature_k'], pressure_pa=a['pressure_pa'],
                relative_humidity=a['relative_humidity'], vapour_percent=100.0 * x_v,
                x_o2=a['dry_fractions']['O2'], x_n2=a['dry_fractions']['N2'],
                dry_molar_mass_kg_mol=a['dry_molar_mass_kg_mol'], gravity_m_s2=a['gravity_m_s2'],
                surface_tension_n_m=surface_tension(a['temperature_k']),
                sound_speed_m_s=float(thunder.sound_speed(a['temperature_k'], mixing_ratio,
                                                          a['dry_molar_mass_kg_mol'])))


@dataclass(frozen=True)
class Scene:
    """The rain on one surface in one place: direct rain, drips and the optional tail, sharing one water budget."""
    target: str
    world: str
    surface: str
    rain_mm_h: float
    direct: Population
    drips: Population | None
    tail: Population | None
    air: dict
    settings: dict = field(default_factory=dict)

    def components(self):
        return [p for p in (self.direct, self.drips, self.tail) if p is not None]

    def ledger(self) -> dict:
        parts = {p.label: p.rain_mm_h for p in self.components()}
        return dict(prescribed_ground_rain_mm_h=self.rain_mm_h, components_mm_h=parts,
                    balance_mm_h=sum(parts.values()) - self.rain_mm_h)

    def describe(self) -> dict:
        return dict(target=self.target, world=self.world, surface=self.surface, rain_mm_h=self.rain_mm_h,
                    air=self.air, settings=self.settings, ledger=self.ledger(),
                    components=[p.describe() for p in self.components()])


def scene(product, target, world, surface, rain_mm_h=EARTH_REFERENCE_MM_H, spectrum=None, drip_fraction=None,
          drip_diameter_mm=DEFAULT_DRIPS['diameter_mm'], drip_height_m=DEFAULT_DRIPS['height_m'],
          tail_fraction=0.0) -> Scene:
    """The scene for one surface. Grass and leaves route drip_fraction of the water (a quarter by default) through
    drips; tail_fraction of the water, at most 5%, arrives as 8-12 mm drops where they reach the ground."""
    rain = float(rain_mm_h)
    if not (math.isfinite(rain) and 0 < rain <= 100):
        raise ValueError('Use a ground rain rate between 0 and 100 mm/h')
    if surface not in ('concrete', 'puddles', 'grass-leaves'):
        raise ValueError('surface must be concrete, puddles or grass-leaves')
    if not 0 <= tail_fraction <= 0.05:
        raise ValueError('The tail carries 0-5% of the water')
    drip_share = (DEFAULT_DRIPS['fraction'] if drip_fraction is None else float(drip_fraction)) \
        if surface == 'grass-leaves' else 0.0
    if not 0 <= drip_share < 1:
        raise ValueError('The drip share lies in [0, 1): some rain falls directly')
    remainder = rain * (1.0 - tail_fraction)
    direct = direct_rain(product, target, world, remainder * (1.0 - drip_share), spectrum)
    dripping = drips(product, target, world, remainder * drip_share, drip_diameter_mm, drip_height_m) \
        if drip_share > 0 else None
    large = tail(product, target, world, rain * tail_fraction) if tail_fraction > 0 else None
    key = tuple(spectrum or (product['central_spectrum']['dm_mm'], product['central_spectrum']['shape_mu']))
    row = next(s for s in _side(product, target, world)['spectra'] if (s['dm_mm'], s['shape_mu']) == key)
    settings = dict(spectrum=list(key),
                    cloud_base_rain_mm_h=direct.rain_mm_h * row['cloud_base_rain_per_ground_rain'],
                    drip_fraction=drip_share, drip_diameter_mm=drip_diameter_mm if drip_share else None,
                    drip_height_m=drip_height_m if drip_share else None, tail_fraction=tail_fraction)
    return Scene(target, world, surface, rain, direct, dripping, large, air(product, target, world), settings)
