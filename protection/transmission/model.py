"""Ultraviolet transmission the optical shield's swarm can hold: seams, micrometeoroid holes, missing cells, defects.

    python -m protection.transmission.model   # writes protection/transmission/results/uv_transmission.json

Below 175 nm the intact titania-silica film passes at most 1.6e-43 of direct
sunlight (spectra/stack_short_wave.json, 10-175 nm). The shield's ultraviolet
transmission is therefore the uncovered share of the aperture, averaged over
time and area: the report's f + (1 - f) T_material (section 8) with T_material
negligible. Four terms make up f.

Seams. Square cells of side a (10 km in the September design) meet along 2/a of
seam per unit area. Neighbours overlap by w and sit at slightly different
distances, h apart, so sunlight slanting at up to the Sun's angular radius
crosses h tan(theta_sun) of the overlap; the overlap left is
w_eff = w - h tan(theta_sun). The neighbours' relative position normal to a seam
scatters with standard deviation sigma, and a gap opens where their separation
exceeds w_eff:

    f_seam = (2/a) [sigma phi(w_eff/sigma) - w_eff (1 - Phi(w_eff/sigma))].

Micrometeoroid holes. The Grun et al. (1985) interplanetary flux at 1 AU gives
the impacts per unit area on each face (both faces count). A particle of density
2.5 g/cm^3 and diameter d at least the film's 11 um perforates a hole of
diameter k d. Penetration experiments in aluminium targets with glass
projectiles of 50 um to 3.2 mm (Horz and colleagues) found the hole shrinking
from about 4 d in massive targets to d in very thin foils. The particles that
dominate the holed area here are near 60 um, five times the film's thickness, so
k is taken over 1-4; the source is known from its abstract, and hypervelocity
tests on the oxide film itself are to come. Holes accumulate until a cell is
replaced or patched, every L years, so the mean holed share is r k^2 L / 2 for a
holing rate r at k = 1.

Missing cells. A cell that fails leaves its area open until a spare or an
overlay covers it: f_missing = lambda * t_cover, for failures per cell-year
lambda and a time to cover t_cover. Planned replacement that places the new cell
before removing the old adds nothing.

Defects. Pinholes from manufacture, a production specification f_defect.

Light through an opening spreads over the Sun's blur before it reaches the
Moon: a missing cell at 78,000 km lights a patch of radius a/2 plus the Sun's
angular radius times the distance, at the share of full sunlight that the
cell's solid angle is of the Sun's.

The three design levels combine choices for each term. They are assumptions to
test, set out so that each term's requirement can be read on its own.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
PROTECTION = HERE.parent
ROOT = PROTECTION.parent
SCHEMA = 'terluna.protection.uv-transmission/1'
YEAR_S = 365.25 * 86400.0
SUN_ANGULAR_RADIUS = 0.00465          # radians at the Moon (report section 6)
CELL_SIDE_M = 10e3                    # the September design's 10 x 10 km cells (report section 8)
SHIELD_DISTANCE_M = 78000e3           # the September design's distance sunward of the Moon (report section 6)
FILM_THICKNESS_M = 11e-6              # 1 um titania and 10 um silica
METEOROID_DENSITY_KG_M3 = 2500.0      # the Grun model's assumed density
HOLE_FACTORS = (1.0, 2.0, 4.0)

LEVELS = {
    'tight': dict(sigma_m=1.0, overlap_m=50.0, layer_gap_m=100.0, patch_years=10.0,
                  failures_per_cell_year=0.01, cover_days=1.0, defect_fraction=1e-6),
    'standard': dict(sigma_m=3.0, overlap_m=50.0, layer_gap_m=100.0, patch_years=20.0,
                     failures_per_cell_year=0.01, cover_days=7.0, defect_fraction=1e-6),
    'relaxed': dict(sigma_m=10.0, overlap_m=20.0, layer_gap_m=100.0, patch_years=100.0,
                    failures_per_cell_year=0.03, cover_days=30.0, defect_fraction=1e-5),
}


def grun_flux(mass_g):
    """Grun et al. (1985) cumulative interplanetary flux at 1 AU: impacts per m^2 per s of particles above mass_g."""
    m = np.asarray(mass_g, dtype=float)
    return ((2.2e3 * m ** 0.306 + 15.0) ** -4.38 + 1.3e-9 * (m + 1e11 * m ** 2 + 1e27 * m ** 4) ** -0.36
            + 1.3e-16 * (m + 1e6 * m ** 2) ** -0.85)


def particle_diameter_m(mass_g):
    return (6.0 * np.asarray(mass_g) * 1e-3 / (math.pi * METEOROID_DENSITY_KG_M3)) ** (1.0 / 3.0)


def holing_rate_per_year(min_diameter_m=FILM_THICKNESS_M, faces=2):
    """Area share perforated per year with holes the particles' own size (k = 1)."""
    m = np.logspace(-15, 2, 20000)
    d = particle_diameter_m(m)
    per_mass = -np.gradient(grun_flux(m), m)
    keep = d >= min_diameter_m
    return faces * float(np.trapezoid((math.pi / 4 * d ** 2 * per_mass)[keep], m[keep])) * YEAR_S


def dominant_particle_diameter_m():
    m = np.logspace(-15, 2, 20000)
    d = particle_diameter_m(m)
    weight = math.pi / 4 * d ** 2 * -np.gradient(grun_flux(m), m) * m     # holed area per logarithmic mass interval
    return float(d[np.argmax(weight)])


def seam_fraction(sigma_m, overlap_m, layer_gap_m, cell_side_m=CELL_SIDE_M):
    w = overlap_m - layer_gap_m * math.tan(SUN_ANGULAR_RADIUS)
    z = w / sigma_m
    mean_gap = sigma_m * norm.pdf(z) - w * norm.sf(z)
    return 2.0 / cell_side_m * max(mean_gap, 0.0)


def hole_fraction(patch_years, hole_factor, rate=None):
    rate = holing_rate_per_year() if rate is None else rate
    return rate * hole_factor ** 2 * patch_years / 2.0


def missing_fraction(failures_per_cell_year, cover_days):
    return failures_per_cell_year * cover_days / 365.25


def level(params, rate):
    terms = dict(seams=seam_fraction(params['sigma_m'], params['overlap_m'], params['layer_gap_m']),
                 missing_cells=missing_fraction(params['failures_per_cell_year'], params['cover_days']),
                 defects=params['defect_fraction'])
    holes = {f'{k:g}': hole_fraction(params['patch_years'], k, rate) for k in HOLE_FACTORS}
    fixed = sum(terms.values())
    return dict(parameters=params, terms=dict(terms, micrometeoroid_holes_by_hole_factor=holes),
                transmission_by_hole_factor={k: fixed + v for k, v in holes.items()},
                overlap_film_mass_share=2.0 * params['overlap_m'] / CELL_SIDE_M)


def missing_cell_patch(cell_side_m=CELL_SIDE_M, distance_m=SHIELD_DISTANCE_M):
    """The patch a single missing cell lights on the Moon, and its share of full sunlight."""
    share = (cell_side_m ** 2 / distance_m ** 2) / (math.pi * SUN_ANGULAR_RADIUS ** 2)
    return dict(radius_km=(cell_side_m / 2 + SUN_ANGULAR_RADIUS * distance_m) / 1e3, share_of_full_sunlight=share)


def film_transmission_max():
    d = json.loads((PROTECTION / 'spectra' / 'stack_short_wave.json').read_text())
    w, t = np.array(d['wavelength_nm']), np.array(d['transmission'])
    keep = (w >= 10.0) & (w <= 175.0)
    return float(t[keep].max())


def sweeps(rate):
    """Each term over its own lever, the other choices at the standard level."""
    s = LEVELS['standard']
    return dict(
        seams_by_overlap_m={f'{w:g}': seam_fraction(s['sigma_m'], w, s['layer_gap_m']) for w in (0, 5, 10, 20, 50, 100)},
        seams_butt_joint_by_sigma_m={f'{x:g}': seam_fraction(x, 0.0, 0.0) for x in (1, 3, 10, 30)},
        holes_by_patch_years={f'{L:g}': {f'{k:g}': hole_fraction(L, k, rate) for k in HOLE_FACTORS}
                              for L in (10, 20, 100, 1000)},
        missing_by_cover_days={f'{t:g}': {f'{x:g}': missing_fraction(x, t) for x in (0.001, 0.01, 0.03)}
                               for t in (1, 7, 30)})


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    rate = holing_rate_per_year()
    levels = {name: level(p, rate) for name, p in LEVELS.items()}
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='protection', files={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16] for p in (
            'protection/transmission/model.py', 'protection/spectra/stack_short_wave.json')}),
        evidence=('The time- and area-averaged share of the aperture left uncovered, for three sets of design choices: '
                  'seams from a Gaussian relative-position error against an overlap, micrometeoroid holes from the '
                  'Grun et al. (1985) flux with holes 1-4 times the particle diameter, cells missing between failure '
                  'and cover, and a manufacturing defect specification. The choices are assumptions to test; no '
                  'formation control, impact test or failure record for such cells exists yet.'),
        reading_rule=('levels[name].transmission_by_hole_factor["k"] is the fraction of direct sunlight below 175 nm that '
                      'reaches the protected region, to set against the loss response\'s allowed transmission '
                      '(its field allowed_leak_fraction). terms gives each contribution; sweeps give each term over '
                      'its lever. The intact film adds film_transmission_max_10_175_nm at most.'),
        film_transmission_max_10_175_nm=film_transmission_max(),
        holing_rate_per_year_hole_factor_1=rate,
        dominant_particle_diameter_m=dominant_particle_diameter_m(),
        cell_side_m=CELL_SIDE_M, missing_cell_patch=missing_cell_patch(), levels=levels, sweeps=sweeps(rate))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'uv_transmission.json').write_text(json.dumps(product, indent=2) + '\n')
    for name, lv in levels.items():
        t = lv['transmission_by_hole_factor']
        terms = lv['terms']
        print(f"{name:8s}: transmission {t['1']:.2e} - {t['4']:.2e}  (seams {terms['seams']:.1e}, holes "
              f"{terms['micrometeoroid_holes_by_hole_factor']['1']:.1e}-{terms['micrometeoroid_holes_by_hole_factor']['4']:.1e}, "
              f"missing cells {terms['missing_cells']:.1e}, defects {terms['defects']:.0e})")
    print(f"holing rate {rate:.2e} per year at k = 1; dominant particle {dominant_particle_diameter_m()*1e6:.0f} um")
    return 0


if __name__ == '__main__':
    sys.exit(main())
