"""How much the lunar cycle's photosynthesis and growth follow the light: the plant model with the sunlight dimmed.

    python -m biosphere.canopy.light_response    # writes results/light_response.json (a few minutes, one thread)

The plant model's stand (plant.py: leaf area index 5, the chosen climate's land air temperatures, the solved sky's
low Sun and twilight) is followed through one solar cycle with the Moon's light, direct, diffuse and twilight alike,
scaled by 1 - d for dimmings d of 1, 5 and 10%. The temperatures stay the chosen climate's, since the dimming exists
to keep them, and stem and root upkeep keeps its Earth calibration. The response is the share of photosynthesis or
growth lost for each share of light taken, the elasticity: 1 where growth follows the light, 0 where it ignores it.

It is the response with water, nutrients and sinks ample. Where they limit growth, the light matters less; the light
limitation study (research/studies/light_limitation) maps where that is.
"""
from __future__ import annotations
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import numpy as np

from biosphere.canopy import model as canopy
from biosphere.canopy import plant

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.biosphere.canopy-light-response/1'
DIMMINGS = (0.01, 0.05, 0.10)
STRATEGIES = ('earth_like', 'idle_quarter')         # the two evergreen ways through the night the register adopts
SIDE = 'near'                                       # the far side differs by 2-3% (plant.py)
LEAF_C = 22.0                                       # the canopy model's daytime leaf temperature, for the Sun-height row
# Parts of the cycle by the Sun's height, degrees: the twilight after sunset and before sunrise, low, middle and high
# Sun.
PARTS = dict(twilight=(-90.0, 0.0), low_sun=(0.0, 10.0), morning_and_afternoon=(10.0, 30.0),
             middle_sun=(30.0, 60.0), high_sun=(60.0, 90.0))


class DimmedTables(plant.Tables):
    """The plant model's photosynthesis tables for the Moon with its light, direct, diffuse and twilight, scaled."""

    def __init__(self, scale):
        lam, light, self.shas = plant.light()
        src = light['moon']
        self.light = dict(moon=dict(src, direct=src['direct'] * scale, diffuse=src['diffuse'] * scale,
                                    par=src['par'] * scale))
        heights = src['elevations']
        positive = heights[heights > 0]
        optics = canopy.Optics(lam, positive)
        co2 = canopy.WORLDS['moon']['co2_pa']
        src = self.light['moon']
        self.gpp = {}
        for lai in plant.LAIS:
            stand = replace(canopy.Stand(), lai=lai)
            resp = optics(stand)
            direct_resp = [resp['direct'][int(np.argmin(np.abs(positive - max(e, positive.min()))))] for e in heights]
            self.gpp['moon', lai] = np.array([[canopy.rates(stand, (direct_resp[i], resp['diffuse']), lam,
                                                            src['direct'][i], src['diffuse'][i],
                                                            src['reflectance_from_below'], co2, t_c=t)['gpp_umol_m2_s']
                                               for i in range(heights.size)] for t in plant.TEMPS_C])
        self.elevations = heights
        self.scale = scale


def elasticity(base, dimmed, d):
    """The share lost for each share of light taken."""
    return (base - dimmed) / base / d if base else 0.0


def cycle_trace(latitude_deg):
    """Hour angle and Sun height at each step of one solar cycle, and the step in seconds."""
    hour = (np.arange(plant.STEPS) + 0.5) / plant.STEPS * 360.0 - 180.0
    elevation = np.degrees(np.arcsin(np.cos(np.radians(latitude_deg)) * np.cos(np.radians(hour))))
    return hour, elevation, plant.MOON_DAY_S / plant.STEPS


def canopy_gpp(tables, lai, latitude_deg, temp):
    """Canopy gross photosynthesis over one solar cycle (g C per m2 of ground), twilight included."""
    _, elevation, dt = cycle_trace(latitude_deg)
    gpp = tables.photosynthesis('moon', lai, elevation, temp)
    return float(np.sum(gpp) * dt * plant.GRAMS)


def by_part(tables_by_d, latitude_deg, temp, d):
    """Where in the cycle the deep canopy's photosynthesis and its loss to a dimming d fall, by the Sun's height."""
    _, elevation, dt = cycle_trace(latitude_deg)
    base = tables_by_d[0.0].photosynthesis('moon', 5.0, elevation, temp)
    dim = tables_by_d[d].photosynthesis('moon', 5.0, elevation, temp)
    total, lost = float(np.sum(base)), float(np.sum(base - dim))
    out = {}
    for name, (lo, hi) in PARTS.items():
        m = (elevation >= lo) & (elevation < hi) if hi < 90.0 else (elevation >= lo)
        b, x = float(np.sum(base[m])), float(np.sum(base[m] - dim[m]))
        out[name] = dict(hours=float(np.sum(m) * dt / 3600.0), share_of_photosynthesis=b / total if total else 0.0,
                         share_of_the_loss=x / lost if lost else 0.0, elasticity=x / b / d if b else 0.0)
    return out


EVIDENCE = ("The plant model's carbon budget (plant.py) recomputed with the Moon's light scaled, at the chosen "
            "climate's temperatures and the Earth calibration of stem and root upkeep: a clear-sky, one-dimensional "
            "stand with ample water, nutrients and sinks for its sugar, and no acclimation to the dimmer light. The "
            "response is a property of the canopy and plant models; where water, nutrients or sinks limit growth "
            "the light matters less.")
READING_RULE = ("Carbon is g C per m2 of ground over one solar cycle, as in plant.json; bands are |latitude| ranges "
                "of land on the near side. elasticity[d] is the share of photosynthesis, growth or store lost per "
                "share of light taken at dimming d (a store elasticity below zero is a larger store). by_lai gives "
                "the canopy's gross photosynthesis alone at each leaf area. by_sun_height is the deep canopy (leaf "
                "area index 5) at 22 C leaf temperature; by_part splits the deep canopy's photosynthesis over the "
                "cycle and its loss at a 5% dimming by the Sun's height, twilight below the horizon.")


def run():
    hour_centres, traces, climate_sha = plant.temperatures()
    base = plant.Tables()
    equator = plant.along_cycle(hour_centres, traces[SIDE, '0'])
    stem_root_25 = plant.calibrate(base, equator)
    tables = {0.0: base, **{d: DimmedTables(1.0 - d) for d in DIMMINGS}}
    bands = {}
    for band in plant.BANDS:
        temp = plant.along_cycle(hour_centres, traces[SIDE, band])
        out = {}
        for name in STRATEGIES:
            night, regrow = plant.STRATEGIES[name]
            runs = {d: plant.cycle(t, 'moon', float(band), temp, plant.MOON_DAY_S, stem_root_25, night=night,
                                   regrow=regrow) for d, t in tables.items()}
            keys = ('gpp', 'twilight_gpp', 'growth', 'store', 'carbon_use_efficiency')
            out[name] = dict(base={k: runs[0.0][k] for k in keys},
                             dimmed={str(d): {k: runs[d][k] for k in keys} for d in DIMMINGS},
                             elasticity={str(d): {k: elasticity(runs[0.0][k], runs[d][k], d)
                                                  for k in ('gpp', 'growth', 'store')} for d in DIMMINGS})
        by_lai = {}
        for lai in plant.LAIS:
            gpp = {d: canopy_gpp(t, lai, float(band), temp) for d, t in tables.items()}
            by_lai[f'{lai:g}'] = dict(gpp=gpp[0.0], elasticity={str(d): elasticity(gpp[0.0], gpp[d], d)
                                                                  for d in DIMMINGS})
        out['by_lai'] = by_lai
        out['by_part_at_5_percent'] = by_part(tables, float(band), temp, 0.05)
        bands[band] = out
    k = int(np.argmin(np.abs(plant.TEMPS_C - LEAF_C)))
    elev = base.elevations
    row = lambda t: t.gpp['moon', 5.0][k]
    sun = dict(elevation_deg=elev.tolist(), par_umol_m2_s=base.light['moon']['par'].tolist(),
               gpp_umol_m2_s=row(base).tolist(),
               elasticity={str(d): [elasticity(a, b, d) for a, b in zip(row(base), row(tables[d]))]
                           for d in DIMMINGS})
    reference = json.loads((RESULTS / 'plant.json').read_text())['results'][SIDE]
    for band, out in bands.items():
        for name in STRATEGIES:
            stored = reference[band][name]['growth']
            if abs(out[name]['base']['growth'] - stored) > 1e-3 * abs(stored):
                raise ValueError(f'undimmed growth at band {band}, {name}, is not plant.json\'s')
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    return dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={p: digest(p) for p in (
            'biosphere/canopy/light_response.py', 'biosphere/canopy/plant.py', 'biosphere/canopy/model.py',
            'biosphere/canopy/radiation.py', 'biosphere/canopy/leaf.py', 'biosphere/canopy/leaf_optics.py',
            'biosphere/long_night.py')},
            inputs=dict(**{f'{k}_sha256': v for k, v in base.shas.items()}, climatology_sha256=climate_sha,
                        climatology_run=plant.CLIMATE_RUN, climatology_years=plant.CLIMATE_YEARS,
                        plant_sha256=digest('biosphere/canopy/results/plant.json'))),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        settings=dict(dimmings=DIMMINGS, strategies=STRATEGIES, side=SIDE, bands_abs_latitude=plant.BANDS,
                      stem_root_respiration_25c=stem_root_25, sun_height_row_leaf_c=float(plant.TEMPS_C[k]),
                      parts_by_sun_height_deg=PARTS),
        bands=bands, by_sun_height=sun)


def main():
    product = canopy._round(run())
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'light_response.json').write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', RESULTS / 'light_response.json')


if __name__ == '__main__':
    main()
