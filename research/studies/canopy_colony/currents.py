"""How fast the south sea's water moves at each depth, and how fast the water at a colony's feeders and at its
trunk's foot moves past a canopy that drifts with the water of its own depth.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.canopy_colony.currents   # seconds

It reads the yearly means of two runs of the south sea's circulation on the research drive (climate/ocean/runs, the
south sea study's): the last year of the 50-year run and the last finished year of the oxygen run under the Moon's
own mixing. It changes neither. A yearly mean holds the steady part of the flow; the currents that turn with the
month and the inertial period are averaged out of it, so the speeds here are low values for the water's motion at any
moment. From the same means it compares the surface water with the water at the foot of a trunk: surface water
brought to the foot's temperature sinks by itself where its salt makes it the denser, and has to be pushed down
where the rivers have freshened it. Writes results/currents.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'currents.json'
SCHEMA = 'terluna.research.canopy-colony-currents/1'
RUNS = (('the 50-year run, model year 50', 'south_sea_1deg'),
        ('the oxygen run under the Moon\'s own mixing, model year 53', 'south_sea_1deg_lunar_mixing'))
CANOPY_M = (20.0, 50.0)
FEEDERS_M = (150.0, 250.0)
FOOT_M = (250.0, 500.0)
FEET_M = ((250.0, 500.0), (250.0, 900.0))       # the trunks' feet the surface water is compared with
SURFACE_M = (0.0, 35.0)
YEAR_S = 31557600.0


def layer_mean(field, volume, mid, span):
    """A field's mean over the levels whose centres lie in span, weighed by each column's open water there; nan where
    a column has none."""
    k = (mid >= span[0]) & (mid <= span[1])
    w = volume[k].sum(axis=0)
    return np.where(w > 0, (field[k] * volume[k]).sum(axis=0) / np.maximum(w, 1e-30), np.nan)


def weighted_percentiles(values, weights, percents):
    ok = np.isfinite(values) & (weights > 0)
    v, w = values[ok], weights[ok]
    order = np.argsort(v)
    v, w = v[order], w[order]
    return [float(np.interp(p / 100.0, np.cumsum(w) / w.sum(), v)) for p in percents]


def summary(values, area):
    p10, p50, p90 = weighted_percentiles(values, area, (10, 50, 90))
    ok = np.isfinite(values)
    return dict(mean=float((values[ok] * area[ok]).sum() / area[ok].sum()), p10=p10, median=p50, p90=p90)


def read_run(name):
    from climate.ocean import mitgcm
    from research.studies.south_sea_ocean import circulation as C
    sea = C.sea_in(name)
    stems = sorted(sea.folder.glob('deep.*.data'))
    deep = mitgcm.read_mds(str(stems[-1])[:-5])
    u, v = C.at_centres(deep['UVEL'], deep['VVEL'])
    return dict(u=u, v=v, volume=C.cell_volumes(sea), mid=np.cumsum(sea.delr) - 0.5 * sea.delr, area=sea.grid.area,
                salt=deep['SALT'], theta=deep['THETA'],
                year=int(stems[-1].name.split('.')[1]) * sea.deltat / C.YEAR_S, file=stems[-1].name)


def trunk_water(run, density, foot_m):
    """Surface water brought to the temperature of the water at a trunk's foot, against that water: its excess
    density (kg/m3; negative where it is the lighter and must be pushed down), over the columns that reach the foot."""
    volume, mid, area = run['volume'], run['mid'], run['area']
    s_top, t_top = (layer_mean(run[k], volume, mid, SURFACE_M) for k in ('salt', 'theta'))
    s_foot, t_foot = (layer_mean(run[k], volume, mid, foot_m) for k in ('salt', 'theta'))
    excess = density(s_top, t_foot) - density(s_foot, t_foot)
    ok = np.isfinite(excess)
    p5, p50, p95 = weighted_percentiles(excess, np.where(ok, area, 0.0), (5, 50, 95))
    mean = lambda x: float((x[ok] * area[ok]).sum() / area[ok].sum())                      # noqa: E731
    return dict(foot_m=list(foot_m), share_of_area_denser=float(area[ok & (excess > 0)].sum() / area[ok].sum()),
                excess_density_kg_m3=dict(p5=p5, median=p50, p95=p95), temperature_step_k=mean(t_top - t_foot),
                surface_at_its_own_warmth_lighter_by_kg_m3=mean(density(s_foot, t_foot) - density(s_top, t_top)))


def reduce(run):
    """Speeds at each level and the relative flow past a canopy."""
    u, v, volume, mid, area = run['u'], run['v'], run['volume'], run['mid'], run['area']
    by_level = []
    for k in range(len(mid)):
        open_ = volume[k] > 0
        if open_.sum() < 20:
            break
        speed = np.where(open_, np.hypot(u[k], v[k]), np.nan)
        by_level.append(dict(depth_m=float(mid[k]), **summary(speed, np.where(open_, area, 0.0))))
    cu, cv = layer_mean(u, volume, mid, CANOPY_M), layer_mean(v, volume, mid, CANOPY_M)
    past = {}
    for label, span in (('feeders', FEEDERS_M), ('foot', FOOT_M)):
        du, dv = layer_mean(u, volume, mid, span) - cu, layer_mean(v, volume, mid, span) - cv
        past[label] = dict(depths_m=list(span), **summary(np.hypot(du, dv), area))
    drift = summary(np.hypot(cu, cv), area)
    return dict(model_year=run['year'], file=run['file'], speed_by_depth_m_s=by_level, canopy_depths_m=list(CANOPY_M),
                canopy_drift_m_s=drift, canopy_drift_km_a_year=drift['median'] * YEAR_S / 1e3,
                water_past_a_drifting_canopy_m_s=past)


def build():
    from research.studies.south_sea_ocean import run as R
    runs = {label: read_run(name) for label, name in RUNS}
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.canopy_colony.currents',
                      code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        evidence=('Read from yearly means of MITgcm runs of the South Pole-Aitken Sea at lunar gravity on 30 km cells '
                  '(the south sea study): the steady part of a modelled flow on the atlas\'s basin, with closure '
                  'constants fitted on Earth. Not validated against any observation.'),
        reading_rule=('runs: for each run, the speed of the yearly-mean flow at each level over the sea\'s open cells '
                      '(area-weighted mean, 10th percentile, median and 90th), the drift of the water at the canopy\'s '
                      'depths, and the speed at which the water at the feeders\' depths and at the trunk\'s foot moves '
                      'relative to the water at the canopy\'s depths in the same column. Monthly and inertial '
                      'currents are averaged out, so these are low values. surface_water_cooled_to_the_foot: in the '
                      'oxygen run\'s year, the density of the top 35 m\'s water at the temperature of the water '
                      'at a trunk\'s foot, less that water\'s (kg/m3; the share of the sea\'s area where it is the '
                      'denser, and the 5th percentile, median and 95th by area).'),
        runs={label: reduce(run) for label, run in runs.items()},
        surface_water_cooled_to_the_foot={f'{foot[0]:g}-{foot[1]:g} m': trunk_water(runs[RUNS[1][0]], R.seawater_density, foot)
                                          for foot in FEET_M})


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    for label, r in out['runs'].items():
        past = r['water_past_a_drifting_canopy_m_s']
        print(f'{label}: canopy drift {r["canopy_drift_m_s"]["median"] * 100:.1f} cm/s; past the feeders '
              f'{past["feeders"]["median"] * 100:.1f} cm/s (10th {past["feeders"]["p10"] * 100:.1f}); past the foot '
              f'{past["foot"]["median"] * 100:.1f} cm/s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
