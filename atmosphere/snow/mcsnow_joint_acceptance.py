"""Pointwise acceptance against the frozen column-density diagnostic.

The contract includes settling, alternate starts and supplied-state rates. Its
10% requirement is a requested calibration criterion, not observational error.
Numerical conservation checks retain their separate, much tighter tolerances.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASELINE = HERE / 'results/mcsnow_column_density.json'
BASELINE_SHA256 = '304b90e52f4c0269cc0959d8157c91c9d0c2aee2def1d8f950bb0a9bce4db251'
BASELINE_COMMIT = 'ded734e'
PRODUCT = HERE / 'results/mcsnow_joint_acceptance.json'
TOLERANCE = .1
# Absolute relative-error slack for floating-point roundoff only.
ROUNDING = 32 * math.ulp(1.)
QUANTITIES = ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m', 'fall_speed_m_s')
RATES = ('mass_rate_kg_s', 'equatorial_log_rate_s_1', 'polar_log_rate_s_1')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_baseline():
    if sha(BASELINE) != BASELINE_SHA256:
        raise ValueError('Frozen non-regression baseline changed')
    return json.loads(BASELINE.read_text())


def point_id(row):
    return json.dumps([row[k] for k in ('family', 'temperature_C', 'start', 'age_min', 'quantity')],
                      separators=(',', ':'), allow_nan=False)


def inventory(warm, cold):
    """Extract actual values and targets; never trust a stored ratio or score."""
    points, invariants = {}, []

    def add(family, temperature, start, age, quantity, value, target, origin):
        if not all(map(math.isfinite, (temperature, age, value, target))) or min(value, target) <= 0:
            raise ValueError('Finite positive physical values and targets required')
        row = dict(family=family, temperature_C=temperature, start=start, age_min=age,
                   quantity=quantity, value=value, target=target, origins=[origin])
        identity = point_id(row)
        if identity in points:
            old = points[identity]
            if not all(math.isclose(old[k], row[k], rel_tol=ROUNDING, abs_tol=0) for k in ('value', 'target')):
                raise ValueError('Inconsistent duplicate point: '+identity)
            old['origins'].append(origin)
        else:
            points[identity] = row

    def assess(a, family, temperature, origin, seed=None):
        if a.get('failed_cases') or a.get('error') or not a.get('rows'):
            raise ValueError('Failed or empty assessment: '+origin)
        for r in a['rows']:
            if r['temperature_C'] != temperature or len(r['values']) != 4 or len(r['target']) != 4:
                raise ValueError('Unexpected assessment shape or temperature')
            start = r['initialization'] if seed is None else 'seed_phi='+format(seed, '.17g')
            for q, v, t in zip(QUANTITIES, r['values'], r['target']):
                add(family, temperature, start, r['age_min'], q, v, t, origin)
        invariants.append(dict(origin=origin, mass_balance_max=a['mass_balance_max'],
            volume_residual_max=a['volume_residual_max'], maximum_density_fraction=a['maximum_density_fraction']))

    for temperature, fit in ((-5.3, warm), (-14.4, cold)):
        for label in ('assessment', 'additional_age_assessment', 'interval_assessment'):
            if label in fit:
                assess(fit[label], 'trajectory', temperature, f'{temperature}/{label}')
        for probe in fit['fixed_state_probes']:
            origin=f'{temperature}/fixed_state_probes'
            for q, v, t in zip(RATES, probe['rates'], probe['target_rates'], strict=True):
                add('supplied_state', temperature, probe['mapping'], probe['age_min'], q, v, t, origin)
            for q in ('effective_axis_ratio', 'deposition_density_kg_m3'):
                add('supplied_state', temperature, probe['mapping'], probe['age_min'], q,
                    probe[q], probe['target_'+q], origin)
        for seed in fit['seed_shape_sensitivity']:
            phi=seed['initial_aspect_ratio']
            assess(seed['assessment'], 'seed_sensitivity', temperature,
                   f'{temperature}/seed_phi={phi}', seed=phi)
    return dict(points=[dict(id=k, **points[k]) for k in sorted(points)], invariants=invariants)


def baseline_inventory():
    p=load_baseline()
    return inventory(next(f for f in p['searches'] if f['name']==p['preferred_search']), p['previous_cold'])


def compare(candidate, baseline=None):
    """Fail closed for missing/extra/duplicate points, altered targets or invariants."""
    if baseline is None:
        baseline=baseline_inventory()
    expected={p['id']:p for p in baseline['points']}
    found={p['id']:p for p in candidate['points']}
    if len(found)!=len(candidate['points']) or set(found)!=set(expected):
        raise ValueError('Candidate must cover every baseline point exactly once')
    rows=[]
    for identity, b in expected.items():
        c=found[identity]
        if point_id(c)!=identity or c['target']!=b['target']:
            raise ValueError('Candidate identity or target changed: '+identity)
        value=c['value']
        if not math.isfinite(value) or value<=0:
            raise ValueError('Invalid candidate value: '+identity)
        before=abs(b['value']/b['target']-1)
        after=abs(value/b['target']-1)
        rows.append(dict(**c, baseline_value=b['value'], baseline_absolute_relative_error=before,
            absolute_relative_error=after, allowed_error=min(TOLERANCE, before),
            within_ten_percent=after<=TOLERANCE+ROUNDING,
            non_regressing=after<=before+ROUNDING))
    expected_origins={r['origin'] for r in baseline['invariants']}
    supplied=candidate['invariants']
    if len(supplied)!=len(expected_origins) or {r['origin'] for r in supplied}!=expected_origins:
        raise ValueError('Candidate lacks complete numerical checks')
    checks=[]
    for r in supplied:
        vals=[r[k] for k in ('mass_balance_max','volume_residual_max','maximum_density_fraction')]
        passed=all(math.isfinite(x) for x in vals) and 0<=vals[0]<=1e-10 and 0<=vals[1]<=1e-10 and 0<vals[2]<=1+1e-10
        checks.append(dict(**r, passed=passed))
    groups=defaultdict(list)
    for r in rows:
        groups[(r['family'],r['temperature_C'],r['quantity'])].append(r)
    summaries=[]
    for (family,temp,q), group in sorted(groups.items()):
        worst=max(group,key=lambda r:r['absolute_relative_error'])
        summaries.append(dict(family=family,temperature_C=temp,quantity=q,points=len(group),
            maximum_absolute_relative_error=worst['absolute_relative_error'],worst_point=worst['id'],
            above_ten_percent=sum(not r['within_ten_percent'] for r in group),
            regressions=sum(not r['non_regressing'] for r in group)))
    return dict(accepted=all(r['within_ten_percent'] and r['non_regressing'] for r in rows) and all(r['passed'] for r in checks),
        point_count=len(rows), above_ten_percent=sum(not r['within_ten_percent'] for r in rows),
        regression_count=sum(not r['non_regressing'] for r in rows),
        summaries=summaries, rows=rows, numerical_checks=checks)


def structural_limits():
    """Necessary constraints that the preceding constant/monotone laws fail.

    Reference densities are derived from fitted histories and envelope mappings.
    These bounds concern matching those diagnostic targets, not empirical truth.
    """
    p=load_baseline();warm=next(f for f in p['searches'] if f['name']==p['preferred_search'])
    cold=p['previous_cold']['fixed_state_probes']
    targets=[r['target_deposition_density_kg_m3'] for r in cold]
    low,high=min(targets),max(targets)
    # Equal relative errors at the extreme targets determine the minimax constant.
    constant=2*low*high/(low+high)
    selected=[r for r in warm['fixed_state_probes'] if r['mapping']=='axes' and r['age_min'] in (3.5,5.5,8.)]
    selected.sort(key=lambda r:r['age_min'])
    intervals=[dict(age_min=r['age_min'],target=r['target_deposition_density_kg_m3'],
        lower=.9*r['target_deposition_density_kg_m3'],upper=1.1*r['target_deposition_density_kg_m3']) for r in selected]
    return dict(
        cold_constant_density=dict(target_minimum_kg_m3=low,target_maximum_kg_m3=high,
            best_constant_kg_m3=constant,minimum_possible_maximum_relative_error=(high-low)/(high+low),
            ten_percent_intersection_lower=.9*high,ten_percent_intersection_upper=1.1*low,
            feasible=.9*high<=1.1*low,
            scope='A single constant porous-deposition density at -14.4 C cannot fit all supplied-state density targets, even before non-regression constraints.'),
        warm_monotone_density=dict(intervals=intervals,
            nondecreasing_possible=all(a['lower']<=b['upper'] for a,b in zip(intervals,intervals[1:])),
            nonincreasing_possible=all(a['upper']>=b['lower'] for a,b in zip(intervals,intervals[1:])),
            scope='The three fitted axes-mapping states have increasing equivalent radius and maximum dimensions above the allowed 80 micrometre compact-stage limit. The compact solid-deposition branch cannot resolve the required decline followed by rise. This excludes a single monotone early-to-late deposition-density transition at these states.'))


def assemble():
    baseline=baseline_inventory()
    return dict(schema='terluna.atmosphere.mcsnow-joint-acceptance/1',
        baseline_path=str(BASELINE.relative_to(HERE.parent.parent)),baseline_sha256=BASELINE_SHA256,
        baseline_commit=BASELINE_COMMIT,relative_error_limit=TOLERANCE,roundoff_slack=ROUNDING,
        acceptance='Every point must have absolute relative error <= min(0.1, its frozen baseline error); complete coverage and numerical invariants are also required.',
        scope='Two current calibrated temperatures (-5.3 and -14.4 C), every recorded trajectory age, both principal starts, alternate volume mapping, five seed shapes and supplied-state growth diagnostics. No claim of other-temperature or lunar validation.',
        evidence=dict(trajectory='Comparison with published population-history fits, including prescribed restart values; fits were used for calibration.',
            seed_sensitivity='Robustness criterion against the same population fit; these seeds were not individually measured.',
            supplied_state='Derivatives and density inferred from piecewise population fits and the declared geometry, not directly observed single-crystal rates or porosity.'),
        baseline_assessment=compare(baseline,baseline),structural_limits=structural_limits(),producer_sha256=sha(__file__))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('assemble',))
    parser.parse_args()
    PRODUCT.write_text(json.dumps(assemble(),indent=2,allow_nan=False)+'\n')
