"""Pumped-vacuum transition and incremental factory continuation screens."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from numbers import Real
from pathlib import Path

from shared.provenance import ROOT
from .industrial_physics import ideal_pump_work

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'industrial_continuation_inputs.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH, root=ROOT):
    s = json.loads(Path(path).read_text())
    if s['schema'] != 'terluna.research.industrial-continuation-inputs/1':
        raise ValueError('Unexpected continuation input schema')
    products = {}
    for name, spec in s['products'].items():
        p = Path(root) / name
        if digest(p) != spec['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        d = json.loads(p.read_text())
        if d['schema'] != spec['schema']:
            raise ValueError(f'Unexpected product schema: {name}')
        products[name] = d
    return s, products


def pumped_chamber(ambient, target, volume, speed, internal_gas_load=0,
                   leak_conductance=0):
    """Constant-speed pump and illustrative molecular leak, all in SI units.

    V dp/dt = Q_internal + C*P_ambient - (S+C)*p. C is held fixed
    only for this diagnostic; real leak conductance depends on flow regime.
    The pump's own ultimate-pressure/backflow contribution is neglected.
    """
    values = (ambient, target, volume, speed, internal_gas_load, leak_conductance)
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Finite nonnegative values required')
    if min(target, volume, speed) <= 0 or ambient < target:
        raise ValueError('Positive target, volume and speed; ambient >= target')
    floor = (internal_gas_load + leak_conductance*ambient)/(speed+leak_conductance)
    reachable = floor < target
    time = (volume/(speed+leak_conductance)*math.log((ambient-floor)/(target-floor))
            if reachable and ambient > target else (0.0 if reachable else None))
    leak_at_target = leak_conductance*(ambient-target)
    return dict(ambient_pressure_Pa=ambient, target_pressure_Pa=target,
        pressure_load_N_m2=ambient-target,
        ideal_reversible_work_J=ideal_pump_work(volume, ambient, target)['ideal_reversible_work_J'],
        no_gas_load_pumpdown_s=volume/speed*math.log(ambient/target),
        fixed_internal_gas_load_Pa_m3_s=internal_gas_load,
        external_leak_at_target_Pa_m3_s=leak_at_target,
        maintenance_pumping_speed_threshold_m3_s=(internal_gas_load+leak_at_target)/target,
        chosen_pump_asymptotic_pressure_Pa=floor,
        target_reachable_with_chosen_pump=reachable,
        loaded_pumpdown_s=time)


def _number(name, value, *, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f'{name} must be a finite number')
    value = float(value)
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f'{name} must be finite' + (' and nonnegative' if nonnegative else ''))
    return value


def _year(name, value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return value


def evaluate_continuation(close_case, long_case, *, real_discount_rate,
                          future_output_kg_year, retrofit_capex_U=0,
                          added_annual_cost_U=0):
    """Return incremental and whole-project NPVs plus signed cost allowances.

    Required source keys are years and unlevered_project_npv_U. If both cases
    provide plant, cost_case, allocation_policy or initial_capital_U, they must
    match. The caller must verify the shared operating history and discount rate:
    summary NPVs alone cannot establish either. Future output is an explicit
    constant scenario, not inferred from a mean over different production years.

    A positive incremental NPV favours continuation for the already-built plant.
    A nonnegative whole-project NPV additionally recovers its original investment
    at the specified rate. Cost limits remain negative when that criterion fails
    before any additional cost; a negative limit is a required benefit, not a
    permissible expenditure. Per-kg limits are None when future output is zero.
    """
    close_year = _year('close_case.years', close_case['years'])
    final_year = _year('long_case.years', long_case['years'])
    if final_year <= close_year:
        raise ValueError('Continuation must finish after the closure decision')
    for key in ('plant', 'cost_case', 'allocation_policy'):
        if key in close_case and key in long_case and close_case[key] != long_case[key]:
            raise ValueError(f'Case identity differs: {key}')
    if 'initial_capital_U' in close_case and 'initial_capital_U' in long_case:
        capital_close = _number('close initial capital', close_case['initial_capital_U'], nonnegative=True)
        capital_long = _number('long initial capital', long_case['initial_capital_U'], nonnegative=True)
        if not math.isclose(capital_close, capital_long, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('Cases must share the same initial capital')
    close_npv = _number('close NPV', close_case['unlevered_project_npv_U'])
    long_npv = _number('long NPV', long_case['unlevered_project_npv_U'])
    rate = _number('real_discount_rate', real_discount_rate, nonnegative=True)
    output = _number('future_output_kg_year', future_output_kg_year, nonnegative=True)
    retrofit = _number('retrofit_capex_U', retrofit_capex_U, nonnegative=True)
    annual_cost = _number('added_annual_cost_U', added_annual_cost_U, nonnegative=True)
    remaining_years = final_year - close_year
    discount_to_initial = math.exp(-close_year * math.log1p(rate))
    if discount_to_initial == 0:
        raise ValueError('Decision-year discount factor exceeds finite numeric range')
    annuity_factor = (-math.expm1(-remaining_years*math.log1p(rate))/rate
                      if rate else float(remaining_years))
    base_increment_initial = long_npv - close_npv
    base_increment_decision = base_increment_initial / discount_to_initial
    whole_project_budget_decision = long_npv / discount_to_initial
    added_cost_decision = retrofit + annual_cost * annuity_factor
    added_cost_initial = added_cost_decision * discount_to_initial
    whole_project_npv = long_npv - added_cost_initial
    incremental_npv = whole_project_npv - close_npv
    continue_capex_limit = base_increment_decision - annual_cost*annuity_factor
    whole_capex_limit = whole_project_budget_decision - annual_cost*annuity_factor
    continue_annual_limit = (base_increment_decision-retrofit)/annuity_factor
    whole_annual_limit = (whole_project_budget_decision-retrofit)/annuity_factor
    results = dict(close_source_npv_U=close_npv,
        continuation_source_npv_U=long_npv,
        base_incremental_npv_U=base_increment_initial,
        base_incremental_value_at_decision_U=base_increment_decision,
        whole_project_added_cost_budget_at_decision_U=whole_project_budget_decision,
        added_cost_pv_at_decision_U=added_cost_decision,
        added_cost_pv_at_initial_year_U=added_cost_initial,
        continuation_project_npv_after_added_cost_U=whole_project_npv,
        incremental_npv_vs_close_U=incremental_npv,
        incremental_value_at_decision_U=incremental_npv/discount_to_initial,
        continuation_retrofit_limit_U=continue_capex_limit,
        whole_project_retrofit_limit_U=whole_capex_limit,
        continuation_annual_cost_limit_U=continue_annual_limit,
        whole_project_annual_cost_limit_U=whole_annual_limit,
        continuation_annual_cost_limit_U_per_kg=continue_annual_limit/output if output else None,
        whole_project_annual_cost_limit_U_per_kg=whole_annual_limit/output if output else None,
        added_annual_cost_U_per_kg=annual_cost/output if output else None,
        continuation_preferred_at_given_costs=incremental_npv > 0,
        whole_project_recovers_at_given_costs=whole_project_npv >= 0)
    if any(isinstance(v, float) and not math.isfinite(v) for v in results.values()):
        raise ValueError('Continuation calculation exceeds finite numeric range')
    return dict(settings=dict(decision_year=close_year, final_year=final_year,
        additional_operating_years=remaining_years, real_discount_rate=rate,
        future_output_kg_year=output, retrofit_capex_at_decision_U=retrofit,
        added_end_year_cost_U=annual_cost,
        annuity_factor_at_decision=annuity_factor,
        decision_to_initial_discount_factor=discount_to_initial),
        results=results,
        interpretation='Incremental value for an already-built plant and recovery of the original investment are separate tests. No financing, competitive repricing, retained orders or pressure-derived cost is established.')


def results():
    s, products = load_inputs()
    physics = products['research/studies/world_economy/results/industrial_physics.json']
    markets = products['research/studies/world_economy/results/industrial_markets.json']
    v = s['vacuum']
    earth = physics['enclosure_cases'][0]['against_vacuum']['internal_pressure_Pa']
    target = v['target_pressure_Pa']
    conductance = v['external_leak_load_Pa_m3_s_at_earth_pressure']/(earth-target)
    loading = next(x for x in physics['loading'] if x['start_year'] == s['loading_example']['source_start_year'])
    final = loading['reference_final_pressure_Pa']
    rows = [pumped_chamber(p, target, v['chamber_volume_m3'],
                          v['effective_pumping_speed_m3_s'],
                          v['fixed_internal_gas_load_Pa_m3_s'], conductance)
            for p in sorted(set(v['ambient_pressures_Pa']+[earth, final]))]
    ref = next(x for x in rows if x['ambient_pressure_Pa'] == earth)
    for row in rows:
        row['initial_gas_relative_to_earth_at_equal_volume_temperature'] = row['ambient_pressure_Pa']/earth
        for k in ('pressure_load_N_m2', 'ideal_reversible_work_J', 'no_gas_load_pumpdown_s'):
            row[k+'_relative_to_earth'] = row[k]/ref[k]
    pressure_years = [dict(year=y, reference_column_pressure_Pa=final*min(1,max(0,(y-loading['start_year'])/loading['loading_years'])))
                     for y in s['loading_example']['sample_years']]
    c = s['continuation']
    close = markets['reference']
    long = next(x for x in markets['scenarios'] if x['plant'] == close['plant']
                and x['cost_case'] == close['cost_case'] and x['years'] == c['continued_year'])
    assert close['years'] == c['close_year']
    extra_mass = long['totals']['produced_kg']-close['totals']['produced_kg']
    assert math.isclose(extra_mass, c['constant_annual_output_kg']*(c['continued_year']-c['close_year']), rel_tol=1e-12)
    finance_rows = [evaluate_continuation(close, long,
        real_discount_rate=c['real_discount_rate'],
        future_output_kg_year=c['constant_annual_output_kg'],
        retrofit_capex_U=k, added_annual_cost_U=annual)
        for k in c['retrofit_capex_U'] for annual in c['added_annual_cost_U']]
    names = ('industrial_continuation.py', 'industrial_continuation_inputs.json', 'industrial_physics.py')
    return dict(schema='terluna.research.industrial-continuation/1',
        producer=dict(files={str((HERE/n).relative_to(ROOT)):digest(HERE/n) for n in names}, input_products=s['products']),
        evidence=s['evidence'], vacuum=rows, vacuum_assumptions=v,
        loading_example=dict(source_start_year=loading['start_year'], end_year=loading['end_year'], samples=pressure_years,
            rule='Historical-column linear retained-gas inventory scaling only; no local-pressure or evolving 301 K simulation.'),
        continuation_assumptions=c, continuation=finance_rows, boundaries=s['boundaries'])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output', type=Path, default=ROOT/'research/runs/world_economy/industrial_continuation.json')
    a=p.parse_args(); a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(a.output)


if __name__ == '__main__':
    main()
