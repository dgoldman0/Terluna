"""Run bounded commercial-rollout and component-life interventions.

The baseline model and its frozen products remain unchanged. Each experiment
supplies its own complete dated orderbook and an explicit engineering override.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import time

from . import bootstrap, model, package, run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ORDER_SCHEMA = 'terluna.humanities.enclosed-bottleneck-orders/1'
VARIANTS = {'baseline', 'demand_25', 'capture_8', 'service_rollout',
            'combined_rollout', 'demand_and_rollout', 'tug_life_30'}
HELPERS = ('model.py', 'run.py', 'bootstrap.py', 'package.py', 'bottlenecks.py')


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_identity():
    return {str((HERE/name).relative_to(ROOT)): digest(HERE/name) for name in HELPERS}


def engineering_for_case(reference, experiment):
    """Apply only the stated active-life experiment, preserving other inputs."""
    variant = experiment['variant']
    if variant not in VARIANTS:
        raise ValueError('Unknown bottleneck variant: '+str(variant))
    overrides = experiment.get('engineering_overrides', {})
    expected = {'asteroid_tug_active_life_years': 30.0} if variant == 'tug_life_30' else {}
    if overrides != expected:
        raise ValueError('Variant and engineering override disagree')
    data = deepcopy(reference)
    original = data['templates']['asteroid_tugs']['life']
    if original != 20.0:
        raise ValueError('Changed reference source-tug active life')
    if overrides:
        data['templates']['asteroid_tugs']['life'] = overrides['asteroid_tug_active_life_years']
    record = dict(original_active_life_years=original,
                  applied_active_life_years=data['templates']['asteroid_tugs']['life'],
                  structural_life_years=data['templates']['asteroid_tugs']['structural_life'],
                  hardware_USD2026_kg=data['templates']['asteroid_tugs']['price'],
                  override=deepcopy(overrides))
    return data, record


def validate_book(book):
    if book.get('schema') != ORDER_SCHEMA:
        raise ValueError('Unexpected diagnostic orderbook schema')
    seen = set()
    for case in book['cases']:
        name = case['name']; experiment = case['experiment']
        if name in seen:
            raise ValueError('Duplicate diagnostic case')
        seen.add(name)
        if name != experiment['baseline_name']+':'+experiment['variant']:
            raise ValueError('Case identity and experiment disagree')
        if case['market'] != 'breakthrough' or case['market_growth_name'] != 'reference':
            raise ValueError('Diagnostic market boundary changed')
    if not seen:
        raise ValueError('Empty diagnostic orderbook')


def cache_identity(case, engineering_sha256, helpers):
    return hashlib.sha256(canonical(dict(case=case, engineering_sha256=engineering_sha256,
                                        helpers=helpers))).hexdigest()


def solve_case(reference, case, engineering_sha256, helpers, cache_dir=None):
    data, override = engineering_for_case(reference, case['experiment'])
    identity = cache_identity(case, engineering_sha256, helpers)
    cache = None if cache_dir is None else Path(cache_dir)/(identity+'.json')
    if cache is not None and cache.exists():
        saved = json.loads(cache.read_text())
        if saved['identity'] != identity or saved['result_sha256'] != hashlib.sha256(canonical(saved['case'])).hexdigest():
            raise ValueError('Changed diagnostic cache entry')
        return saved['case'], True
    try:
        result = model.solve(data, case['annual'], max_iterations=50)
        out = run.compact_case(case, result, 'base')
    except (ValueError, OverflowError) as error:
        out = {k:case[k] for k in ('name','portfolio','market','market_growth','market_growth_name')}
        out.update(seed=None, annual=[], convergence=dict(converged=False,error=str(error)),
                   first_capacity_shortfall=case['annual'][0]['year'])
    out['physical_scenario'] = case['experiment']['variant']
    out['experiment'] = deepcopy(case['experiment'])
    out['engineering_intervention'] = override
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(canonical(dict(identity=identity, result_sha256=hashlib.sha256(canonical(out)).hexdigest(), case=out)))
    return out, False


def product(reference, book, order_path, *, cases=None, cache_dir=None):
    validate_book(book)
    order_path = Path(order_path).resolve()
    engineering_path = HERE/'inputs/engineering.json'
    identity = json.loads((HERE/'inputs/identity.json').read_text())
    engineering_hash = digest(engineering_path)
    if engineering_hash != identity['engineering_sha256']:
        raise ValueError('Changed reference engineering snapshot')
    if canonical(reference) != canonical(json.loads(engineering_path.read_text())):
        raise ValueError('Unrecorded reference engineering modification')
    if canonical(book) != canonical(json.loads(order_path.read_text())):
        raise ValueError('Orderbook bytes and supplied cases disagree')
    helpers = source_identity(); outputs=[]
    for case in (book['cases'] if cases is None else cases):
        started = time.monotonic()
        out, cached = solve_case(reference, case, engineering_hash, helpers, cache_dir)
        outputs.append(out)
        print(case['name'], 'verified-cache' if cached else 'computed',
              'iterations',out['convergence'].get('iterations'),
              'first_shortfall',out['first_capacity_shortfall'],
              'seconds',round(time.monotonic()-started,2), flush=True)
    try:
        order_key = str(order_path.relative_to(ROOT))
    except ValueError:
        order_key = str(order_path)
    result = dict(schema=run.SCHEMA,
        producer=dict(path=str(Path(__file__).relative_to(ROOT)),sha256=digest(__file__),
                      model_path='provisioning/enclosed_industry/model.py',model_sha256=digest(HERE/'model.py')),
        dependencies={**helpers, 'provisioning/enclosed_industry/inputs/engineering.json':engineering_hash,
                      'provisioning/enclosed_industry/inputs/identity.json':digest(HERE/'inputs/identity.json'),
                      order_key:digest(order_path)},
        orderbook=dict(schema=book['schema'],path=order_key,sha256=digest(order_path),
                       producer=book.get('producer'),dependencies=book.get('dependencies'),
                       experiment_spec=book.get('experiment_spec',book.get('experiments'))),
        evidence='Conditional intervention experiment on the unchanged demand-sized physical model. Commercial requests and longer component life are scenarios, not demonstrated market growth or durability improvements. Every intervention pays its resulting assets, replacement, freight and seed fuel.',
        units='kg, W, persons, constant USD2041. Annual physical flows are per year.',
        reading_rule='Case names exactly match the diagnostic orders. Only converged cases without dated capacity or workforce failures qualify for a completed financial account. Seed cost/feed precede annual year8 and are included once. Capital and operating totals already include resource provisions. The source-tug experiment changes active mean life20 to30 years, keeping structural life50, price and all other technology unchanged.',
        shared_constants={k:getattr(model,k) for k in ('JULIAN_DAY','JULIAN_YEAR_DAYS','STANDARD_GRAVITY','STEFAN_BOLTZMANN')},
        scenarios=dict(lunar=reference['industrial']['lunar'],source=reference['industrial']['source'],
                       lifetime_shape=reference.get('lifetime_shape',3.),
                       case_interventions={x['name']:x['engineering_intervention'] for x in outputs}),
        rules=[
            'The baseline physical model, engineering snapshot and original products remain unchanged. Each case is solved independently from its complete annual order schedule.',
            'Market variants change requested activity; the physical solver prices and sizes actual machinery and transport. A failed workforce or capacity screen remains a failure, not a funded delivery.',
            'Longer active source-tug life is an unqualified durability sensitivity, with no claimed improvement in production, price, propellant consumption per trip or structural lifetime.',
            'Caches require identical case content, engineering bytes and every model/helper/diagnostic-runner hash. Cached content is also checked against its recorded result hash.',
            'Seed and annual direct costs retain the baseline cost boundary; installation and contingency allowances belong to the separately identified financial consumer.'
        ], cases=outputs)
    before_hash = hashlib.sha256(canonical(result)).hexdigest()
    result = bootstrap.apply(result, reference)
    result['postprocessing'] = dict(path='provisioning/enclosed_industry/bootstrap.py',sha256=digest(HERE/'bootstrap.py'),
        source_product_sha256=before_hash,source_serialization='Canonical sorted compact JSON plus newline before bootstrap.apply',
        engineering_sha256=engineering_hash,scope='Initial external route-fuel procurement included once; annual physics unchanged.')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orders',type=Path,default=HERE/'inputs/bottleneck_orders.json')
    parser.add_argument('--case',action='append')
    parser.add_argument('--output',type=Path,default=HERE/'results/bottlenecks.json.gz')
    parser.add_argument('--no-cache',action='store_true')
    args=parser.parse_args()
    reference=json.loads((HERE/'inputs/engineering.json').read_text())
    book=json.loads(args.orders.read_text())
    chosen=[c for c in book['cases'] if not args.case or c['name'] in args.case]
    if not chosen:raise ValueError('No matching diagnostic cases')
    cache=None if args.no_cache else ROOT/'research/runs/enclosed_industry/bottlenecks/cache'
    result=product(reference,book,args.orders,cases=chosen,cache_dir=cache)
    raw=ROOT/'research/runs/enclosed_industry/bottlenecks/product.json'
    raw.parent.mkdir(parents=True,exist_ok=True);raw.write_bytes(canonical(result))
    manifest=package.package(raw,args.output)
    print(args.output);print(manifest['sha256'])


if __name__=='__main__':main()
