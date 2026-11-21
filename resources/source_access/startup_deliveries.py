"""Conditional material received on Earth by the end of the eight-year founding build."""
import gzip
import json
import math
from pathlib import Path

from . import coproducts as c

HERE = Path(__file__).resolve().parent
SOURCE_METRICS = ('PGM6', 'Pt', 'Pd', 'Rh', 'Ru', 'Ir', 'Os', 'Au', 'Ag', 'Ni', 'Co', 'Ga', 'Ge', 'In', 'REY', 'NdPr', 'Sc')
LUNAR_METRICS = ('REY', 'NdPr', 'Sc', 'Dy', 'Tb')


def allocate_feed(record, weights):
    if len(weights) != 8 or any(not math.isfinite(v) or v < 0 for v in weights) or sum(weights) <= 0:
        raise ValueError('Eight nonnegative weights with a positive total required')
    if not math.isclose(record['base_kg']+record['late_years_6_7_kg'], record['total_kg'], rel_tol=1e-12):
        raise ValueError('Seed inventory does not close')
    annual = [record['base_kg']*v/sum(weights) for v in weights]
    for i in (6, 7):
        annual[i] += record['late_years_6_7_kg']/2
    return annual


def deliveries(feed, shares, grade_ppm, recovery, export_share, delay, loss):
    if len(feed) != 8 or len(shares) != 8 or not math.isfinite(delay) or delay < 0:
        raise ValueError('Eight production intervals and a nonnegative delivery lag required')
    c.bounded(export_share, 0, 1)
    c.bounded(loss, 0, 1)
    refined = [c.element_balance(f, grade_ppm, s, recovery)['refined_kg'] for f, s in zip(feed, shares)]
    attempted_arrivals = [0.]*8
    for origin, amount in enumerate(refined):
        # Uniform production within each assumed annual interval, translated by
        # the delivery lag. Fractional years do not round an entire cargo away.
        for arrival in range(8):
            overlap = max(0., min(origin+1+delay, arrival+1)-max(origin+delay, arrival))
            attempted_arrivals[arrival] += amount*export_share*overlap
    attempted = math.fsum(attempted_arrivals)
    total = math.fsum(refined)
    return {'refined_total_kg': total,
            'retained_off_Earth_kg': total*(1-export_share),
            'awaiting_or_in_transit_at_year_8_kg': max(0., total*export_share-attempted),
            'lost_on_Earth_return_by_year_8_kg': attempted*loss,
            'received_on_Earth_by_year_8_kg': attempted*(1-loss),
            'annual_Earth_receipts_kg': [v*(1-loss) for v in attempted_arrivals]}


def calculate(case, config, base, assays, final_share, weights=None, source_delay=None, export_share=None):
    weights = config['feed_weights'] if weights is None else weights
    export_share = config['Earth_allocation_fraction'] if export_share is None else export_share
    result = {'source': {}, 'lunar': {}}
    for origin, field, profile, metrics, share, lag in (
        ('source', 'source_feed', 'ci', SOURCE_METRICS, final_share,
         config['source_delivery_lag_years'] if source_delay is None else source_delay),
        ('lunar', 'lunar_raw', 'kreep15386', LUNAR_METRICS,
         config['lunar_selected_final_treated_share'], config['lunar_delivery_lag_years']),
    ):
        feed = allocate_feed(case[field], weights)
        treated = [share*v for v in config['refinery_ramp']]
        result[origin+'_feed_kg_year'] = feed
        result[origin+'_treated_feed_kg_year'] = [f*s for f, s in zip(feed, treated)]
        for metric in metrics:
            ppm = c.grade(assays['profiles'][profile], metric, assays['group_definitions'])
            result[origin][metric] = {name: deliveries(feed, treated, ppm*factor,
                base['recovery_scenarios'][name], export_share, lag, config['return_loss_fraction'])
                for name, factor in zip(c.SCENARIOS, base['grade_transfer_factors'][profile])}
    return result


def load():
    config = c.read(HERE/'startup_delivery_inputs.json')
    if config['schema'] != 'terluna.resources.startup-delivery-inputs/1':
        raise ValueError('Unexpected startup input schema')
    for name, digest in config['dependencies'].items():
        if c.sha(HERE/name) != digest:
            raise ValueError('Changed startup dependency: '+name)
    feeds = c.read(HERE/'inputs/startup_feed.json')
    if feeds['schema'] != 'terluna.resources.startup-feed/1' or len(feeds['cases']) != 81:
        raise ValueError('Expected the existing 81-case founding product')
    if feeds['producer']['sha256'] != c.sha(HERE/'export_startup_feed.py'):
        raise ValueError('Changed seed exporter')
    if any(row['opening_population'] != 12500 for row in feeds['cases']):
        raise ValueError('Opening population changed')
    if len(config['refinery_ramp']) != 8:
        raise ValueError('Eight refinery fractions required')
    for value in config['refinery_ramp']:
        c.bounded(value, 0, 1)
    base, assays, _ = c.load()
    return config, feeds, base, assays


def run():
    config, feeds, base, assays = load()
    refcase = next(r for r in feeds['cases'] if r['name'] == config['reference_case'])
    details = {}
    summaries = []
    for final in config['source_final_treated_shares']:
        for case in feeds['cases']:
            r = calculate(case, config, base, assays, final)
            summaries.append({'name': case['name'], 'source_final_treated_share': final,
                'PGM6_received_kg_low_middle_high': [r['source']['PGM6'][s]['received_on_Earth_by_year_8_kg'] for s in c.SCENARIOS],
                'REY_received_kg_low_middle_high': [sum(r[o]['REY'][s]['received_on_Earth_by_year_8_kg'] for o in ('source', 'lunar')) for s in c.SCENARIOS]})
            if case['name'] == config['reference_case']:
                details[str(final)] = r
    sensitivity = []
    for shape, weights in config['feed_shape_sensitivities'].items():
        for lag in config['source_delivery_lag_sensitivities_years']:
            r = calculate(refcase, config, base, assays, .1, weights=weights, source_delay=lag)
            sensitivity.append({'feed_shape': shape, 'source_delivery_lag_years': lag,
                'PGM6_received_kg_middle': r['source']['PGM6']['middle']['received_on_Earth_by_year_8_kg']})
    return {
        'schema': 'terluna.resources.startup-deliveries/1',
        'producer': {'path': 'resources/source_access/startup_deliveries.py', 'sha256': c.sha(__file__)},
        'dependencies': {'startup_delivery_inputs.json': c.sha(HERE/'startup_delivery_inputs.json'), **config['dependencies']},
        'evidence': 'Conditional first-eight-year physical receipts from the existing 12,500-resident founding feed. Assumed feed timing, new refining capacity, Earth allocation and delivery delays; not a qualified deployment or contracted sales forecast.',
        'units': 'kg of elemental product; annual receipts cover elapsed intervals [0,1) through [7,8).',
        'reading_rule': 'Cases and refinery scales are alternatives. PGM6 includes its six listed elements; NdPr is included in REY. Only source and lunar inventories of the same element may be added. Refining losses already occur in the recovery chain; the return loss is applied once after shipment.',
        'reference_seed': refcase, 'reference_by_final_source_share': details,
        'case_summaries': summaries, 'timing_sensitivities': sensitivity,
        'no_refinery_receipts_kg': 0., 'incremental_revenue_USD2041': None,
        'qualification': 'The source provides an eight-year total, not a buildable annual sequence. The assumed ramp peaks at 1.78 times average base feed; extra peak capacity and coproduct-plant resources are not recoupled into the feed fixed point. Earth allocation is illustrative, not estimated domestic consumption or guaranteed market absorption.',
    }


if __name__ == '__main__':
    destination = HERE/'results/startup_deliveries.json.gz'
    destination.write_bytes(gzip.compress(json.dumps(run(), sort_keys=True, separators=(',', ':')).encode(), mtime=0))
    print(destination)
