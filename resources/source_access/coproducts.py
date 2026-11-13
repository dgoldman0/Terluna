"""Conditional trace-element coproduct yields from existing fresh-feed scenarios.

Assay uncertainty, analogue transfer, plant recovery and selected feed share are
separate evidence states. Output envelopes are scenarios, not probability bounds.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

from . import protection

HERE = Path(__file__).resolve().parent
SCENARIOS = ('low', 'middle', 'high')
CI_METRICS = ('PGM6', 'Au', 'Ag', 'Ni', 'Co', 'REY', 'Sc', 'Ga', 'Ge', 'In')
LUNAR_METRICS = ('REY', 'NdPr', 'Sc', 'Zr', 'Nb', 'Th', 'U')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def bounded(value, lower, upper):
    if not math.isfinite(value) or not lower <= value <= upper:
        raise ValueError('Value outside finite physical bounds')
    return value


def element_balance(feed_kg, grade_ppm, treated_share, recovery):
    """One element: contained inventory partitioned exactly once along a process."""
    bounded(feed_kg, 0, math.inf)
    bounded(grade_ppm, 0, 1e6)
    bounded(treated_share, 0, 1)
    if set(recovery) != {'retention', 'separation', 'finishing'}:
        raise ValueError('Three named recovery stages are required')
    r, s, f = [bounded(recovery[k], 0, 1) for k in ('retention', 'separation', 'finishing')]
    contained = feed_kg * grade_ppm * 1e-6
    selected = contained * treated_share
    retained = selected * r
    separated = retained * s
    refined = separated * f
    return {'contained_kg': contained, 'untreated_kg': contained-selected,
            'not_retained_in_selected_stream_kg': selected-retained,
            'separation_residue_kg': retained-separated,
            'finishing_residue_kg': separated-refined, 'refined_kg': refined}


def grade(profile, metric, groups):
    """Missing assays remain unknown; a missing member invalidates a group sum."""
    elements = profile['elements']
    names = groups.get(metric, [metric])
    if any(e not in elements for e in names):
        return None
    return math.fsum(elements[e]['ppm'] for e in names)


def envelope(feed_kg, profile, metric, share, config, groups, factors):
    concentration = grade(profile, metric, groups)
    if concentration is None:
        return None
    return [element_balance(feed_kg, concentration*factor, share,
                            config['recovery_scenarios'][scenario])['refined_kg']
            for scenario, factor in zip(SCENARIOS, factors)]


def break_even(grade_ppm, recovery, marginal_cost):
    bounded(marginal_cost, 0, math.inf)
    yield_kg = element_balance(1, grade_ppm, 1, recovery)['refined_kg']
    return marginal_cost/yield_kg if yield_kg else None


def load():
    config = read(HERE/'coproduct_inputs.json')
    if config['schema'] != 'terluna.resources.coproduct-inputs/1':
        raise ValueError('Unexpected coproduct input schema')
    for path, expected in config['dependencies'].items():
        if sha(HERE/path) != expected:
            raise ValueError('Changed coproduct dependency: '+path)
    assays = read(HERE/'coproduct_assays.json')
    flows = read(HERE/'inputs/coproduct_flows.json')
    if assays['schema'] != 'terluna.resources.coproduct-assays/1' or flows['schema'] != 'terluna.resources.coproduct-flows/1':
        raise ValueError('Unexpected assay or flow schema')
    if flows['producer']['sha256'] != sha(HERE/'export_coproduct_flows.py'):
        raise ValueError('Changed flow exporter')
    if len(flows['cases']) != 81 or len({c['name'] for c in flows['cases']}) != 81:
        raise ValueError('Expected 81 distinct physical cases')
    for case in flows['cases']:
        if [r['year'] for r in case['checkpoints']] != config['checkpoints']:
            raise ValueError('Physical checkpoints changed')
    for profile in assays['profiles'].values():
        if math.fsum(e['ppm'] for e in profile['elements'].values()) > 1e6:
            raise ValueError('Transcribed elemental mass exceeds feed')
    for factors in config['grade_transfer_factors'].values():
        if len(factors) != 3 or sorted(factors) != factors:
            raise ValueError('Three ordered analogue-transfer factors required')
    return config, assays, flows


def run():
    config, assays, flows = load()
    profiles, groups = assays['profiles'], assays['group_definitions']
    per_mt = {}
    for name, profile in profiles.items():
        factors = config['grade_transfer_factors'].get(name, [1, 1, 1])
        metrics = {e: e for e in profile['elements']}
        metrics.update({g: g for g in groups})
        rows = {}
        for metric in metrics:
            ppm = grade(profile, metric, groups)
            if ppm is not None:
                refined = envelope(1e9, profile, metric, 1, config, groups, factors)
                rows[metric] = {'reference_contained_tonnes_per_Mt': ppm,
                                'refined_tonnes_per_Mt_low_middle_high': [x/1000 for x in refined]}
        per_mt[name] = rows
    # Only fresh-feed fields are multiplied. Tails and generic/recycled metals
    # remain separate inherited diagnostics, never extra raw inputs.
    annual = []
    reference = []
    for case in flows['cases']:
        for row in case['checkpoints']:
            output = {'name': case['name'], 'portfolio': case['portfolio'],
                      'admitted_completed_path': case['admitted_completed_path'], **row}
            for name, field, share, metrics in (
                ('ci', 'source_feed_kg_year', .01, CI_METRICS),
                ('kreep15386', 'lunar_raw_kg_year', .1, LUNAR_METRICS),
            ):
                output[name] = {'selected_treated_share': share,
                    'refined_kg_year_low_middle_high': {
                        metric: envelope(row[field], profiles[name], metric, share, config, groups,
                                         config['grade_transfer_factors'][name]) for metric in metrics}}
            annual.append(output)
            if case['name'] == config['reference_case']:
                for share in config['feed_treatment_shares']:
                    reference.append({'year': row['year'], 'selected_treated_share': share,
                        'source_PGM6_refined_kg_year_low_middle_high': envelope(
                            row['source_feed_kg_year'], profiles['ci'], 'PGM6', share, config, groups,
                            config['grade_transfer_factors']['ci']),
                        'lunar_REY_refined_kg_year_low_middle_high': envelope(
                            row['lunar_raw_kg_year'], profiles['kreep15386'], 'REY', share, config, groups,
                            config['grade_transfer_factors']['kreep15386'])})
    # Patch arithmetic explicitly couples a hypothetical rich lithology with a
    # separately computed bulk column. No geographical occurrence is inferred.
    lunar = read(HERE/'results/lunar_potential.json.gz')
    depth = next(r for r in lunar['bulk_inventory'] if r['depth_m'] == config['lunar_patch_depth_m'])
    patch = next(r for r in depth['footprints'] if r['area_km2'] == config['lunar_patch_area_km2'])
    masses = patch['bulk_mass_kg_p05_p50_p95']
    patch_grade = grade(profiles['kreep15386'], 'REY', groups)
    patch_factors = config['grade_transfer_factors']['kreep15386']
    patch_out = {**patch, 'depth_m': config['lunar_patch_depth_m'],
        'contained_REY_kg_low_middle_high': [mass*patch_grade*factor*1e-6 for mass, factor in zip(masses, patch_factors)],
        'refined_REY_kg_low_middle_high': [element_balance(mass, patch_grade*factor, 1,
            config['recovery_scenarios'][scenario])['refined_kg'] for mass, factor, scenario in zip(masses, patch_factors, SCENARIOS)],
        'evidence': config['lunar_patch_evidence']}
    thresholds = {}
    for name, metric in [('ci', 'PGM6'), ('ci', 'Ni'), ('kreep15386', 'REY'), ('kreep15386', 'NdPr')]:
        thresholds[name+':'+metric] = [break_even(grade(profiles[name], metric, groups),
            config['recovery_scenarios']['middle'], c) for c in config['marginal_cost_USD2041_per_kg_treated_feed']]
    # Named observations are a different evidence tier from sampled meteorites.
    candidates = read(HERE/'inputs/coproduct_candidates.json')
    limits = read(HERE/'../enclosed_mining/inputs.json')['wider_geometry_screen']['limits']
    target_rows = []
    for record in candidates['records']:
        release = protection.evaluate_release({'name': record['full_name'].strip(),
            'source_id': str(record['spkid']), 'region': 'near_earth'})
        failures = [k for k, threshold in [('a', limits['maximum_semimajor_axis_AU']),
                    ('e', limits['maximum_eccentricity']), ('i', limits['maximum_inclination_deg'])]
                    if float(record[k]) > threshold]
        target_rows.append({'record': record, 'failed_wider_orbit_axes': failures,
            'in_wider_orbital_screen': not failures,
            'release_status': release['release_status'], 'unknown_conditions': release['unknown_conditions'],
            'source_holds': release['source_holds'], 'extraction_authorized': False,
            'precious_metal_assay_available': False, 'qualified_stock_kg': None})
    return {
        'schema': 'terluna.resources.coproducts/1',
        'producer': {'path': 'resources/source_access/coproducts.py', 'sha256': sha(__file__)},
        'dependencies': {'coproduct_inputs.json': sha(HERE/'coproduct_inputs.json'), **config['dependencies']},
        'evidence': 'Sample-based composition analogues and explicit geological/engineering sensitivities. Low/middle/high are joint scenarios, not probability quantiles. No mine, plant, product specification or incremental sales is qualified.',
        'units': 'Elemental tonnes per million tonnes of original fresh feed; annual results explicitly kg/year. Monetary thresholds in constant USD2041 per kg of the specified elemental product or product group.',
        'reading_rule': 'Only distinct elements can be summed. PGM6 and REY are group totals; NdPr is within REY. Profiles, treatment shares, recoveries and physical orderbooks are alternatives. Refined outputs assume installed plant capacity, specification and material access; saleable/export quantities require further conditions.',
        'per_Mt': per_mt, 'annual_cases': annual, 'reference_share_sensitivity': reference,
        'hypothetical_kreep_patch': patch_out,
        'marginal_cost_USD2041_per_kg_treated_feed': config['marginal_cost_USD2041_per_kg_treated_feed'],
        'single_product_break_even_USD2041_per_kg_refined': thresholds,
        'metal_rich_Pt_sensitivity': [{'Pt_ppm': ppm, 'refined_Pt_tonnes_per_Mt_at_middle_recovery':
            element_balance(1e9, ppm, 1, config['recovery_scenarios']['middle'])['refined_kg']/1000}
            for ppm in config['metal_rich_pt_ppm_sensitivity']],
        'named_candidates': target_rows, 'extraction_authorized': False,
        'additional_saleable_output_kg': None, 'additional_export_revenue_USD2041': None,
        'accounting_rule': 'Recoverable potential is not automatically additional to existing generic metal or manufactured exports. Value each element once after domestic use, substitution, product qualification, logistics, market absorption and existing output credits have been reconciled.',
    }


def main():
    result = run()
    path = HERE/'results/coproducts.json.gz'
    path.write_bytes(gzip.compress(json.dumps(result, sort_keys=True, separators=(',', ':')).encode(), mtime=0))
    print(f'{len(result["annual_cases"])} checkpoint cases written to {path}')


if __name__ == '__main__':
    main()
