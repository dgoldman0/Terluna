"""Material recovery and equipment refurbishment in the vacuum buildout's industry.

Two measures act on the industry model at its three points. Closure raises the recovery of water, household goods
and retired structure, and adds the equipment and energy that sorting and treating them take. Refurbishment
collects retired machines, inspects them a year later, and returns the accepted ones to service with new repair
parts in place of a share of their mass, in workshops at the place the machine works. Each measure is costed inside
the model's own balance: its equipment, energy, consumables and labour are sized with everything else.

The accounts at the end report what the measures did: the water streams with their makeup, the retired structure
and goods with what was collected, accepted and used, and the machines withdrawn for repair.
"""
from __future__ import annotations
from copy import deepcopy
import math
from . import model

# The named sets of measures: (closure, refurbishment).
OPTIONS = {'none': (False, False), 'closure': (True, False), 'refurbishment': (False, True), 'combined': (True, True)}


def check(spec):
    if spec['refurbishment']['inspection_delay_years'] < 1:
        raise ValueError('Refurbishment needs at least a year between retirement and inspection')
    for group in ('goods', 'structural', 'refurbishment'):
        for key in ('collection', 'quality', 'repair_part_mass_fraction'):
            if key in spec[group] and not 0 <= spec[group][key] <= 1:
                raise ValueError('Recovery fractions lie between 0 and 1')
    return spec


def water_account(spec, population, guest_years=0., improved=False):
    """The water streams of the residents and their guests, with what each recovers and what must be made up."""
    scale = (population + guest_years) / spec['normalization_population']
    streams = {}
    for name, w in spec['water'].items():
        gross = w['gross_kg_year'] * scale
        r = w['improved_recovery' if improved else 'recovery']
        streams[name] = dict(gross_kg=gross, recovered_kg=gross * r, makeup_kg=gross * (1 - r), loss_kg=gross * (1 - r), recovery=r)
    extra = spec['landscape_agriculture_allowance_kg_year'] * scale
    return dict(streams=streams, landscape_agriculture_makeup_and_loss_kg=extra,
                gross_processed_kg=sum(v['gross_kg'] for v in streams.values()),
                recovered_kg=sum(v['recovered_kg'] for v in streams.values()),
                operational_makeup_kg=extra + sum(v['makeup_kg'] for v in streams.values()))


def configure(data, spec, closure, refurbishment):
    """The engineering coefficients with the measures' recoveries and their equipment templates."""
    d = deepcopy(data)
    if closure:
        d['demand']['yearly_water'] = water_account(spec, spec['normalization_population'], improved=True)['operational_makeup_kg']
        g = spec['goods']
        factor = (1 - g['collection'] * g['quality']) / (1 - g['baseline_recovery'])
        d['demand']['goods_virgin'] *= factor
        d['settlement']['operating_import_kg_year']['everyday_carbon_feed'] *= factor
        d['industrial']['lunar']['recovery'] = d['economy']['production']['source_recovery'] = spec['structural']['collection'] * spec['structural']['quality']

    def template(name, dest, life, structure_life):
        a = deepcopy(d['templates']['fabrication_cells'])
        a.update(name=name, mass=0., destination=dest, life=life, structural_life=structure_life)
        d['templates'][name] = a
    if closure:
        w = spec['water_equipment']
        template('loop_water_equipment', 'surface', w['active_life'], w['structure_life'])
        for dest in ('surface', 'source'):
            template('loop_sort_' + dest, dest, 15, 30)
    if refurbishment:
        r = spec['refurbishment']
        for dest in ('surface', 'source', 'depot'):
            template('loop_repair_' + dest, dest, r['workshop_active_life'], r['workshop_structural_life'])
    return d


def retirements(data, rows, seed, events):
    """The mass of each component that expires in each year, from the seed cohorts and the dated orders."""
    cats = model.catalogue(data)
    start = rows[0]['year']
    end = rows[-1]['year']
    shape = data.get('lifetime_shape', 3.)
    out = {y: {} for y in range(start, end + 1)}
    births = {}
    for n, a in cats.items():
        for part, share, life in [('active', 1 - a['structural_share'], a['life']), ('structure', a['structural_share'], a['structural_life'])]:
            if not share:
                continue
            phases = max(1, min(8, math.ceil(life)))
            den = sum(model.survival(j + 1, life, shape) for j in range(phases))
            births[n, part] = [(start - j - 1, seed.get(n, 0.) * share / den) for j in range(phases)]
    for ev in events:
        births[ev['asset'], ev['component']].append((ev['commission_year'], ev['mass_kg']))
    for (n, part), cohorts in births.items():
        life = cats[n]['life' if part == 'active' else 'structural_life']
        curve = [model.survival(a, life, shape) for a in range(end + 2)]
        values = [0.] * (end + 1)
        for birth, mass in cohorts:
            for y in range(max(start + 1, birth + 1), end + 1):
                age = y - birth
                values[y] += mass * (curve[age - 1] - curve[age])
        for y in range(start, end + 1):
            out[y][n, part] = values[y]
    return out


def empty_repair():
    return dict(retired_kg=0., collected_kg=0., uncollected_kg=0., inspection_input_kg=0., quality_reject_kg=0., accepted_kg=0.,
                pending_inspection_kg=0., qualified_inventory_kg=0., withdrawn_for_repair_kg=0., retained_old_parts_kg=0.,
                discarded_old_parts_kg=0., repair_parts_imported_kg=0., refurbished_commissioned_kg=0., repair_WIP_kg=0.,
                imported_mass_avoided_kg=0., import_purchase_avoided_USD2026=0., mass_balance_residual_kg=0.)


def refurbish(cats, spec, rows, flows, events, retired):
    """Fill orders for eligible machines from retired, accepted machines at the same place, and lower the imports.

    cats is the model's catalogue of assets. The flows and the events are changed in place; the return is the yearly
    account of each destination's workshop.
    """
    cfg = spec['refurbishment']
    start = rows[0]['year']
    end = rows[-1]['year']
    lag = cfg['inspection_delay_years']
    audits = {y: {k: empty_repair() for k in model.DESTS} for y in range(start - 2, end + 1)}
    for y, f in flows.items():
        f['loop_inspected'] = {k: 0. for k in model.DESTS}
        f['loop_repaired'] = {k: 0. for k in model.DESTS}
        f['loop_consume_kg'] = {k: 0. for k in model.DESTS}
        f['loop_import_avoided'] = {}
    selected = {n: [] for n in cfg['eligible_assets']}
    for ev in events:
        if ev['component'] == 'active' and ev['asset'] in selected:
            selected[ev['asset']].append(ev)
    for n, orders in selected.items():
        a = cats[n]
        dest = a['destination']
        byyear = {}
        for ev in orders:
            byyear.setdefault(ev['order_year'], []).append(ev)
        bank = 0.
        pending = 0.
        wip = 0.
        for y in range(start - 2, end + 1):
            v = audits[y][dest]
            opening_bank = bank
            opening_pending = pending
            opening_wip = wip
            retired_now = retired.get(y, {}).get((n, 'active'), 0.)
            inspected = retired.get(y - lag, {}).get((n, 'active'), 0.) * cfg['collection']
            accepted = inspected * cfg['quality']
            rejected = inspected - accepted
            collected = retired_now * cfg['collection']
            bank += accepted
            pending += collected - inspected
            used = 0.
            for ev in byyear.get(y, []):
                use = min(bank, ev['mass_kg'])
                bank -= use
                used += use
                k = cfg['repair_part_mass_fraction']
                avoided = use * (1 - k)
                growth = (ev['mass_kg'] - ev['replacement_kg']) / ev['mass_kg']
                # The order would have been bought from Earth and made locally in these shares; repair spares both.
                made = ev.get('locally_made_share', 0.)
                bought = 1 - made
                labour = ev.get('local_labour_fraction', 0.)
                for wy in (y, y + 1):
                    if wy not in flows:
                        continue
                    f = flows[wy]
                    m = avoided / 2
                    money = m * a['price'] * ev.get('price_index', 1.0)
                    f['imports'][dest] -= m * bought
                    f['capital_imports'][dest] -= m * growth * bought
                    f['external_capital'] -= money * growth * bought
                    f['external_replacement'] -= money * (1 - growth) * bought
                    if a['owner'] == 'industry':
                        f['industry_replacement_import_kg'] -= m * (1 - growth) * bought
                    for owner in ('service', 'hospitality'):
                        if a['owner'] == owner:
                            f[owner + '_external_capital'] -= money * growth * bought
                            f[owner + '_external_replacement'] -= money * (1 - growth) * bought
                    f['local'][dest] -= m * made
                    f['capital_local'][dest] -= m * growth * made
                    f['local_working_capital'] -= money * growth * made * labour
                    f['local_working_replacement'] -= money * (1 - growth) * made * labour
                    f['loop_repaired'][dest] += use / 2
                    f['loop_import_avoided'][n] = f['loop_import_avoided'].get(n, 0.) + m * bought
                    av = audits[wy][dest]
                    av['repair_parts_imported_kg'] += use * k / 2 * bought
                    av['imported_mass_avoided_kg'] += m * bought
                    av['import_purchase_avoided_USD2026'] += money * bought
                cy = y + 2
                if cy in audits:
                    audits[cy][dest]['refurbished_commissioned_kg'] += use
                ev['refurbished_mass_kg'] = use
                ev['new_imported_mass_kg'] = (ev['mass_kg'] - avoided) * bought
            # Arrivals are scheduled for this asset alone; a destination's totals hold other assets as well.
            kit = sum(ev.get('refurbished_mass_kg', 0.) * cfg['repair_part_mass_fraction'] / 2 for ev in orders if ev['order_year'] in (y, y - 1))
            commissioned = sum(ev.get('refurbished_mass_kg', 0.) for ev in orders if ev['commission_year'] == y)
            retained = used * (1 - cfg['repair_part_mass_fraction'])
            discarded = used - retained
            wip += retained + kit - commissioned
            values = dict(retired_kg=retired_now, collected_kg=collected, uncollected_kg=retired_now - collected, inspection_input_kg=inspected,
                          quality_reject_kg=rejected, accepted_kg=accepted, pending_inspection_kg=pending, qualified_inventory_kg=bank,
                          withdrawn_for_repair_kg=used, retained_old_parts_kg=retained, discarded_old_parts_kg=discarded, repair_WIP_kg=wip)
            for key, value in values.items():
                v[key] += value
            v['mass_balance_residual_kg'] += opening_bank + opening_pending + opening_wip + retired_now + kit - (bank + pending + wip + retired_now - collected + rejected + discarded + commissioned)
            flows[y]['loop_inspected'][dest] += inspected
    for y, f in flows.items():
        for dest in model.DESTS:
            mass = (f['loop_inspected'][dest] + f['loop_repaired'][dest]) * cfg['consumables_fraction_of_processed_mass']
            f['loop_consume_kg'][dest] = mass
            f['imports'][dest] += mass
    return audits


class Recovery:
    """The measures the industry model calls at its three points, and the record of what they did."""

    def __init__(self, spec, closure=True, refurbishment=True, residual_water_allowance=True):
        self.spec = check(deepcopy(spec))
        if not residual_water_allowance:
            self.spec['landscape_agriculture_allowance_kg_year'] = 0.
        self.closure = closure
        self.refurbishment = refurbishment
        self.flows = self.retired = self.audit = self.last_physical = None

    @property
    def active(self):
        return self.closure or self.refurbishment

    def configure(self, data):
        return configure(data, self.spec, self.closure, self.refurbishment)

    def after_ledger(self, data, rows, seed, flows, events):
        retired = retirements(data, rows, seed, events)
        cats = model.catalogue(data)
        audit = refurbish(cats, self.spec, rows, flows, events, retired) if self.refurbishment else None
        for y, f in flows.items():
            f['loop_structural_retired'] = {k: 0. for k in (*model.DESTS, 'bulk')}
            for (n, part), mass in retired.get(y, {}).items():
                if part == 'structure':
                    a = cats[n]
                    key = 'bulk' if a['material'] == 'bulk' else a['destination']
                    f['loop_structural_retired'][key] += mass
        self.flows = flows
        self.retired = retired
        self.audit = audit

    def after_physical(self, data, row, flow, out):
        if not self.active:
            self.last_physical = deepcopy(out)
            return
        extra = {k: 0. for k in model.DESTS}
        spec = self.spec
        if self.closure:
            water = water_account(spec, row['population'], row.get('visitor_nights_year', 0.) / data['constants']['JULIAN_YEAR_DAYS'], True)
            extra['surface'] += water['gross_processed_kg'] * spec['water_equipment']['extra_kWh_per_kg_gross']
            retired = flow.get('loop_structural_retired', {})
            s = spec['structural']
            extra['surface'] += retired.get('bulk', 0.) * s['bulk_extra_sort_kWh_kg']
            for dest in model.DESTS:
                extra[dest] += retired.get(dest, 0.) * s['metal_extra_sort_kWh_kg']
            # Structural metal in orbit already travels to the source; its sorting is done there.
            extra['source'] += extra['depot'] + extra['leo']
            extra['depot'] = extra['leo'] = 0.
            g = spec['goods']
            extra['surface'] += row['population'] * g['gross_kg_person_year'] * (g['collection'] * g['quality'] - g['baseline_recovery'])
        if self.refurbishment:
            r = spec['refurbishment']
            for dest in model.DESTS:
                extra[dest] += flow.get('loop_inspected', {}).get(dest, 0.) * r['inspection_kWh_kg'] + flow.get('loop_repaired', {}).get(dest, 0.) * r['remanufacture_kWh_kg']
            out['operating_material_USD2026'] += sum(flow.get('loop_consume_kg', {}).values()) * data['industrial']['prices']['other_special_USD2026_kg']
        watts = {k: v * 1000 / model.HOURS for k, v in extra.items()}
        out['lunar_process_W'] += watts['surface']
        out['ground_firm_W'] += watts['surface'] * data['industrial']['economics']['capacity_margin'] / data['industrial']['lunar']['process_availability']
        out['source_process_W'] += watts['source']
        out['source_heat_W'] += watts['source']
        out['depot_W'] += watts['depot']
        out['loop_extra_process_kWh_year'] = extra
        out['loop_repair_input_kg_year'] = {k: flow.get('loop_inspected', {}).get(k, 0.) + flow.get('loop_repaired', {}).get(k, 0.) for k in model.DESTS}
        out['loop_structural_retired_kg_year'] = flow.get('loop_structural_retired', {})
        out['loop_repair_consumables_kg_year'] = dict(flow.get('loop_consume_kg', {}))
        self.last_physical = deepcopy(out)

    def after_requirements(self, data, row, physical, out):
        if not self.active:
            return
        l = data['industrial']['lunar']
        margin = data['industrial']['economics']['capacity_margin']

        def machine(name, flow, rate, mass):
            out[name] = math.ceil(flow * margin / (rate * model.HOURS * l['process_availability'])) * mass
        if self.closure:
            w = self.spec['water_equipment']
            wa = water_account(self.spec, row['population'], row.get('visitor_nights_year', 0.) / data['constants']['JULIAN_YEAR_DAYS'], True)
            machine('loop_water_equipment', wa['gross_processed_kg'], w['train_kg_hr'], w['hardware_kg'])
            ret = physical.get('loop_structural_retired_kg_year', {})
            # Bulk is sorted by a stated share of a bulk train; metal by fabrication cells, as a stand-in for a sorter.
            surface = ret.get('surface', 0.) + ret.get('bulk', 0.) * self.spec['structural']['bulk_processing_fraction_of_existing_train'] * l['fabrication_cell_kg_hr'] / l['bulk_train_kg_hr']
            g = self.spec['goods']
            surface += row['population'] * g['gross_kg_person_year'] * (g['collection'] * g['quality'] - g['baseline_recovery'])
            machine('loop_sort_surface', surface, l['fabrication_cell_kg_hr'], l['fabrication_cell_hardware_kg'])
            machine('loop_sort_source', sum(ret.get(k, 0.) for k in ('source', 'depot', 'leo')), l['fabrication_cell_kg_hr'], l['fabrication_cell_hardware_kg'])
        if self.refurbishment:
            r = self.spec['refurbishment']
            for dest in ('surface', 'source', 'depot'):
                machine('loop_repair_' + dest, physical['loop_repair_input_kg_year'][dest], r['workshop_kg_hr'], r['workshop_hardware_kg'])
            # A depot workshop's electricity all leaves as heat; the depot's chemical plant keeps part as chemical energy.
            kw = physical.get('loop_extra_process_kWh_year', {}).get('depot', 0.) * 1000 / model.HOURS
            p = data['industrial']['power']
            z = data['industrial']['source']
            flux = p['emissivity'] * model.STEFAN_BOLTZMANN * p['radiator_K']**4
            out['depot_radiators'] += kw * z['hydrogen_HHV_kWh_kg'] / z['electrolysis_liquefaction_kWh_kg_H2'] * margin / flux * p['radiator_kg_m2']

    def material_account(self, data, row, previous_bank):
        """One year's retired structure and goods, and the machines sent to repair."""
        year = row['year']
        f = self.flows[year]
        phys = row['physical']
        ret = self.retired[year]
        cats = model.catalogue(data)
        structural = {}
        endbank = phys['recycled_inventory_kg']
        cfg = self.spec['structural']
        for dest in (*model.DESTS, 'bulk'):
            gross = sum(m for (n, part), m in ret.items() if part == 'structure' and ('bulk' if cats[n]['material'] == 'bulk' else cats[n]['destination']) == dest)
            recovered = f['bulk_recycled'] if dest == 'bulk' else f['recycled'][dest]
            collected = gross * cfg['collection'] if self.closure else recovered
            opening = previous_bank.get(dest, 0.)
            closing = endbank[dest]
            used = opening + recovered - closing
            structural[dest] = dict(retired_kg=gross, collected_kg=collected, recovered_kg=recovered,
                                    collection_loss_kg=gross - collected, quality_reject_kg=collected - recovered,
                                    opening_inventory_kg=opening, used_kg=used, closing_inventory_kg=closing,
                                    mass_balance_residual_kg=opening + gross - (closing + used + gross - recovered))
        g = self.spec['goods']
        gross = g['gross_kg_person_year'] * row['population']
        recovery = g['collection'] * g['quality'] if self.closure else g['baseline_recovery']
        collected = gross * g['collection'] if self.closure else gross * recovery
        goods = dict(gross_kg=gross, collected_kg=collected, recovered_kg=gross * recovery, collection_loss_kg=gross - collected,
                     quality_reject_kg=collected - gross * recovery, loss_and_virgin_makeup_kg=gross * (1 - recovery),
                     imported_carbon_feed_kg=gross * (1 - recovery) * g['carbon_share'], local_makeup_kg=gross * (1 - recovery) * (1 - g['carbon_share']),
                     interpretation='Inherited steady annual goods-return assumption; distinct from capital scrap stocks.')
        active_total = sum(m for (n, p), m in ret.items() if p == 'active')
        refurbishment_audit = self.audit[year] if self.refurbishment else {k: empty_repair() for k in model.DESTS}
        eligible_retired = sum(v['retired_kg'] for v in refurbishment_audit.values())
        return dict(structural=structural, goods=goods, active_retired_total_kg=active_total,
                    active_outside_refurbishment_scope_kg=active_total - eligible_retired, refurbishment=refurbishment_audit,
                    active_rejects_are_uncredited=True)
