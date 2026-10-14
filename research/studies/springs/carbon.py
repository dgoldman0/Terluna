"""What the groundwater means for the carbon cycle: the carbon the filling crust takes down, the weathering the springs
carry and the calcite they lay in the rock, against the carbon cycle's own figures.

    python -m research.studies.springs.carbon       # -> results/springs_carbon.json, seconds

A light reader of the springs product, stage 1's uptake and the carbon cycle, kept apart from run.py so that it can be
rebuilt in seconds whenever the carbon cycle's product is. While the crust fills, the water it draws from the land, the
lakes and the seas carries their dissolved inorganic carbon down with it: the land's at the water table (the springs
product, band by band), the lakes' and seas' as the carbon cycle's central held run has them through time. Beyond the
carbon cycle's 10,000 years its last waters are held.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'springs_carbon.json'
SCHEMA = 'terluna.research.springs-carbon/1'
SPRINGS = 'research/studies/springs/results/springs.json'
STAGE1 = 'geography/results/groundwater_uptake.json'
CYCLE = 'research/studies/carbon_cycle/results/carbon_cycle.json'
INPUTS = {SPRINGS: 'terluna.research.springs/1', STAGE1: 'terluna.geography.groundwater-uptake/1',
          CYCLE: 'terluna.research.carbon-cycle/1'}
CODE = ('research/studies/springs/carbon.py',)
RUN = 'central/held/natural'
M_C = 12.011
SPANS = (('0', '10'), ('10', '100'), ('100', '1000'), ('1000', '10000'), ('10000', '100000'))
STATE_OF_SPAN = {'10': 'young_1000_yr', '100': 'young_1000_yr', '1000': 'young_1000_yr',
                 '10000': 'young_10000_yr', '100000': 'mature'}
EVIDENCE = ('Products multiplied together: stage 1\'s uptake by setting and band, the springs product\'s water at the '
            'water table, and the carbon cycle\'s central held run for the lakes\' and seas\' waters and its sinks. The '
            'lakes\' and seas\' chemistry is the carbon cycle\'s, which does not yet count what the floors take; the '
            'land\'s water before 1,000 years is the 1,000-year state\'s.')
READING_RULE = ('Spans of years from the moment the planned climate and the water are in place. Water in km3 a year; '
                'carbon in Gt C a year; alkalinity in eq a year. The carbon cycle\'s rates are means over each span of '
                'its central held run (methane held, magnesium dissolved), its last values held past 10,000 years.')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rnd(x, sig=4):
    from geography.groundwater_uptake import rnd as r
    return r(x, sig)


def span_mean(years, values, a, b):
    """The mean over [a, b] of a series linear between its years and held beyond its last."""
    t = np.linspace(a, b, 4001)
    return float(np.interp(t, years, values).mean())


def results():
    data = {}
    for path, schema in INPUTS.items():
        data[path] = json.loads((ROOT / path).read_text())
        if data[path].get('schema') != schema:
            raise ValueError(f'{path}: unexpected schema {data[path].get("schema")}')
    springs, stage1, cycle = data[SPRINGS], data[STAGE1], data[CYCLE]
    series = cycle['runs'][RUN]['series']
    years = np.array(series['year'], float)
    waters = series['waters']
    rates = series['rates_gt_c_yr']
    bands = {b['band']: b for b in stage1['setting']['bands']}
    out = {}
    for case, v in stage1['cases'].items():
        rows = {}
        for a, b in SPANS:
            state = STATE_OF_SPAN[b]
            sp_bands = springs['cases'][case]['states'][state]['variants']['central']['by_band']
            water_km3 = dict(seas=0.0, lakes=0.0, land=0.0)
            carbon = dict(seas=0.0, lakes=0.0, land=0.0)
            alk_seas = 0.0
            sea_dic = span_mean(years, waters['seas']['dic_mmol_kg'], float(a), float(b)) * 1e-3
            lake_dic = span_mean(years, waters['lakes']['dic_mmol_kg'], float(a), float(b)) * 1e-3
            sea_ta = span_mean(years, waters['seas']['alkalinity_meq_kg'], float(a), float(b)) * 1e-3
            for name, band in v['bands'].items():
                cols = band['columns']
                for setting, key, area in (('seas', 'sea_floor', 'sea_m2'), ('lakes', 'lake_floor', 'lake_m2'),
                                           ('land', 'land', 'land_m2')):
                    if key not in cols:
                        continue
                    m3 = cols[key]['by_year'][b]['uptake_m_yr_since_last'] * bands[name][area]
                    water_km3[setting] += m3 / 1e9
                    if setting == 'land':
                        dic = sp_bands[name]['water_table_entry']['dic_mmol_kg'] * 1e-3 if name in sp_bands else 0.0
                    else:
                        dic = sea_dic if setting == 'seas' else lake_dic
                    carbon[setting] += m3 * 1e3 * dic * M_C / 1e15
                    if setting == 'seas':
                        alk_seas += m3 * 1e3 * sea_ta
            total = sum(carbon.values())
            net = span_mean(years, rates['net_sink'], float(a), float(b))
            imp = span_mean(years, rates['import'], float(a), float(b))
            weathering = span_mean(years, rates['weathering'], float(a), float(b))
            rows[f'{a}-{b}'] = dict(
                water_km3_yr=water_km3, carbon_gt_c_yr=dict(carbon, total=total),
                carbon_cycle=dict(net_sink_gt_c_yr=net, import_gt_c_yr=imp, weathering_gt_c_yr=weathering),
                over_carbon_cycle_import=total / imp if imp > 0 else None,
                over_carbon_cycle_net_sink=total / net if net > 0 else None,
                seas_alkalinity_taken_eq_yr=alk_seas,
                seas_alkalinity_over_weathering=alk_seas / (weathering * 1e15 / M_C) if weathering > 0 else None,
                land_water_from_state=state)
        out[case] = rows
    soil_carbonate = cycle['corners']['central']['soil_carbonate_gt_c_yr']
    weathering = {}
    for state, row in springs['cases']['central']['states'].items():
        m = row['variants']['central']['moon']
        weathering[state] = dict(
            land_export_gt_c_yr=m['weathering']['total_gt_c_yr'],
            over_carbon_cycle_catchment=m['weathering']['over_carbon_cycle_catchment'],
            groundwater_share=m['weathering']['groundwater_share'],
            calcite_laid_gt_c_yr=m['calcite_laid_gt_c_yr'],
            rock_calcite_over_soil_carbonate=(m['calcite_laid_gt_c_yr']['saturated_rock']
                                              + m['calcite_laid_gt_c_yr']['unsaturated_rock']) / soil_carbonate)
    return dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
                setting=dict(carbon_cycle_run=RUN, carbon_cycle_soil_carbonate_gt_c_yr=soil_carbonate,
                             carbon_cycle_catchment_gt_c_yr=cycle['setting']['weathering_gt_c_yr_at_catchment_rate']),
                filling_crust=out, weathering_and_calcite=weathering)


def main(argv=None) -> int:
    product = rnd(results())
    product['producer'] = dict(study='research/studies/springs', files={p: sha256(ROOT / p) for p in CODE},
                               inputs={p: sha256(ROOT / p) for p in INPUTS})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    for span, row in product['filling_crust']['central'].items():
        print(span, row['carbon_gt_c_yr']['total'], 'Gt C/yr;', row['over_carbon_cycle_import'], 'of the import')
    return 0


if __name__ == '__main__':
    sys.exit(main())
