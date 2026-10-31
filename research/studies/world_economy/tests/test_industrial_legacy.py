"""Cross-product scale, equal-service stock and ownership controls."""
import json
import pytest
from research.studies.world_economy import industrial_legacy as m

@pytest.fixture(scope='module')
def product():return m.results()

def test_product_reproduces(product):
 assert product==json.loads((m.HERE/'results/industrial_legacy.json').read_text())

def test_wrong_source_bytes_are_rejected(tmp_path):
 s,_=m.load_inputs();s['products'][next(iter(s['products']))]['sha256']='0'*64
 p=tmp_path/'bad.json';p.write_text(json.dumps(s))
 with pytest.raises(ValueError,match='Input changed'):m.load_inputs(p)

def test_common_service_and_no_repeat_manufacturing_product(product):
 for row in product['normalized_scenarios']:
  comparator=[c['totals']['pv_all_just_in_time_U'] for c in row['cases']]
  assert max(comparator)==pytest.approx(min(comparator))
  for c in row['cases']:
   t=c['totals']
   assert t['stock_withdrawals_kg']+t['just_in_time_shortfall_purchases_kg']==pytest.approx(50)
   assert t['original_manufacture_gva_recounted_on_withdrawal_U']==0
   assert t['terminal_liquidation_credit_in_pv_U']==0

def test_durable_and_obsolescent_reserves_have_different_useful_horizons(product):
 r=product['normalized_reference']
 assert r['durable_feed']['cheapest_tested_reserve_years']==50
 assert r['modular_parts']['cheapest_tested_reserve_years']==20
 assert r['rapidly_changing_hardware']['cheapest_tested_reserve_years']==5
 high=next(c for c in r['rapidly_changing_hardware']['cases'] if c['reserve_years']==50)
 assert high['totals']['terminal_incompatible_kg']>0
 assert high['totals']['pv_saving_vs_just_in_time_U']<0

def test_price_competition_can_remove_the_reason_to_stockpile(product):
 row=next(r for r in product['normalized_scenarios'] if r['stock_class']=='durable_feed' and r['future_buy_cost_per_kg']==1 and r['annual_price_decline']==0 and r['real_discount_rate']==.005)
 assert row['cheapest_tested_reserve_years']==0

def test_film_feed_is_unrecovered_fraction_not_all_turnover(product):
 rows=product['film_fresh_feed_reserves']
 for r in rows:
  assert r['annual_virgin_feed_kg']==pytest.approx(r['gross_renewal_Gt_year']*1e12*r['fresh_fraction'])
  assert r['loss_free_50_year_kg']==pytest.approx(r['annual_virgin_feed_kg']*50)
  assert r['loss_free_50_year_kg']<r['attrition_adjusted_initial_kg']<r['loss_free_50_year_kg']*1.02
 assert min(r['attrition_adjusted_initial_kg'] for r in rows)==pytest.approx(.11693482970010536e12)
 assert max(r['attrition_adjusted_initial_kg'] for r in rows)==pytest.approx(2.348825206362446e12)

def test_asset_return_deducts_replacement_before_owner_income(product):
 for r in product['installed_asset_examples']:
  assert r['gross_operating_surplus_U_year']-r['cfc_and_replacement_U_year']==pytest.approx(r['net_operating_surplus_before_financing_U_year'])
  assert r['lunar_owner_income_before_financing_U_year']<=r['net_operating_surplus_before_financing_U_year']

def test_construction_saving_applies_only_to_eligible_package(product):
 rows=product['construction_package_sensitivities']
 assert any(r['net_cost_reduction_U']<0 for r in rows)
 for r in rows:
  assert r['gross_avoided_cost_U']==pytest.approx(r['programme_cost_U']*r['eligible_share']*r['eligible_package_saving_fraction'])
  assert r['net_cost_reduction_U']==pytest.approx(r['gross_avoided_cost_U']-r['extra_enclosure_transition_cost_U'])
  assert r['independent_export_receipt_U']==0
