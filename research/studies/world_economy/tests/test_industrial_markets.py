"""Demand, resource, competitive-cost and lifetime-finance controls."""
import copy
import json
import pytest
from research.studies.world_economy import industrial_markets as m

@pytest.fixture(scope='module')
def product():return m.results()

@pytest.fixture
def setup():
 s,p=m.load_inputs();return s,m.anchors(p),s['plants'][1],s['cost_cases']['modest_site_advantage']

def test_product_reproduces(product):
 assert product==json.loads((m.HERE/'results/industrial_markets.json').read_text())

def test_bound_source_change_rejected(tmp_path):
 s,_=m.load_inputs();s['products'][next(iter(s['products']))]['sha256']='0'*64
 p=tmp_path/'bad.json';p.write_text(json.dumps(s))
 with pytest.raises(ValueError,match='Input changed'):m.load_inputs(p)

def test_offer_must_beat_both_competitors(setup):
 s,b,p,c=setup
 for f in s['families']:
  for buyer in s['buyers']:
   r=m.offer(s,f,c,buyer)
   if r['cost_competitive']:
    assert r['delivered_lunar_cost_U_kg']<r['offered_delivered_price_U_kg']<min(r['competitor_costs_U_kg'].values())
   assert r['producer_receipt_U_kg']+r['delivery_U_kg']==pytest.approx(r['offered_delivered_price_U_kg'])

def test_no_advantage_and_orbital_learning_have_no_cost_based_orders(product):
 for r in product['scenarios']:
  if r['cost_case'] in ('parity','orbital_learning'):
   assert r['totals']['produced_kg']==0
   assert r['totals']['cash_assigned_U']==0
   assert r['totals']['unfunded_cash_gap_U']>0
   assert r['unlevered_project_npv_U']<0

def test_reference_budgets_orders_and_finance_close(product,setup):
 s,b,p,c=setup;r=product['reference'];normal=r['normal_saleable_capacity_kg_year']
 for row in r['annual']:
  assert row['produced_kg']<=normal+1e-3
  assert row['produced_kg']<=p['outbound_route_kg_year']+1e-3
  assert row['purchased_energy_J']<=p['power_W']*m.JULIAN_DAY*m.JULIAN_YEAR_DAYS+1
  assert abs(row['cash_closure_U'])<r['initial_capital_U']*1e-12
  for buyer,spend in row['buyer_spending_U'].items():
   assert spend<=r['buyer_budgets_U_year'][buyer]+1
  assert row['independent_cash_assigned_U']<=row['cash_assigned_U']
 for buyer,stock in r['remaining_eligible_buyer_stock_kg'].items():
  for key,value in stock.items():
   sold=sum(x['buyer_allocations_kg'][buyer][key] for x in r['annual'])
   assert value+sold==pytest.approx(r['initial_eligible_buyer_stock_kg'][buyer][key])
 assert r['totals']['local_feed_kg']+r['local_feed_remaining_kg']==pytest.approx(s['local_feed_allowance_kg'])

def test_lifetime_cash_including_replacement_and_closure(setup):
 s,b,p,c=setup;r=m.run_case(s,b,p,c,50,keep_annual=True);t=r['totals'];close=r['closure']
 variable=sum(sum(row['made_by_family_kg'][f['id']]*f['base_production_cost_U_kg']*c['moon']*.7 for f in s['families']) for row in r['annual'])
 sources=r['initial_sponsor_equity_U']+r['initial_borrowing_U']+t['producer_receipts_U']+t['unfunded_cash_gap_U']+close['assumed_salvage_U']
 uses=r['initial_capital_U']+variable+r['fixed_annual_operations_U']*50+t['replacement_capex_U']+t['loan_interest_U']+t['loan_principal_U']+close['remaining_loan_to_settle_U']+close['decommission_cost_U']+t['cash_assigned_U']+close['terminal_cash_after_obligations_U']
 assert sources==pytest.approx(uses,rel=1e-12)
 assert t['replacement_capex_U']==pytest.approx(r['initial_capital_U'])
 assert close['replacement_reserve_released_U']==pytest.approx(r['initial_capital_U']/4)

def test_positive_interim_cash_is_not_positive_npv(product):
 r=product['reference'];assert r['totals']['cash_assigned_U']>0
 assert r['unlevered_project_npv_U']<0
 assert 0<r['closure']['zero_npv_terminal_recovery_fraction_of_remaining_book']<1

def test_independent_buyer_priority_retains_budgets_and_changes_financing_mix(product):
 private=product['reference'];public=product['independent_buyer_priority']
 assert public['totals']['independent_cash_assigned_U']>private['totals']['independent_cash_assigned_U']
 assert public['totals']['cash_assigned_U']<=private['totals']['cash_assigned_U']
 for buyer,row in public['buyer_totals'].items():
  assert row['spending_U']<=public['buyer_budgets_U_year'][buyer]*public['years']+1

def test_oversized_factory_is_demand_limited(product):
 r=next(x for x in product['scenarios'] if x['plant']=='large' and x['cost_case']=='modest_site_advantage' and x['years']==50)
 assert r['actual_normal_capacity_utilization']<.2
 assert r['totals']['unfunded_cash_gap_U']>0

def test_power_route_and_reserves_each_bind(setup):
 s,b,p,c=setup
 low_power=dict(p,power_W=1e6);r=m.run_case(s,b,low_power,c,10,keep_annual=True)
 assert all(x['purchased_energy_J']<=1e6*m.JULIAN_DAY*m.JULIAN_YEAR_DAYS+1 for x in r['annual'])
 low_route=dict(p,outbound_route_kg_year=1e6);r=m.run_case(s,b,low_route,c,10)
 assert r['totals']['produced_kg']<=1e7+1
 s['local_feed_allowance_kg']=1e6;r=m.run_case(s,b,p,c,10)
 assert r['totals']['local_feed_kg']<=1e6+1e-6

def test_no_replacement_at_known_shutdown_and_no_expired_debt_payment(setup):
 s,b,p,c=setup;r=m.run_case(s,b,p,c,40,keep_annual=True)
 assert r['totals']['replacement_capex_U']==0
 assert r['closure']['remaining_book_U']==0
 r=m.run_case(s,b,p,c,101,keep_annual=True)
 assert r['annual'][-1]['loan_interest_U']==r['annual'][-1]['loan_principal_U']==0
 assert r['annual'][-1]['replacement_reserve_U']>0

def test_terminal_recovery_sensitivity_does_not_change_past_production(product):
 rows=product['terminal_recovery_sensitivities']
 assert [r['unlevered_npv_U'] for r in rows]==sorted(r['unlevered_npv_U'] for r in rows)
 assert rows[0]['unlevered_npv_U']<0<rows[-1]['unlevered_npv_U']

def test_post_atmosphere_cost_changes_are_product_specific(product):
 r=product['post_atmosphere_cost_thresholds']['modest_site_advantage']
 assert not r['structures']['array']['cost_competitive']
 assert r['sealed_specialty']['array']['cost_competitive']
