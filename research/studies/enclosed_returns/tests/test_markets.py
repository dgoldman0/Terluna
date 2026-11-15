from copy import deepcopy
import json
import pytest
from research.studies.enclosed_returns import markets as m

@pytest.fixture(scope="module")
def data():return m.inputs()


def test_schedule_closes_at_cap_and_stops_construction(data):
    cfg,*_=data
    assert m.population(8,cfg)==12500
    assert m.population(200,cfg)==pytest.approx(cfg['cap_population'])
    assert m.population(220,cfg)==m.population(200,cfg)
    construction=sum(m.population(t+2,cfg)-m.population(t+1,cfg) for t in range(8,199))
    assert construction+12500==pytest.approx(cfg['cap_population'])


def test_materials_below_known_return_cost_do_not_generate_revenue_or_feed(data):
    cfg,q,d,b,*_=data
    for name in ('raw_resource_price','commodity_refined'):
        r=m.annual_market(cfg,q,d,b,name,'breakthrough',.005,200,52,0,0)
        assert r['goods']==r['source_goods']==0
        assert r['revenue_by_sector_USD2041']['materials']==r['revenue_by_sector_USD2041']['source_products']==0
        assert r['refused_export_offers']


def test_exports_follow_customers_and_capture_not_resident_multiplier(data):
    cfg,q,d,b,*_=data
    one=m.annual_market(cfg,q,d,b,'specialty','breakthrough',.005,100,30,0,0)
    larger=dict(cfg,seed_population=cfg['seed_population']*2,cap_population=cfg['cap_population']*2)
    two=m.annual_market(larger,q,d,b,'specialty','breakthrough',.005,100,30,0,0)
    for field in ('goods','source_goods','compute','propellant_sales'):
        assert one[field]==two[field]


def test_customer_budgets_share_transactions_once(data):
    cfg,q,d,b,*_=data
    for name in q['strategies']:
        r=m.annual_market(cfg,q,d,b,name,'breakthrough',.01,200,120,0,0)
        v=r['revenue_by_sector_USD2041'];tol=1.000000001
        assert v['materials']+v['source_products']<=r['industrial_budget_USD2041']*tol
        assert v['sports']+v['entertainment']+v['licenses']+v['destination']<=r['household_budget_USD2041']*tol
        assert v['propellant']+v['mission_support']+r['mission_research_USD2041']<=r['mission_budget_USD2041']*tol
        assert v['research']-r['mission_research_USD2041']+v['professional']+v['computing']<=r['institutional_budget_USD2041']*tol
        assert sum(v.values())==pytest.approx(r['gross_USD2041'])


def test_service_growth_builds_extra_facilities_and_capital(data):
    cfg,q,d,b,*_=data
    c=m.case(cfg,q,d,b,'services','breakthrough','reference',.005)
    end=next(r for r in c['annual'] if r['year']==200)
    assert end['service_kit_count']>20
    assert c['annual'][0]['service_kit_count']==1
    for one,two in zip(c['annual'],c['annual'][1:]):
        assert 0<=two['service_kit_count']-one['service_kit_count']<=1
    assert end['service_installed_capital_USD2041']==b['cultural_platform']['capital_USD2041']+end['service_kit_count']*b['service_kit']['capital_USD2041']


def test_sustaining_check_holds_markets_and_population_after_completion(data):
    cfg,q,d,b,*_=data
    c=m.case(cfg,q,d,b,'diversified','breakthrough','reference',.005)
    rows={r['year']:r for r in c['annual']}
    for field in ('population','goods','source_goods','propellant_sales','compute','service_kit_count','gross_USD2041'):
        assert rows[220][field]==pytest.approx(rows[200][field])
    assert rows[220]['new_residential_people']==rows[200]['new_residential_people']==0


def test_product_reproduces_sources_and_all_cases():
    saved=json.loads((m.HERE/'results/orders.json').read_text())
    assert m.product()==saved
    assert len(saved['cases'])==81
    assert saved['producer']['sha256']==m.digest(m.__file__)
    for path,h in saved['dependencies'].items():assert m.digest(m.HERE/path)==h
