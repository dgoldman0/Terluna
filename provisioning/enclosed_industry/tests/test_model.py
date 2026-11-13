"""Independent transport controls and component/stock accounting properties."""
from copy import deepcopy
import json
import math
from pathlib import Path
import pytest
from provisioning.enclosed_industry import model as m

HERE=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def data():return json.loads((HERE/'inputs/engineering.json').read_text())

@pytest.mark.parametrize('route,trips,burned',[
    ('Earth',529.7385319490899,645960632.1573431),
    ('Moon',1388.5098881093122,293462569.4859547),
    ('source',12.777403176831827,149396259.74518904)])
def test_transport_matches_pinned_producer_controls(data,route,trips,burned):
    # Independent controls obtained from146a5582 revenue.ferry with the recorded
    # original vehicle/tank specifications and100e6 outward/20e6 return kg/year.
    got=m.ferry(1e8,2e7,m.routes(data)[route])
    assert got['trips']==pytest.approx(trips,rel=1e-10)
    assert got['burned']==pytest.approx(burned,rel=1e-10)


def test_zero_ferry_has_no_phantom_trip(data):
    assert m.ferry(0,0,m.routes(data)['source'])['burned']==0


def test_lifetime_distribution_preserves_mean():
    life=20.
    mean=sum(m.survival(age,life) for age in range(300))
    assert mean==pytest.approx(life,rel=1e-12)
    assert m.survival(20,20,None)==0


def test_two_year_orders_and_idle_retirement(data):
    rows=[dict(year=y,population=0.) for y in range(8,20)]
    targets={r['year']:({'source_computer':100.} if r['year']==10 else {}) for r in rows}
    flows,stocks,events=m.cohort_ledger(data,rows,targets,{})
    chosen=[e for e in events if e['asset']=='source_computer']
    assert len(chosen)==1
    assert chosen[0]['order_year']==8 and chosen[0]['commission_year']==10
    assert flows[8]['imports']['source']==pytest.approx(50.)
    assert flows[9]['imports']['source']==pytest.approx(50.)
    assert stocks[10]['source_computer']==pytest.approx(100.)
    assert 0<stocks[19]['source_computer']<stocks[11]['source_computer']<100


def test_initial_cohort_gross_purchases_cover_survivors(data):
    life=20.;shape=3.;phases=8;stock=100.
    denominator=sum(m.survival(j+1,life,shape) for j in range(phases))
    purchased=stock*phases/denominator
    assert purchased>stock
    assert sum(purchased/phases*m.survival(j+1,life,shape) for j in range(phases))==pytest.approx(stock)


def test_source_coproduct_feed_is_maximum(data):
    row=dict(year=8,population=12500.,source_goods=1e7,new_residential_people=0.)
    f=m.physical(data,row,m.blank_flow());z=data['industrial']['source']
    assert f['source_feed_kg_year']==pytest.approx(max(f['source_water_kg_year']/(z['water_grade']*z['water_recovery']),f['source_nitrogen_kg_year']/(z['nitrogen_grade']*z['nitrogen_recovery']),f['source_virgin_metal_kg_year']/(z['metal_grade']*z['metal_recovery'])))
    assert sum(f['source_cargo_breakdown_kg_year'].values())==pytest.approx(f['source_cargo_kg_year'])


def test_poor_grade_changes_feed_bound_heat(data):
    row=dict(year=8,population=12500.,new_residential_people=0.)
    base=m.physical(data,row,m.blank_flow());bad=deepcopy(data);bad['industrial']['source']['water_grade']/=2
    stress=m.physical(bad,row,m.blank_flow())
    assert stress['source_feed_kg_year']>base['source_feed_kg_year']
    assert stress['source_feed_process_W']>base['source_feed_process_W']
    assert stress['source_product_refining_W']==base['source_product_refining_W']


def test_settlement_target_does_not_copy_industry(data):
    assert all(data['templates'].get(n,{'owner':'settlement'})['owner']=='settlement' for n in m.residential_targets(data,250000.))
    row=dict(year=8,population=250000.,new_residential_people=0.)
    flow=m.physical(data,row,m.blank_flow());req=m.requirements(data,row,flow)
    assert req['asteroid_tugs']<data['original_architecture']['assets']['asteroid_tugs']['mass_kg']


def test_recycled_inventory_closes_material_balance(data):
    flow=m.blank_flow();flow['local']['source']=100.;flow['recycled']['source']=150.
    p=m.physical(data,dict(year=8,population=0.),flow)
    assert p['recycled_inventory_kg']['source']==50.
    assert p['source_virgin_metal_kg_year']==0.


def test_deterministic_packaging_round_trip(tmp_path):
    import gzip
    import hashlib
    from provisioning.enclosed_industry.package import package
    source=tmp_path/'source.json';a=tmp_path/'a.json.gz';b=tmp_path/'b.json.gz'
    source.write_text(json.dumps(dict(schema='test/1',producer={},dependencies={},annual=[dict(year=8,flow=3.5)])))
    manifest=package(source,a);package(source,b)
    assert a.read_bytes()==b.read_bytes()
    decoded=gzip.decompress(a.read_bytes())
    assert json.loads(decoded)==json.loads(source.read_text())
    assert hashlib.sha256(decoded).hexdigest()==manifest['decompressed_sha256']


def test_computing_adds_power_but_not_chemical_equipment(data):
    row=dict(year=8,population=12500.,new_residential_people=0.)
    normal=m.physical(data,row,m.blank_flow());compute_row=dict(row,compute=1e9)
    digital=m.physical(data,compute_row,m.blank_flow())
    a=m.requirements(data,row,normal);b=m.requirements(data,compute_row,digital)
    assert b['source_chemistry']==pytest.approx(a['source_chemistry'])
    assert b['source_solar']>a['source_solar']
    assert b['source_radiators']>a['source_radiators']


def test_seed_route_fuel_is_priced_once_without_source_feed(data):
    from provisioning.enclosed_industry.bootstrap import apply
    old=data['original_architecture'];unit=old['assets']['Earth_hull_tanks']['mass_kg']/old['engineering']['counts']['Earth_route_ships']
    p=dict(rules=[],cases=[dict(seed=dict(cost=dict(total_USD2041=10.),source_feed_kg=123.),annual=[dict(asset_stocks_kg=dict(Earth_hull_tanks=3*unit))])])
    got=apply(p,data);seed=got['cases'][0]['seed'];fuel=seed['priced_initial_route_fuel']
    assert fuel['Earth_route_ship_count']==3
    assert seed['source_feed_kg']==123.
    assert seed['cost']['total_USD2041']==pytest.approx(10.+fuel['external_propellant_purchase_USD2041']+fuel['Earth_launch_USD2041'])
    with pytest.raises(ValueError,match='already priced'):apply(got,data)
