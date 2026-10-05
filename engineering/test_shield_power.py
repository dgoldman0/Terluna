import pytest

from engineering.shield_power import routed_power, optical_requirement, exhaust_rate


def test_energy_partition_and_operation_before_delivery():
    r=routed_power(1000.,.5,.4,.8,70.,10.)
    assert r['captured_optical_W']==500.
    assert r['gross_electric_W']==200.
    assert r['ideal_collector_waste_heat_W']==300.
    assert r['delivered_W']==96.
    assert r['final_link_loss_W']==pytest.approx(24.)
    assert sum(r[k] for k in ['uncaptured_optical_W','ideal_collector_waste_heat_W',
        'propulsion_W','other_operations_W','delivered_W','final_link_loss_W'])==1000.
    assert optical_requirement(96.,.5,.4,.8,70.,10.)==1000.


def test_deficit_is_explicit_and_propellant_is_separate():
    r=routed_power(100.,.1,.2,.8,20.)
    assert r['deficit_W']==18.
    assert r['delivered_W']==0.
    assert exhaust_rate(23.835e12,30000.,.7)==pytest.approx(37076.666666666664)
