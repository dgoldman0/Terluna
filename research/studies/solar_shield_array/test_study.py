import numpy as np
import pytest
import hashlib
import json
from pathlib import Path
from shared import constants as K
from .run import constant_coeff, CORE_DIVERTED, average
from .aperture import disk_nodes
from .run import HERE, ROOT
from shared.provenance import constants_used,constants_changed


def test_climate_momentum_budget_uses_bolometric_flux():
    assert CORE_DIVERTED==pytest.approx((1361-1235*0.95)/1361)
    assert constant_coeff(78000).shape==(3,7)


def test_area_quadrature_integrates_area_and_second_moment():
    r=7300e3;c=2100e3
    xy,w=disk_nodes(r,c,4,12)
    assert w.sum()==pytest.approx(1)
    np.testing.assert_allclose(w@xy,0,atol=1e-9)
    assert w@np.sum(xy*xy,axis=1)==pytest.approx(r*r/2)
    assert w[:48].sum()==pytest.approx((c/r)**2)


def test_time_average_preserves_irregular_samples():
    t=np.array([0.,1.,2.5,7.]);f=2*t+3
    assert average(f,t)==pytest.approx(10)


def test_exported_snapshots_bind_current_sources_and_validation():
    path=HERE/"results/holding.json"
    holding=json.loads(path.read_text())
    validation=json.loads((HERE/"results/validation.json").read_text())
    assert validation["holding_product_sha256"]==hashlib.sha256(path.read_bytes()).hexdigest()
    for relative,expected in holding["producer"]["source_hashes"].items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,relative
    assert "shared/constants.json" not in holding["producer"]["source_hashes"]
    assert not constants_changed(holding["producer"]["constants"])
    assert holding["producer"]["constants"]==constants_used(holding["producer"]["source_hashes"])
    assert holding["producer"]["exporter_sha256"]==hashlib.sha256((HERE/"publish.py").read_bytes()).hexdigest()
    assert validation["producer"]["runner_sha256"]==hashlib.sha256((HERE/"validate.py").read_bytes()).hexdigest()


def test_final_results_resolve_coverage_control_and_numerical_convergence():
    holding=json.loads((HERE/"results/holding.json").read_text())
    validation=json.loads((HERE/"results/validation.json").read_text())
    q=holding["aperture_quadrature"]
    assert q["variable_distance_finer"]["mean_power_TW"]<q["baseline"]["mean_power_TW"]<q["baseline_without_solar"]["mean_power_TW"]
    assert abs(q["variable_distance_finer"]["mean_power_TW"]/q["variable_distance"]["mean_power_TW"]-1)<0.001
    assert validation["propagation"]["feedback"]["minimum_full_sun_coverage_margin_m"]>0
    assert validation["propagation"]["solar_only"]["first_coverage_failure_hours"] is not None
    assert validation["nodal_eclipses"]["events"]
    assert validation["worst_eclipse"]["umbra_duration_hours"]>0
