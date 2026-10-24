"""Independent objective, selection and native-replay checks for the tuner."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

from atmosphere.snow import mcsnow_tuning as tuning


def synthetic_batch():
    entries = [dict(temperature_C=t, initial_age_min=3., axis_scale=1.) for t in (-14.4,-5.3)]
    batch = SimpleNamespace(entries=entries, rows=[(i,a) for i in range(2) for a in (3.,5.,10.,20.,30.)])
    batch.values = np.column_stack([tuning.targets(batch), np.zeros((10,4))])
    batch.evaluate = lambda x: dict(rows=batch.values.tolist(), error=None)
    return batch


def test_objective_excludes_prescribed_start_withheld_ages_and_speed():
    batch=synthetic_batch()
    assert np.max(abs(tuning.residual(batch,tuning.BASE))) == 0
    # Alter every withheld-age value, every supplied initial state and all speeds.
    for j,(_,age) in enumerate(batch.rows):
        if age in (3.,30.): batch.values[j,:3] *= 100
    batch.values[:,3] *= 20
    assert np.max(abs(tuning.residual(batch,tuning.BASE))) == 0


def test_objective_uses_log_ratios_and_equal_case_and_quantity_weight():
    batch=synthetic_batch()
    for j,(i,age) in enumerate(batch.rows):
        if i==0 and age in tuning.TRAIN_AGES: batch.values[j,0] *= 2
    loss=np.sum(tuning.residual(batch,tuning.BASE)**2)
    assert loss == pytest.approx(np.log(2)**2/6)
    # Equal multiplicative over- and underprediction have equal objective weight.
    for j,(i,age) in enumerate(batch.rows):
        if i==0 and age in tuning.TRAIN_AGES: batch.values[j,0] /= 4
    assert np.sum(tuning.residual(batch,tuning.BASE)**2) == pytest.approx(loss)


def test_failed_trajectory_cannot_improve_the_objective():
    batch=synthetic_batch()
    batch.evaluate=lambda x:dict(rows=[],error='native failure')
    assert np.all(tuning.residual(batch,tuning.BASE)==100)


@pytest.mark.parametrize('bad',[[0]*5,[0]*7,[float('nan')]*6,[-.1]*6,[1.1]*6])
def test_nonphysical_search_coordinates_rejected(bad):
    with pytest.raises(ValueError): tuning.physical(bad)


def test_parameter_mapping_and_neutral_controls():
    assert tuning.physical(tuning.BASE) == pytest.approx(np.ones(6))
    for x in [np.zeros(6),np.ones(6),np.linspace(.1,.9,6)]:
        assert tuning.coordinates(tuning.physical(x)) == pytest.approx(x,abs=1e-15)
    no_direction=tuning.BASE.copy();no_direction[2]=0
    assert tuning.physical(no_direction)[2] == 0


def test_geometry_conversion_keeps_mass_and_speed():
    batch=synthetic_batch();batch.entries[0]['axis_scale']=2
    result=tuning.reported(batch,batch.evaluate(tuning.BASE))
    assert result[:5,0] == pytest.approx(batch.values[:5,0])
    assert result[:5,3] == pytest.approx(batch.values[:5,3])
    assert result[:5,1:3] == pytest.approx(batch.values[:5,1:3]/2)


def test_native_adapter_parity_and_numerical_refinements():
    product=json.loads(tuning.PRODUCT.read_text())
    assert product['checks'] and all(c['passed'] for c in product['checks'])
    assert product['mass_balance_max'] < 1e-10
    assert product['volume_residual_max'] < 1e-10
    assert product['failed_evaluations'] < product['evaluations']


def test_search_respects_bounds_and_saved_training_selection():
    product=json.loads(tuning.PRODUCT.read_text())
    for fit in product['searches']:
        assert tuning.physical(fit['coordinates']) == pytest.approx(list(fit['parameters'].values()))
        assert fit['training_loss'] <= fit['baseline_loss']+1e-8
        for i,value in fit['specification']['fixed'].items():
            assert fit['coordinates'][int(i)] == value
    for r in product['assessments']['shared']['rows']:
        if r['group']=='training':
            assert r['temperature_C'] in tuning.CORE_TEMPERATURES
            assert r['initialization'] in ('sphere','axes')
            assert r['age_min'] in tuning.TRAIN_AGES
