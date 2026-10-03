"""Sampling and mass conservation for the saved regional cloud survey."""
import numpy as np
import pytest
from climate.crm.evening_columns import snapshot_numbers, column_summary
from shared.constants import SYNODIC_MONTH_DAYS


def test_second_complete_lunar_cycle_excludes_spinup():
    record={'configuration':{'output_s':10800.,'days':59.}}
    numbers=snapshot_numbers(record)
    days=(numbers-1)*10800/86400
    assert len(days)==236
    assert days[0]==29.625
    assert days[-1]==59.
    assert days.min()>=SYNODIC_MONTH_DAYS
    assert days.max()<=2*SYNODIC_MONTH_DAYS


def test_cloud_paths_partition_partial_layers_and_respect_terrain():
    edges=np.array([[500.,0.],[1500.,1000.],[2500.,2000.],[40500.,40000.]])
    centres=(edges[:-1]+edges[1:])/2
    fields={s:np.full((3,2),2e-5 if s=='qc' else 0.) for s in ('qc','qr','qi','qs','qg')}
    result=column_summary(fields,np.full((3,2),.5),centres,edges)
    np.testing.assert_allclose(result['qc_path_kg_m2'],.4)
    np.testing.assert_allclose(result['cloud_path_above_20000_kg_m2'],[.205,.2])
    np.testing.assert_allclose(result['cloud_top_m'],centres[-1])
    assert result['low_cloud'].all()
    edges[1]=edges[0]
    with pytest.raises(ValueError,match='thickness'):
        column_summary(fields,np.ones((3,2)),centres,edges)
