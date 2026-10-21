"""Meaningful processing controls separate from the frozen trajectories."""
import numpy as np
import pytest

from atmosphere.snow.welss_archive_comparison import mean_lattice,profile_metrics


def test_mean_lattice_distinguishes_31_from_four_without_claiming_uniqueness():
    # Means of31 three-digit, order-1000 samples, rounded to archive float32.
    sums=np.array([12345,13211,14217,15231,16127])*10.
    values=(sums/31).astype(np.float32)
    assert mean_lattice(values,31)['compatible']==len(values)
    assert mean_lattice(values,62)['compatible']==len(values)
    assert mean_lattice(values,4)['compatible']<len(values)


def test_regional_metrics_keep_nucleation_discrepancy_visible():
    z=np.array([0.,20.,2700.,2980.,3000.])
    reference=np.array([2.,2.,1.,1.,np.nan])
    model=np.array([2.,2.,1.5,1.5,0.])
    result=profile_metrics(model,reference,z)
    assert result['all_populated']['relative_L1']==pytest.approx(1/6)
    assert result['below_nucleation']['relative_L1']==0.
    assert result['nucleation_belt']['relative_L1']==.5
    assert result['nucleation_belt']['maximum_absolute_relative_difference']==.5
