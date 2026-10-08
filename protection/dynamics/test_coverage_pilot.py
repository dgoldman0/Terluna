import numpy as np
import pytest
from shared import constants as K
from .coverage_pilot import (receiver_cells, capacity_rows, cover_capacity,
                            conditional_cost, rotation_step_angle)
from .test_cycling import static_sample


def test_partition_and_finite_square_capacity_conserve_area():
    cells=receiver_cells(4*K.MOON_RADIUS)
    assert abs(sum(c.area for c in cells)/(np.pi*(4*K.MOON_RADIUS)**2)-1)<3e-5
    sample=static_sample()
    # Use an off-axis tile away from cell seams and the target boundary.
    state=np.array([[15e6,1e6,2e6,0,0,0.]])
    normals=np.array([[-1.,0,0]])
    rows,labels=capacity_rows(state,normals,sample,cells,1,1,limb_count=4)
    for source in range(5):
        idx=[i for i,l in enumerate(labels) if l[0]==source]
        total=sum(rows[i,0]*cells[labels[i][1]].area/1e6 for i in idx)
        assert total==pytest.approx(9890.**2*(K.AU/(K.AU-15e6))**2,rel=1e-8)


def test_covering_dual_and_missing_cell_witness():
    a=np.array([[1.,0],[0,2.],[.5,1.]])
    r=cover_capacity(a)
    assert r['success']
    assert r['inventory_tiles']==pytest.approx(1.5e6)
    assert r['dual_lower_bound_million_tiles']==pytest.approx(1.5)
    assert r['dual_max_column_price']<=1+1e-8
    failure=cover_capacity(a,np.array([True,False]))
    assert failure['empty_rows']==[1]


def test_duplicate_geometry_does_not_gain_capacity_per_tile():
    cells=receiver_cells(4*K.MOON_RADIUS)
    state=np.array([[15e6,1e6,2e6,0,0,0.]])
    normal=np.array([[-1.,0,0]])
    one,_=capacity_rows(state,normal,static_sample(),cells,1,1)
    two,_=capacity_rows(np.repeat(state,2,axis=0),np.repeat(normal,2,axis=0),static_sample(),cells,1,2)
    np.testing.assert_allclose(one,two)


def test_cost_mass_closure_and_rotation_units():
    prop=dict(cant_deg=45.,exhaust_velocity_m_s=30000.,efficiency=.7,
              peak_margin_factor=1.25,specific_power_W_kg=300.,propellant_buffer_days=7.)
    low=conditional_cost(1e7,1.,prop,duty=.1)
    high=conditional_cost(1e7,5.,prop,duty=.1)
    assert low['closed'] and high['closed']
    assert high['propulsion_TW']>5*low['propulsion_TW']
    assert low['total_scenario_mass_kg']>low['optical_mass_kg']
    angle=.3
    turn=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    assert rotation_step_angle(np.eye(3),turn)==pytest.approx(angle)
    assert rotation_step_angle(np.eye(3),np.diag([1.,-1.,-1.]))==pytest.approx(0.)
    assert rotation_step_angle(np.eye(3),np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]]))==pytest.approx(0.)


def test_phase_grouping_preserves_uniform_population_capacity():
    cells=receiver_cells(4*K.MOON_RADIUS)
    state=np.array([[15e6,sy*1e6,sz*2e6,0,0,0.] for sy,sz in [(1,1),(1,-1),(-1,1),(-1,-1)]])
    normal=np.tile([-1.,0.,0.],(4,1))
    one,_=capacity_rows(state,normal,static_sample(),cells,1,4)
    grouped,_=capacity_rows(state,normal,static_sample(),cells,2,2)
    np.testing.assert_allclose(one[:,0],grouped.mean(axis=1),atol=1e-12)
