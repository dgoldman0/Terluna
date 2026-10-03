"""Physical invariants of the atmospheric momentum transfer and geographic remap."""
import numpy as np
import pytest
import struct

from climate.gcm.wave_snapshots import edit_namelist, scalar
from climate.waves.weather import apply_weights, drag_coefficient, ocean_weights, stress_equivalent_wind
from climate.waves.weather_checks import compare, read_grid_records


def test_stress_coupling_conserves_signed_momentum_in_both_drag_regimes():
    speed = np.array([0., .1, 4., 7.5, 8., 20., 60.])
    angle = np.linspace(-np.pi, np.pi, len(speed))
    u, v = speed*np.cos(angle), speed*np.sin(angle)
    density = 1.41
    factor = density * drag_coefficient(speed) * speed
    found = stress_equivalent_wind(factor*u, factor*v, density)
    np.testing.assert_allclose(found, [u,v], atol=3e-14, rtol=2e-14)
    # A different wave-solver density preserves the same dimensional stress.
    east, north = stress_equivalent_wind(factor*u, factor*v, 1.2)
    scale = 1.2 * drag_coefficient(np.hypot(east,north)) * np.hypot(east,north)
    np.testing.assert_allclose([scale*east,scale*north], [factor*u,factor*v], atol=1e-14)


def test_ocean_remap_preserves_affine_fields_and_removes_land_stress():
    lon, lat = np.array([0.,2.,4.]), np.array([-2.,0.,2.])
    xx, yy = np.meshgrid(lon,lat)
    x, y = np.array([[1.,3.]]), np.array([[-1.,1.]])
    water = np.ones((3,3),bool)
    indices, weights, fallback = ocean_weights(lon,lat,water,x,y)
    field = 2*xx-3*yy+5
    np.testing.assert_allclose(apply_weights(field,indices,weights),2*x-3*y+5)
    assert not fallback.any()
    water[1,1] = False
    field[:] = 7
    field[1,1] = 1e12
    indices, weights, _ = ocean_weights(lon,lat,water,x,y)
    np.testing.assert_allclose(apply_weights(field,indices,weights),7)
    np.testing.assert_allclose(weights.sum(axis=-1),1)


def test_ocean_remap_fills_empty_stencil_from_closest_ocean():
    lon, lat = np.array([0.,2.,4.,6.]), np.array([-4.,-2.,0.,2.])
    water = np.zeros((4,4),bool)
    water[0,0], water[3,3] = True, True
    indices, weights, fallback = ocean_weights(lon,lat,water,np.array([[3.1]]),np.array([[-.9]]))
    assert fallback.item()
    assert indices[0,0,0] == 15
    field = np.arange(16).reshape(4,4)
    np.testing.assert_array_equal(apply_weights(np.stack([field,field*2]),indices,weights)[:,0,0],[15,30])


def test_surface_coupling_refuses_missing_or_invalid_inputs():
    with pytest.raises(ValueError,match="density"):
        stress_equivalent_wind(1.,2.,0.)
    with pytest.raises(ValueError,match="finite"):
        stress_equivalent_wind(np.nan,0.,1.)
    with pytest.raises(ValueError,match="outside"):
        ocean_weights(np.array([0,1]),np.array([0,1]),np.ones((2,2),bool),np.array([[2]]),np.array([[0]]))


def test_restart_output_edits_leave_physical_settings_intact():
    original = " &plasim_nl\n N_RUN_STEPS = 17160\n NSNAPSHOT = 0\n NSTPS = 480\n PSURF = 121590.\n /END\n"
    edited = edit_namelist(original,dict(N_RUN_STEPS=1446,NSNAPSHOT=1,NSTPS=6))
    assert scalar(edited,"PSURF") == scalar(original,"PSURF")
    assert scalar(edited,"NSTPS") == 6
    with pytest.raises(ValueError,match="Expected one"):
        edit_namelist(original,dict(NHCADENCE=1))
    with pytest.raises(ValueError,match="found 2"):
        edit_namelist(original+" NSTPS = 12\n",dict(NSTPS=6))


def test_independent_raw_reader_preserves_order_sign_and_detects_truncation(tmp_path):
    def record(data):
        marker = struct.pack("<i",len(data))
        return marker+data+marker
    field = np.array([[1,-2,3],[-4,5,-6]],dtype="<f4")
    content = b"".join(record(struct.pack("<8i",180,0,1,0,3,2,t,0))+record((field*t).tobytes()) for t in (5,11))
    path=tmp_path/"snapshot"
    path.write_bytes(content)
    times,values=read_grid_records(path,[180])[180]
    np.testing.assert_array_equal(times,[5,11])
    np.testing.assert_array_equal(values,np.stack((field*5,field*11)))
    path.write_bytes(content[:-1])
    with pytest.raises(ValueError,match="Incomplete"):
        read_grid_records(path,[180])


def test_wave_history_comparison_aligns_the_same_weather_window():
    shorter=np.ones((49,2,9))
    shorter[...,0]=np.arange(49)[:,None]*3600
    shorter[...,3]=2
    shorter[...,5]=10
    longer=np.ones((97,2,9))
    longer[...,0]=np.arange(97)[:,None]*3600
    longer[48:,:,1:]=shorter[...,1:]
    a,b=dict(values=shorter.tolist(),hours=48),dict(values=longer.tolist(),hours=96)
    result=compare(a,b,start_hour=24,offset_b=48)
    assert result["within_five_percent"]
    assert result["paired_wet_samples"]==50
    longer[-1,0,3]=1
    b["values"]=longer.tolist()
    result=compare(a,b,start_hour=24,offset_b=48)
    assert result["hs_relative_max"]==1
    assert not result["within_five_percent"]
