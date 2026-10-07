"""Independent geometry and extended-source checks for the calendar reader."""
import json
from pathlib import Path

import numpy as np
import pytest

from geography.lunar_ephemeris import earth_and_sun
from illumination.calendar.disk import samples
from illumination.calendar.geometry import geometry, local_frame
from shared.constants import MOON_RADIUS


def test_topocentric_earth_parallax_and_local_frame():
    jd=2465476.3704272127
    g=earth_and_sun(np.array(jd))
    frame=local_frame(90.,0.)
    np.testing.assert_allclose(frame@frame.T,np.eye(3),atol=1e-14)
    centre=np.degrees(np.arcsin(g['earth']@frame[2]))
    top=geometry(jd,90.,0.)
    # The observer is closer to the limb than the centre, depressing Earth's elevation.
    assert .24 < centre-top['earth_elevation_deg'] < .29
    distance=np.linalg.norm(g['earth']*g['earth_distance_m']-MOON_RADIUS*frame[2])
    assert top['earth_distance_m']==pytest.approx(distance)


def test_vector_dates_and_scalar_dates_agree():
    dates=np.array([2451544.5,2465476.37,2634531.47])
    values=geometry(dates,-73.875,25.125)
    for i,date in enumerate(dates):
        for key,value in geometry(date,-73.875,25.125).items():
            assert values[key][i]==pytest.approx(value,rel=2e-15,abs=1e-9)


def test_uniform_solar_disk_horizontal_flux_in_vacuum():
    # The horizontal flux of a symmetric small disk is sin(elevation) times
    # its mean direction cosine. It is positive at the geometric horizon.
    for elevation in [0.,30.,90.]:
        rays,weights=samples(elevation,.267,order=32)
        assert weights.sum()==pytest.approx(1.)
        actual=weights@np.maximum(0.,np.sin(np.radians(rays)))
        if elevation==0:
            expected=2*np.radians(.267)/(3*np.pi)
            assert actual==pytest.approx(expected,rel=.002)
        else:
            assert actual==pytest.approx(np.sin(np.radians(elevation)),rel=1e-5)


def test_earth_phase_orients_the_illuminated_half_without_double_phase_factor():
    full,w=samples(0.,1.,0.,order=32)
    upper,wu=samples(0.,1.,90.,0.,order=32)
    lower,wl=samples(0.,1.,90.,np.pi,order=32)
    assert np.dot(w,full)==pytest.approx(0.,abs=1e-12)
    assert np.dot(wu,upper)>0.5
    assert np.dot(wl,lower)<-0.5
    assert wu.sum()==pytest.approx(1.)
    assert wl.sum()==pytest.approx(1.)
    assert samples(0.,1.,180.)[1].sum()==0.


def test_point_source_grazing_and_zenith_flux_in_vacuum():
    pytest.importorskip('numba')
    from illumination.calendar.transport import point_table
    edges=np.array([1000.,1100.])
    angles=np.radians([-1.,0.,30.,90.])
    value=point_table(edges[:1],angles,edges,np.zeros((1,1)))
    np.testing.assert_allclose(value[0,:,0,1],[0.,0.,.5,1.],atol=1e-15)
    np.testing.assert_allclose(value[0,:,0,0],[0.,1.,1.,1.],atol=1e-15)


def test_committed_astronomy_geometry_matches_domain():
    product=json.loads((Path(__file__).parent/'results/astronomy.json').read_text())
    for case in product['golden_geometry']:
        values=geometry(case['jd_tt'],case['longitude'],case['latitude'])
        for key,value in values.items():
            assert case['values'][key]==pytest.approx(value,rel=1e-12,abs=1e-10)
