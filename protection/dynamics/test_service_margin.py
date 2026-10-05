import numpy as np

from shared import constants as K
from .service_margin import service_margin


def test_circular_window_margin_includes_perspective_and_polygon_allowance():
    sample=dict(positions={'sun':np.array([0.,0.,K.AU]),'earth':np.array([K.AU,0.,0.])},
        sun_v=np.zeros(3),earth_v=np.zeros(3),sun_a=np.zeros(3),earth_a=np.zeros(3),moon_v=np.zeros(3))
    state=np.array([[0.,0.,15e6,0.,0.,0.]])
    result=service_margin(state,np.eye(3),sample,state[0],[0.,0.],radius=1000.)
    expected=4945.*K.AU/(K.AU-15e6)-1000.
    assert result['covered_fraction']==1.
    assert abs(result['margin_m']+result['circle_sagitta_allowance_m']-expected)<1e-5
    shifted=state.copy();shifted[0,0]=10000.
    missing=service_margin(shifted,np.eye(3),sample,state[0],[0.,0.],radius=1000.)
    assert abs(missing['covered_fraction'])<1e-12
    assert missing['margin_m']<0.
