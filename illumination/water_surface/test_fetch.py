import numpy as np
import pytest
from illumination.water_surface.fetch import inverse_wave_age


def test_fetch_closure_scale_invariance_and_fully_developed_limit():
    # Keeping g*F/U**2 unchanged must keep wave age unchanged.
    age = inverse_wave_age(16000., 4., 1.62)
    assert age == pytest.approx(inverse_wave_age(96000., 8., 1.08))
    assert inverse_wave_age(1e10, 4., 1.62) == pytest.approx(.84)
    ages = [inverse_wave_age(f, 4., 1.62) for f in np.geomspace(1e3, 1e8, 30)]
    assert np.all(np.diff(ages) <= 0)
    assert .84 < age < 5


@pytest.mark.parametrize('args', [(0,4,1.62), (10,-1,1.62), (10,4,np.nan)])
def test_invalid_fetch_inputs_fail(args):
    with pytest.raises(ValueError):
        inverse_wave_age(*args)
