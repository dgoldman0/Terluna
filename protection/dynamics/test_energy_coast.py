import json
import os
import subprocess
import sys

import numpy as np
import pytest

from .energy_coast import coast_command, schedules
from .attitude_load import rigid_square_load
from .test_pattern_departure import inertial_environment


def test_partial_turn_impulse_scales_with_angle_and_inverse_duration():
    env=inertial_environment();times=np.arange(0.,21600.+1,2.)
    values=[]
    for angle,duration in [(90.,3600.),(30.,14400.)]:
        command=coast_command(env,0.,21600.,angle,3600.,duration)
        value=np.trapezoid(rigid_square_load(command,times)['force_per_mass'],times)
        expected=10000.*np.deg2rad(angle)*1.875/(3*duration)
        assert value==pytest.approx(expected,rel=3e-4)
        values.append(value)
    assert values[1]/values[0]==pytest.approx(1/12,rel=3e-4)
    assert len(schedules())==45
    assert sum(s['angle_deg']==0 for s in schedules())==1


def test_signed_control_can_use_an_inward_impulse_and_keeps_the_norm_cap():
    # HiGHS retains its first process-wide thread configuration. Earlier LP
    # tests use the default; this resource-bounded study uses one thread in
    # fresh processes. Exercise that actual execution context independently.
    script='''
import json
import numpy as np
from protection.dynamics.coast_control import signed_impulse
x,result=signed_impulse(np.eye(3)[None,:,:],np.zeros(3),
    [dict(gradient=[-1.,0.,0.],rhs=.5)],7200.)
print(json.dumps([x.tolist(),result]))
'''
    run=subprocess.run([sys.executable,'-c',script],check=True,capture_output=True,
        text=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    x,result=json.loads(run.stdout)
    assert result['success'], result
    np.testing.assert_allclose(x,[-.5,0.,0.],atol=1e-6)
    assert result['mean_delta_v_m_s']==pytest.approx(.5,abs=1e-6)
    assert result['maximum_acceleration_m_s2']<=.001
