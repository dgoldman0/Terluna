"""Smooth partial-sail coast schedules for an energy-first local search."""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline

from .pattern_departure import sun_frame


def coast_command(env, start, end, angle_deg, delay_s, duration_s, axis=0, sign=-1):
    if not 0 <= angle_deg <= 90 or duration_s <= 0 or delay_s < 0:
        raise ValueError('Invalid partial-sail schedule')
    if start+delay_s+duration_s > end or axis not in (0, 1) or sign not in (-1, 1):
        raise ValueError('Turn must fit inside the tested interval')
    times = np.unique(np.r_[np.arange(start, end, 60.), end,
                            start+delay_s, start+delay_s+duration_s])
    f = np.clip((times-start-delay_s)/duration_s, 0., 1.)
    angle = sign*np.deg2rad(angle_deg)*(10*f**3-15*f**4+6*f**5)
    vector = np.zeros((len(times), 3)); vector[:, axis] = angle
    frames = np.array([sun_frame(env, t) for t in times])@Rotation.from_rotvec(vector).as_matrix()
    return RotationSpline(times, Rotation.from_matrix(frames))


def schedules():
    """Forty-five distinct schedules; zero turn appears only once."""
    result = [dict(angle_deg=0., delay_s=0., duration_s=3600., axis=0, sign=-1)]
    for angle in [15., 30., 60., 90.]:
        for hours in ([1., 4.] if angle == 90. else [1., 2., 4.]):
            for axis in [0, 1]:
                for sign in [-1, 1]:
                    result.append(dict(angle_deg=angle, delay_s=min(7200., (5-hours)*3600.),
                        duration_s=hours*3600., axis=axis, sign=sign))
    return result
