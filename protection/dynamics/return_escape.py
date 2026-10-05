"""Smooth in-plane return freedoms, preserving the existing optical normal.

Changing the square edges changes shadows and finite-area forces. These
commands require coupled propagation even though their normal is unchanged.
"""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline


def quintic_fraction(t, start, duration):
    if duration <= 0:
        raise ValueError('Positive transition duration required')
    x = np.clip((np.asarray(t)-start)/duration, 0., 1.)
    return x**3*(10.-15.*x+6.*x*x)


def rolled_return(base, start, end, angle_deg=0., turn_s=5400.,
                  recovery_start_s=64800., recovery_duration_s=14400.):
    """Quintic roll, hold and recovery; shared normal and bounded smooth edges."""
    if not angle_deg:
        return base
    if start+turn_s > recovery_start_s or recovery_start_s+recovery_duration_s > end:
        raise ValueError('Roll and recovery must fit without overlap')
    t = np.unique(np.r_[np.arange(start-240., end+241., 60.), start, start+turn_s,
                        recovery_start_s, recovery_start_s+recovery_duration_s, end])
    angle = np.deg2rad(angle_deg)*(quintic_fraction(t, start, turn_s)-
        quintic_fraction(t, recovery_start_s, recovery_duration_s))
    spin = Rotation.from_rotvec(np.c_[np.zeros((len(t), 2)), angle])
    return RotationSpline(t, base(t)*spin)


def applied_impulse(controls, windows, start, end):
    """Per-member sum of burn-vector magnitudes over an executed interval."""
    widths = np.diff(windows, axis=1).ravel()
    if np.any(widths <= 0) or end < start:
        raise ValueError('Invalid burn or integration interval')
    def primitive(t):
        x = np.clip((t-windows[:, 0])/widths, 0., 1.)
        return x-np.sin(2*np.pi*x)/(2*np.pi)
    return np.linalg.norm(controls, axis=-1)@(primitive(end)-primitive(start))
