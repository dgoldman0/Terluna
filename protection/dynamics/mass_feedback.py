"""Distributed additional hardware mass on the existing parallel-square model."""
import numpy as np

from .parallel_shadow import ParallelPattern


def allocated_mass_ratio(optical_specific_peak_W_kg,prop,fixed_mass_ratio=1.2):
    """Scalar power feedback for each member, under proportional-load sizing.

    This allocation must be checked against actual loads after propagation.
    Fixed mass and power mass share the optical square's uniform distribution.
    """
    fraction=np.asarray(optical_specific_peak_W_kg)*prop['peak_margin_factor']/prop['specific_power_W_kg']
    if fixed_mass_ratio<1 or np.any(fraction>=1) or np.any(fraction<0):
        raise ValueError('Unclosed mass allocation')
    return fixed_mass_ratio/(1-fraction)


class MassiveParallelPattern(ParallelPattern):
    """Unchanged optical area; gravity per mass unchanged, sail force diluted."""
    def __init__(self,env,command,mass_ratio,**kwargs):
        super().__init__(env,command,**kwargs)
        self.mass_ratio=np.asarray(mass_ratio)
        if np.any(self.mass_ratio<1):raise ValueError('Mass cannot remove baseline optics')

    def acceleration(self,t,state,diagnostics=False):
        total,light,sail=super().acceleration(t,state,diagnostics=True)
        changed=sail/np.broadcast_to(self.mass_ratio,(len(state),))[:,None]
        force=total-sail+changed
        if diagnostics:return force,light,changed
        return force
