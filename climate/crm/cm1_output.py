"""How CM1 r22.0 writes what the ring and box analyses read: its constants, its temperature and its radiative
tendencies, which under the Bryan-Fritsch equation set (eqtset = 2, every lunar case's) it stores at constant volume.

Just after the RRTMG calls, radiation_driver.F turns RRTMG's heating rate hr, the rate at constant pressure that every
radiation code gives, into CM1's temperature tendency, lwten = hr cp/cvm and swten likewise, with
cvm = cv + cvv qv + cpl (qc + qr) + cpi (qi + qs + qg). The potential temperature then takes (swten + lwten) cv/cp/pi.
So lwten and swten as written run cp/cvm, about 1.38 in the lunar air, above RRTMG's own rates, while the air takes
RRTMG's heating. rrtmg_rate() gives back RRTMG's own rate, for comparisons with other radiation codes; applied_rate()
gives the temperature tendency at constant pressure that the potential temperature takes (pi times its tendency), as
budgets count CM1's other potential-temperature tendencies. CM1 hands RRTMG its own cp and, since cm1r21.1, the
ground's skin temperature as the air's temperature at the ground (dum2(i,j,1) = tsk); offline columns built to match
CM1's call take both. Under eqtset = 1 the tendencies are stored as RRTMG gives them.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# CM1 r22.0's own constants (constants.F, its defaults), J/kg/K; its cp is also the one it hands RRTMG
RD, CP, CPV, RV, CPL, CPI = 287.04, 1005.7, 1870.0, 461.5, 4190.0, 2106.0
CV, CVV = CP - RD, CPV - RV
P00 = 1.0e5


def eqtset(case: Path) -> int:
    """A case's equation set, from its namelist.input: 1, the traditional set; 2, Bryan and Fritsch's."""
    m = re.search(r'^\s*eqtset\s*=\s*(\d+)', (Path(case) / 'namelist.input').read_text(), re.M)
    if not m:
        raise KeyError('eqtset')
    return int(m.group(1))


def temperature(d):
    """Temperature (K) from a snapshot's potential temperature and pressure, through CM1's own Exner function."""
    return d['th'] * (d['prs'] / P00) ** (RD / CP)


def cp_over_cvm(d):
    """The factor CM1 multiplies RRTMG's heating by under eqtset = 2, from a snapshot's own water (kg/kg), as
    radiation_driver.F reckons it; a mixing ratio the snapshot lacks counts as none."""
    none = np.zeros_like(np.asarray(d['qv'], float))

    def q(*names):
        return np.maximum(0.0, sum((np.asarray(d[k], float) for k in names if k in d), none))
    return CP / (CV + CVV * q('qv') + CPL * q('qc', 'qr') + CPI * q('qi', 'qs', 'qg'))


def rrtmg_rate(tendency, d, equation_set: int):
    """RRTMG's own heating rate from a stored lwten or swten (in the tendency's units): the tendency times cvm/cp
    under eqtset = 2, as stored under eqtset = 1."""
    return np.asarray(tendency) / cp_over_cvm(d) if equation_set == 2 else np.asarray(tendency)


def applied_rate(tendency, equation_set: int):
    """The temperature tendency at constant pressure that CM1's potential temperature takes from a stored lwten or
    swten, pi times its tendency: the tendency times cv/cp under eqtset = 2, as stored under eqtset = 1."""
    return np.asarray(tendency) * (CV / CP) if equation_set == 2 else np.asarray(tendency)
