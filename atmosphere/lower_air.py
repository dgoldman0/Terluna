"""Lower-air inputs for mechanics; prescribed thermodynamics, SI units.

The molecular constants match atmosphere/thermal_column.py. A1Profile consumes
existing open-moon-atmospheric-profile/1 exports without importing a renderer.
A1 uses its own documented Rd/Rv closure; it is not silently relabelled N2/O2.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import json
import numpy as np
from numpy.polynomial.legendre import leggauss

R_MOON = 1_737_400.0
GM_MOON = 4.902800118e12
R_UNIVERSAL = 1.380649e-23 * 6.02214076e23
_GL_X, _GL_W = leggauss(16)


def _heights(z, top):
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)) or np.any(z < 0) or np.any(z > top):
        raise ValueError('Height must be finite and inside the supplied profile')
    return z


@dataclass(frozen=True)
class AirConfig:
    pressure_pa: float = 121590.0
    temperature_k: float = 288.0
    oxygen_mole_fraction: float = .175
    vapor_kg_per_kg_dry: float = 0.0
    lapse_k_per_m: float = 0.0
    radius_m: float = R_MOON
    surface_gravity_m_s2: float = GM_MOON / R_MOON**2
    top_m: float = 5000.0

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Finite air parameters required')
        if min(self.pressure_pa, self.radius_m, self.surface_gravity_m_s2,
               self.top_m) <= 0 or not 0 <= self.oxygen_mole_fraction <= 1:
            raise ValueError('Invalid pressure, radius, gravity, top or composition')
        if not 0 <= self.vapor_kg_per_kg_dry <= .1:
            raise ValueError('Vapor mixing ratio outside this model domain')
        ends = [self.temperature_k, self.temperature_k-self.lapse_k_per_m*self.top_m]
        if min(ends) < 150 or max(ends) > 350:
            raise ValueError('Prescribed temperature must remain in 150--350 K')


class LowerAir:
    """Hydrostatic ideal mixture at prescribed T(z), fixed vapor and spherical g."""
    def __init__(self, config=AirConfig()):
        config.validate()
        self.config = config
        molar = .0280134*(1-config.oxygen_mole_fraction)+.031998*config.oxygen_mole_fraction
        rd, rv = R_UNIVERSAL/molar, R_UNIVERSAL/.01801528
        r = config.vapor_kg_per_kg_dry
        self.gas_constant = (rd+r*rv)/(1+r)
        self.metadata = dict(kind='prescribed_hydrostatic_lower_air', **asdict(config))

    def sample(self, z):
        c = self.config
        z = _heights(z, c.top_m)
        t = c.temperature_k-c.lapse_k_per_m*z
        # Gauss quadrature is independent for each query. In the isothermal
        # control the exact spherical potential is used instead.
        if c.lapse_k_per_m == 0:
            integral = c.surface_gravity_m_s2*c.radius_m*z/(c.radius_m+z)/self.gas_constant/c.temperature_k
        else:
            h = z[..., None]*(1+_GL_X)/2
            g = c.surface_gravity_m_s2*(c.radius_m/(c.radius_m+h))**2
            integral = z/2*np.sum(_GL_W*g/(self.gas_constant*(c.temperature_k-c.lapse_k_per_m*h)), axis=-1)
        p = c.pressure_pa*np.exp(-integral)
        rho = p/(self.gas_constant*t)
        g = c.surface_gravity_m_s2*(c.radius_m/(c.radius_m+z))**2
        return dict(z_m=z, pressure_pa=p, temperature_k=t, density_kg_m3=rho,
                    gravity_m_s2=g, vapor_kg_per_kg_dry=np.full_like(z, c.vapor_kg_per_kg_dry))


class A1Profile:
    """Read-only adapter matching A1 log-p, linear-T/r sampling and moist EOS.

    Only thermodynamics is consumed. Optical continuations, selected wind inputs
    and cloud morphology are not climate predictions. No extrapolation is allowed.
    """
    def __init__(self, state, source_sha256=None):
        if state.get('schema') != 'open-moon-atmospheric-profile/1':
            raise ValueError('Expected open-moon-atmospheric-profile/1')
        try:
            rows = state['rows']
            arr = np.array([[r[k] for k in ('z', 'p', 'T', 'r')] for r in rows], float)
            self.radius = float(state['planet']['radius'])
            self.g0 = float(state['planet']['g0'])
            self.top = float(state['upper']['top_m'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('Incomplete A1 thermodynamic state') from exc
        if arr.ndim != 2 or arr.shape[0] < 2 or arr.shape[1] != 4 or not np.all(np.isfinite(arr)):
            raise ValueError('At least two finite A1 rows required')
        if (not np.all(np.diff(arr[:, 0]) > 0) or not np.all(np.diff(arr[:, 1]) < 0)
                or np.any(arr[:, 1:3] <= 0) or np.any(arr[:, 3] < 0)
                or arr[0, 0] != 0 or arr[-1, 0] != self.top
                or not all(np.isfinite(x) and x > 0 for x in (self.radius, self.g0, self.top))):
            raise ValueError('Invalid A1 profile domain or thermodynamics')
        self.rows = arr.copy()
        self.rows.flags.writeable = False
        self.metadata = dict(kind='borrowed_A1_thermodynamics', source_sha256=source_sha256,
                             world=state.get('world'), regime=state.get('regime'),
                             dry_gas_constant_j_kg_k=287.05, vapor_gas_constant_j_kg_k=461.5,
                             extrapolation=False)

    @classmethod
    def from_file(cls, path):
        data = Path(path).read_bytes()  # A missing export is an error; no fallback.
        return cls(json.loads(data), hashlib.sha256(data).hexdigest())

    def sample(self, z):
        z = _heights(z, self.top)
        p = np.exp(np.interp(z, self.rows[:, 0], np.log(self.rows[:, 1])))
        t = np.interp(z, self.rows[:, 0], self.rows[:, 2])
        r = np.interp(z, self.rows[:, 0], self.rows[:, 3])
        tv = t*(1+r/(287.05/461.5))/(1+r)
        return dict(z_m=z, pressure_pa=p, temperature_k=t,
                    density_kg_m3=p/(287.05*tv),
                    gravity_m_s2=self.g0*(self.radius/(self.radius+z))**2,
                    vapor_kg_per_kg_dry=r)
