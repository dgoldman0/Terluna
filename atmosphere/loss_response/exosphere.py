"""What the sunlit exosphere beyond the shield's shadow loses, with and without a lunar magnetosphere.

    python -m atmosphere.loss_response.exosphere   # after model.py; writes results/exosphere_loss.json

absorption.py counts the ions sunlight makes in the exosphere outside the
shield's shadow cylinder. This step follows them, adds what else sunlight and
the solar wind do to the same molecules, and asks how much leaves the Moon. It
is a screening step: each process is a physical bound or a comparison of
timescales with stated ranges, and no plasma is simulated.

The exosphere. Chamberlain's ballistic density for N2 and O2 separately, from
the thermal column's exobase density, temperature and radius (absorption.py:
no satellite orbits, two-body, counted to a third of the Hill radius). The
column keeps the air well mixed to the exobase, 17.5% O2. Above a homopause,
molecular diffusion separates the gases; with the homopause at 0.01 Pa, where
molecular diffusion of O2 in N2 overtakes eddy mixing of about 100 m^2/s, O2
falls to a few percent at the exobase. Both are carried.

Sunlight. Ionization at the repository's rate (absorption.ionization_rate).
Dissociation of N2 through its predissociating absorption at 80-100 nm
(2e-21 m^2, the extreme-ultraviolet band's cross-section) and of O2 in the
Schumann-Runge continuum (3e-22 m^2, flux-weighted) and at Lyman-alpha, on
the WHI 2008 spectrum: 4.4e-7 and 2.3e-6 per second at quiet Sun, against
4.1e-7 and 2.4e-6 from Heays et al. (2017, A&A 602, A105; their solar field
times 37,700 for 1 AU). Solar maximum raises N2's rate by the loss response's
ultraviolet factor and O2's by 1.3 (assumed: light above 130 nm varies by tens
of percent over the cycle).

Fragments. A dissociated molecule's atoms leave with 0.1-1.7 eV each (N) and
mostly 0.2-1.5 eV (O), above the 0.1-0.2 eV an atom needs to escape from two to
five lunar radii. An atom leaves unless its path meets the exobase sphere,
where the collisional air stops it. For an atom of speed v launched
isotropically at radius r above an exobase at r_c the escaping share is
(1 + sqrt(1 - s^2)) / 2, with s^2 = (r_c/r)^2 (1 + 2GM(1/r_c - 1/r)/v^2) for
paths bent by gravity (0.8 eV for N, 0.5 eV for O). The share of atoms
energetic enough to escape is 0.75-1 for O (the continuum's weak long-wave
end, near 175 nm, gives slow atoms) and 0.9-1 for N.

Charge exchange. Solar-wind protons (5 cm^-3 at 400 km/s) take an electron
from exospheric molecules with a cross-section of 1.5e-19 m^2, the size
measured for protons on N2 and O2 near 1 keV (Lindsay and Stebbings 2005
compile the measurements; the value was not re-read for this step). Along
straight paths parallel to the Sun-Moon line, down to the exobase, each
exchange puts an ion into the flowing wind, inside the optical shadow as well:
the swarm stops the wind that meets it, but the wake behind it refills over
about eight shield radii, within the 78,000 km from the September screen to the
Moon. A loaded wind is deflected before it reaches the dense exosphere, so this
is an upper bound; the range takes 0.1-1 of it, times the share of the orbit
outside Earth's magnetotail. The ring fleet's screen, about 20,000 km out, lies
under three screen radii from the Moon, so the Moon may sit in its unrefilled
wake. A sensitivity takes the wake to close linearly over eight screen radii
and the wind to be absent on paths within its core at the Moon (2.6 lunar radii
behind a 4-radius screen): the "wake" rows, which rest on that refill length
and are no plasma calculation.

No magnetosphere. In the flowing wind the motional electric field (2 mV/m)
pulls an N2+ ion about 40,000 times harder than lunar gravity at three lunar
radii; gravity holds an ion only where the flow slows to a few metres per
second, at an obstacle's surface. The exosphere's ions are therefore carried
off, as at comets, where ion production likewise far exceeds the solar wind's
mass flux through the region that makes the ions: the loaded wind slows and
its interaction region grows (here to tens of lunar radii) instead of refusing
the ions. This replaces the loss response's mass-loading cap. The pickup ions
that curve back into the air, 0.1-0.5 of them, stay but sputter 1-10 molecules
each (the loss response's ranges; the high yields belong to keV ions, and mass
loading lowers the energies near the Moon).

A lunar magnetosphere. A dipole of moment M along the spin axis holds the wind
off at r_mp = (mu0 M^2 / (8 pi^2 rho v^2))^(1/6), pressure balance with the
field doubled at the boundary: 10.0 lunar radii for the September design's
1.5e21 A m^2 in the typical wind (its report's vacuum-pressure diagnostic gives
7.2 at 2 nPa). The boundary has the form r_mp (2 / (1 + cos theta))^0.6
(Shue et al. 1997), open down the tail. Ions made outside it, and molecules
charge-exchanged there, are carried off. Inside, a point whose dipole shell
L = r / cos^2(latitude) exceeds 1-2.5 r_mp lies on an open field line into the
tail, and its ions leave (Earth's polar-cap boundary, at 70-80 degrees
magnetic latitude, maps to 0.85-3.3 times its stand-off). On closed field lines
three rates compete. Draining along the field into the collisional air, where
the ions recombine and stay, at kappa sqrt(2kT/m + GM(1/r_c - 1/r)), kappa
0.5-1 for the ambipolar field and the plasma's own pressure. Convection driven
by the solar wind, which carries closed flux to the dayside boundary to be
opened, at eps E_sw / B over a path of r_mp: eps is 0.05-0.15, the share of the
wind's electric field that drives convection (about a tenth at Earth;
Mercury's dayside reconnection rate is 0.15, three times Earth's, DiBraccio et
al. 2013), and the Moon turns too slowly to hold a corotating plasmasphere.
And dissociative recombination inside the magnetosphere (N2+: 2.2e-13
(300/T_e)^0.39 m^3/s, T_e 1000 K), whose fast atoms (0.5-3 eV, 1.5 taken)
escape unless aimed at the exobase. The ion density is the local balance of
production and these losses, with the rates for N2+ applied to both ions. The
fate of the ions a magnetosphere keeps from the wind is thus a comparison of
timescales, not a plasma model: entry at the cusps and the sputtering it
brings (taken as zero here), the magnetosphere's own plasma pressure, the
passages through Earth's magnetotail and the extra escape a weak field can
drive (Egan et al. 2019; Gunell et al. 2018) are left out.

Outputs. The largest UV transmission each budget allows once the sunlit
exosphere's central loss is added to the loss response's ultraviolet-driven
loss, for protected radii of 3-10 lunar radii, with no magnetosphere (and with
the ring fleet's wake) and with the September moment. The ultraviolet-driven
state at each radius is the limb tracing's (traced.py): the UV transmission
passes gaps over that radius's aperture and unfiltered light falls beyond it.
A radius counts only where it covers the thermosphere's heating as the design
study asks: light beyond the aperture gives at most a tenth of the heat. And,
at the loss response's own allowed transmissions with the ring fleet's 4-radius
shadow, the ions, fragments and charge exchange, the losses low, central and
high, and the shadow radius at which the loss falls to a tenth of the budget.
assess() also gives the loss against the shadow radius and against the dipole
moment (1e19-1e23 A m^2) for any column.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from shared.constants import AVOGADRO, BOLTZMANN, ELEMENTARY_CHARGE, MOON_GM, MOON_RADIUS, VACUUM_PERMEABILITY
from atmosphere.thermal_column import solve_column
from atmosphere.middle_atmosphere import escape
from atmosphere.loss_response import absorption as ab
from atmosphere.loss_response import model as lr
from atmosphere.loss_response import traced

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.exosphere-loss/2'
LEVELS = ('low', 'central', 'high')

MOLECULE_KG = dict(N2=0.0280134 / AVOGADRO, O2=0.0319988 / AVOGADRO)
ATOM_KG = dict(N2=0.0140067 / AVOGADRO, O2=0.0159994 / AVOGADRO)
HOMOPAUSE_PA = 0.01
OUTER_R = ab.OUTER['third_hill']

N2_PREDISSOCIATION_M2 = 2e-21                    # 80-100 nm
O2_CONTINUUM_M2 = 3e-22                          # Schumann-Runge continuum, 122.5-175 nm, flux-weighted
O2_LYMAN_ALPHA_M2 = 1e-24
FAR_UV_ACTIVITY = dict(quiet=1.0, solar_maximum_stress=1.3)
HEAYS_1AU_S = dict(N2=1.1e-11 * 37700, O2=6.4e-11 * 37700)   # Heays et al. 2017, solar field scaled to 1 AU

FRAGMENT_EV = dict(N2=0.8, O2=0.5)
RECOMBINATION_EV = 1.5
ENERGETIC_SHARE = dict(N2=(0.9, 1.0), O2=(0.75, 1.0))

CHARGE_EXCHANGE_M2 = 1.5e-19
CHARGE_EXCHANGE_SHARE = (0.1, 1.0)
IMF_T = 5e-9

SEPTEMBER_MOMENT = 1.5e21
FLARING = 0.6
MAGNETOSPHERE = dict(open_over_standoff=(2.5, 1.0),      # (low loss, high loss)
                     convection_efficiency=(0.05, 0.15),
                     drain_factor=(1.0, 0.5))
ELECTRON_K = 1000.0
MOMENTS = tuple(float(m) for m in np.geomspace(1e19, 1e23, 17))
RADII = tuple(float(x) for x in np.round(np.arange(1.5, 11.51, 0.25), 2))
SHADOWS_R = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0)
RING_FLEET_DISTANCE_KM = 20000.0                 # the ring fleet's screen (decisions.md, 2026-10-07)
WAKE_REFILL_SCREEN_RADII = 8.0                   # the wind closes in behind an absorbing screen over about this many radii
MAGNETOSPHERE_KEYS = ('none', 'september', 'none_wake')

EVIDENCE = ('Screening step. Chamberlain exospheres for N2 and O2 from the thermal column\'s exobase; photoionization '
            'and photodissociation rates on the WHI 2008 spectrum, the latter checked against Heays et al. (2017); '
            'fragment escape from kinematics; charge exchange with the solar wind along straight paths (an upper '
            'bound); ion fates with no magnetosphere from the dominance of the motional electric field over gravity, '
            'and with a lunar dipole from pressure balance, dipole field-line geometry and a comparison of draining, '
            'convection and recombination timescales. No plasma is simulated.')
READING_RULE = ('Loss rates are kg/s of atmosphere. allowed_summary[magnetosphere][shield][protected radius][budget] '
                'is the largest UV transmission through gaps whose ultraviolet-driven loss at the traced state for '
                'that radius plus the sunlit exosphere\'s central loss stays within the budget, as a range over the '
                'upper-air treatments, quiet Sun and solar maximum, with the number of those cases that cannot meet '
                'it; stress_summary[case] gives the same at each stress case alone (cycle 19\'s year, the strongest '
                'on record, and the escape model\'s 2.5 on the ultraviolet); '
                'magnetosphere is none, none_wake (no magnetosphere, the wind absent within the ring fleet\'s '
                'unrefilled wake, a sensitivity) or september (1.5e21 A m^2). These replace the loss response\'s '
                'allowed transmissions, its capped pickup and its ion sputtering. cases evaluates the loss response\'s '
                'own allowed transmissions with the ring fleet\'s 4-radius shadow: sunlit carries the ions made, the '
                'fragments escaping and the charge-exchange upper bound, with and without the wake; '
                'no_magnetosphere, no_magnetosphere_wake and september_magnetosphere the low, central and high loss; '
                'heating_radius_R the smallest protected radius that covers the heating; '
                'radius_for_a_tenth_of_budget_R the smallest shadow radius at which the central loss is a tenth of '
                'the budget (null if none up to 11.5 lunar radii).')


def pick(pair, level):
    """A range's low-loss end, geometric middle or high-loss end."""
    if level == 'low':
        return pair[0]
    if level == 'high':
        return pair[1]
    return math.sqrt(pair[0] * pair[1])


def photon_flux(lo_nm, hi_nm, measured=None):
    """WHI 2008 photons per m^2 per s between two wavelengths; with a measured activity, each band weighted by its
    own factor (cycle.py)."""
    w, f = escape.whi_quiet_sun()
    sel = (w >= lo_nm) & (w < hi_nm)
    photons = f[sel] * 0.1 * w[sel] * 1e-9 / (lr.PLANCK * lr.LIGHT)
    if measured is not None:
        photons = photons * traced.band_factors(w[sel], measured)
    return float(photons.sum())


def rates(activity='quiet'):
    """Per-molecule ionization and dissociation rates (s^-1) in full sunlight at 1 AU: for a named activity the escape
    model's factor on the ultraviolet up to 100 nm and FAR_UV_ACTIVITY above, for a measured one each band's own."""
    ion = ab.ionization_rate(activity)
    act = lr.activity_of(activity)
    if 'bands' in act:
        n2 = N2_PREDISSOCIATION_M2 * photon_flux(80.0, 100.0, act)
        o2 = O2_CONTINUUM_M2 * photon_flux(122.5, 175.0, act) + O2_LYMAN_ALPHA_M2 * photon_flux(121.0, 122.5, act)
    else:
        n2 = N2_PREDISSOCIATION_M2 * photon_flux(80.0, 100.0) * act['uv']
        o2 = (O2_CONTINUUM_M2 * photon_flux(122.5, 175.0)
              + O2_LYMAN_ALPHA_M2 * photon_flux(121.0, 122.5)) * FAR_UV_ACTIVITY[activity]
    return dict(N2=dict(ionization=ion, dissociation=n2), O2=dict(ionization=ion, dissociation=o2))


def exobase_shares(profile, homopause_pa=None):
    """Mole fractions at the exobase: the column's well-mixed air, or diffusive separation above a homopause."""
    x = ab.O2_SHARE
    p_exo = float(profile['pressure_Pa'][-1])
    if homopause_pa is not None and p_exo < homopause_pa:
        ratio = x / (1 - x) * (p_exo / homopause_pa) ** ((MOLECULE_KG['O2'] - MOLECULE_KG['N2']) / ab.AIR_KG)
        x = ratio / (1 + ratio)
    return dict(N2=1.0 - x, O2=x)


def escaping_share(r_R, r_c, energy_ev, atom_kg):
    """Share of atoms launched isotropically at radius r (lunar radii) that escape without meeting the exobase sphere."""
    r = np.asarray(r_R, dtype=float)
    v2 = 2 * energy_ev * ELEMENTARY_CHARGE / atom_kg
    gm = MOON_GM / MOON_RADIUS
    s2 = (r_c / r) ** 2 * (1 + 2 * gm / v2 * (1 / r_c - 1 / r))
    share = 0.5 * (1 + np.sqrt(np.clip(1 - s2, 0.0, None)))
    return np.where(v2 < 2 * gm / r, 0.0, share)


def sunlit(profile, activity, shadows_R, homopause_pa=None, outer_R=OUTER_R, points=6000):
    """Ions made and fragments escaping (kg/s) per species in the ballistic exosphere outside shadows of radii shadows_R."""
    n_c, t_c, r_c = ab.exobase_state(profile)
    r = np.geomspace(r_c, outer_R, points)
    shares = exobase_shares(profile, homopause_pa)
    j = rates(activity)
    shadows = np.atleast_1d(np.asarray(shadows_R, dtype=float))
    lit = np.sqrt(np.clip(1.0 - (shadows[:, None] / r[None, :]) ** 2, 0.0, None))
    out = {}
    for s, kg in MOLECULE_KG.items():
        bal, _ = ab.exosphere_density(r, shares[s] * n_c, t_c, r_c, mass_kg=kg)
        mass = kg * bal * 4 * math.pi * r ** 2 * MOON_RADIUS ** 3
        esc = escaping_share(r, r_c, FRAGMENT_EV[s], ATOM_KG[s])
        out[s] = dict(ions=j[s]['ionization'] * np.trapezoid(lit * mass, r, axis=1),
                      fragments=j[s]['dissociation'] * np.trapezoid(lit * mass * esc, r, axis=1))
    return out


def wake_core_R(protected_R, distance_km=RING_FLEET_DISTANCE_KM, refill=WAKE_REFILL_SCREEN_RADII):
    """Radius (lunar radii) of the screen's wake that the solar wind has not refilled at the Moon: the aperture of the
    protected radius at the screen's distance, closing linearly over `refill` aperture radii."""
    screen = protected_R + ab.SUN_ANGULAR_RADIUS * distance_km * 1e3 / MOON_RADIUS
    return max(0.0, screen * (1.0 - distance_km * 1e3 / MOON_RADIUS / (refill * screen)))


def charge_exchange(profile, homopause_pa=None, outer_R=OUTER_R, b_points=500, z_points=600, wind=lr.SOLAR_WIND,
                    sheltered_R=0.0):
    """Ions (kg/s) the solar wind makes by charge exchange with ballistic molecules on straight paths, at full exposure;
    paths within sheltered_R (lunar radii) of the Sun-Moon line, in a screen's unrefilled wake, carry no wind."""
    n_c, t_c, r_c = ab.exobase_state(profile)
    shares = exobase_shares(profile, homopause_pa)
    table = np.geomspace(r_c, outer_R, 4000)
    density = {}
    for s, kg in MOLECULE_KG.items():
        bal, esc = ab.exosphere_density(table, shares[s] * n_c, t_c, r_c, mass_kg=kg)
        density[s] = (np.log(np.maximum(bal, 1e-300)), np.log(np.maximum(bal + esc, 1e-300)))
    b = np.unique(np.r_[np.linspace(0.0, r_c, b_points // 2, endpoint=False),
                        r_c + (outer_R - r_c) * np.linspace(0.0, 1.0, b_points - b_points // 2) ** 2])
    z0 = np.sqrt(np.clip(r_c ** 2 - b ** 2, 0.0, None))
    z1 = np.sqrt(np.clip(outer_R ** 2 - b ** 2, 0.0, None))
    u = np.linspace(0.0, 1.0, z_points) ** 2
    z = z0[:, None] + (z1 - z0)[:, None] * u[None, :]
    r = np.clip(np.sqrt(b[:, None] ** 2 + z ** 2), r_c, outer_R)
    sides = np.where(b < r_c, 1.0, 2.0)
    total = np.zeros_like(b)
    ballistic_mass = np.zeros_like(b)
    for s, kg in MOLECULE_KG.items():
        ln_bal, ln_all = density[s]
        col_bal = sides * np.trapezoid(np.exp(np.interp(r, table, ln_bal)), z, axis=1) * MOON_RADIUS
        col_all = sides * np.trapezoid(np.exp(np.interp(r, table, ln_all)), z, axis=1) * MOON_RADIUS
        total += col_all
        ballistic_mass += col_bal * kg
    flux = wind['density_m3'] * wind['speed_m_s']
    exchanged = -np.expm1(-CHARGE_EXCHANGE_M2 * total)
    per_path = np.where((total > 0) & (b >= sheltered_R), exchanged * ballistic_mass / np.maximum(total, 1e-300), 0.0)
    return float(flux * np.trapezoid(per_path * 2 * math.pi * b, b) * MOON_RADIUS ** 2)


def standoff_R(moment, wind=lr.SOLAR_WIND):
    """Magnetopause stand-off (lunar radii) of a dipole of the given moment (A m^2), field doubled at the boundary."""
    pressure = wind['density_m3'] * lr.PROTON_KG * wind['speed_m_s'] ** 2
    return (VACUUM_PERMEABILITY * moment ** 2 / (8 * math.pi ** 2 * pressure)) ** (1 / 6) / MOON_RADIUS


def _field_path(u):
    """Integral of sqrt(1 + 3 u^2) du: a dipole field line's length per L from the equator to latitude asin(u)."""
    return 0.5 * u * np.sqrt(1 + 3 * u ** 2) + np.arcsinh(math.sqrt(3) * u) / (2 * math.sqrt(3))


class Exosphere:
    """The sunlit ballistic exosphere of one column on a spherical grid, for the magnetosphere's partition."""

    def __init__(self, profile, activity, homopause_pa=None, outer_R=OUTER_R, r_points=220, mu_points=240, phi_points=12):
        self.n_c, self.t_c, self.r_c = ab.exobase_state(profile)
        self.activity = activity
        shares = exobase_shares(profile, homopause_pa)
        j = rates(activity)
        self.r = np.geomspace(self.r_c, outer_R, r_points)
        mu = (np.arange(mu_points) + 0.5) / mu_points * 2 - 1
        phi = (np.arange(phi_points) + 0.5) / phi_points * (math.pi / 2)   # a quarter; the rest by symmetry
        self.mu, phi = np.meshgrid(mu, phi, indexing='ij')
        self.sin_t = np.sqrt(1 - self.mu ** 2)
        self.sin_lat = self.sin_t * np.cos(phi)
        self.cell = (2 / mu_points) * (math.pi / 2 / phi_points) * 4    # solid angle per cell, for all four quarters
        self.ion_mass = np.zeros_like(self.r)        # kg/s per unit lunar radius per steradian, before sunlight
        self.ion_number = np.zeros_like(self.r)      # ions per m^3 per s, before sunlight
        self.molecule_mass = np.zeros_like(self.r)   # kg per m^3 of ballistic molecules
        for s, kg in MOLECULE_KG.items():
            bal, _ = ab.exosphere_density(self.r, shares[s] * self.n_c, self.t_c, self.r_c, mass_kg=kg)
            self.ion_mass += j[s]['ionization'] * kg * bal * self.r ** 2 * MOON_RADIUS ** 3
            self.ion_number += j[s]['ionization'] * bal
            self.molecule_mass += kg * bal
        self.recombined_escape = escaping_share(self.r, self.r_c, RECOMBINATION_EV, ATOM_KG['N2'])
        self.alpha = 2.2e-13 * (300.0 / ELECTRON_K) ** 0.39

    def fates(self, shadow_R, moment, level='central', wind=lr.SOLAR_WIND):
        """Ion production (kg/s) outside the magnetopause, on open field lines, and on closed ones by fate."""
        r_mp = standoff_R(moment, wind) if moment > 0 else 0.0
        opening = pick(MAGNETOSPHERE['open_over_standoff'], level)
        eps = pick(MAGNETOSPHERE['convection_efficiency'], level)
        kappa = pick(MAGNETOSPHERE['drain_factor'], level)
        e_sw = wind['speed_m_s'] * IMF_T
        gm = MOON_GM / MOON_RADIUS
        thermal = 2 * BOLTZMANN * self.t_c / MOLECULE_KG['N2']
        with np.errstate(divide='ignore'):
            boundary = r_mp * (2.0 / (1.0 + self.mu)) ** FLARING
        cos2 = 1 - self.sin_lat ** 2
        keys = ('outside', 'open', 'convected', 'recombined_escaping', 'recombined_staying', 'drained')
        shell = {k: np.zeros_like(self.r) for k in keys}
        exchange = np.zeros_like(self.r)
        for i, ri in enumerate(self.r):
            lit = ri * self.sin_t > shadow_R
            outside = ri > boundary if r_mp > self.r_c else np.ones_like(lit)
            L = ri / cos2
            opened = ~outside & (L > opening * r_mp)
            closed = ~(outside | opened)
            s_d = np.maximum(L * (_field_path(np.sqrt(np.clip(1 - self.r_c / L, 0, None)))
                                  - _field_path(np.abs(self.sin_lat))) * MOON_RADIUS, 1.0)
            nu_d = kappa * math.sqrt(thermal + gm * (1 / self.r_c - 1 / ri)) / s_d
            b_field = (VACUUM_PERMEABILITY * moment / (4 * math.pi * (ri * MOON_RADIUS) ** 3)
                       * np.sqrt(1 + 3 * self.sin_lat ** 2))
            nu_c = eps * e_sw / np.maximum(b_field, 1e-30) / max(r_mp * MOON_RADIUS, 1.0)
            q = self.ion_number[i] * lit
            nu = nu_d + nu_c
            n_e = 2 * q / (nu + np.sqrt(nu ** 2 + 4 * self.alpha * q))
            nu_r = self.alpha * n_e
            total = nu + nu_r
            f_esc = self.recombined_escape[i]
            w = lit * self.cell
            shell['outside'][i] = np.sum(w * outside)
            shell['open'][i] = np.sum(w * opened)
            shell['convected'][i] = np.sum(w * closed * nu_c / total)
            shell['recombined_escaping'][i] = np.sum(w * closed * nu_r * f_esc / total)
            shell['recombined_staying'][i] = np.sum(w * closed * nu_r * (1 - f_esc) / total)
            shell['drained'][i] = np.sum(w * closed * nu_d / total)
            exchange[i] = np.sum(self.cell * outside)
        out = {k: float(np.trapezoid(v * self.ion_mass, self.r)) for k, v in shell.items()}
        flux = wind['density_m3'] * wind['speed_m_s']
        out['charge_exchange_outside'] = float(flux * CHARGE_EXCHANGE_M2 * np.trapezoid(
            exchange * self.molecule_mass * self.r ** 2 * MOON_RADIUS ** 3, self.r))
        out['standoff_R'] = r_mp
        return out


def no_magnetosphere(ions, exchange, level):
    """Loss (kg/s) with the wind reaching the exosphere: escaping ions and the sputtering of those that return."""
    source = ions + exchange * pick(CHARGE_EXCHANGE_SHARE, level) * pick(lr.SW_RANGES['exposure'], level)
    back = pick(lr.SW_RANGES['reimpact'], level)
    sputtered = source * back * pick(lr.SW_RANGES['ion_yield'], level)
    return dict(ions_escaping_kg_s=source * (1 - back), ion_sputtering_kg_s=sputtered,
                plasma_kg_s=source * (1 - back) + sputtered)


def with_magnetosphere(fates, level):
    """Loss (kg/s) with a lunar dipole: ions outside, on open field lines, convected away or recombined to fast atoms."""
    exchange = (fates['charge_exchange_outside'] * pick(CHARGE_EXCHANGE_SHARE, level)
                * pick(lr.SW_RANGES['exposure'], level))
    escaping = fates['outside'] + fates['open'] + fates['convected'] + fates['recombined_escaping']
    return dict(ions_escaping_kg_s=escaping + exchange, charge_exchange_kg_s=exchange,
                ions_kept_kg_s=fates['drained'] + fates['recombined_staying'], plasma_kg_s=escaping + exchange,
                parts_kg_s={k: fates[k] for k in ('outside', 'open', 'convected', 'recombined_escaping',
                                                    'recombined_staying', 'drained')},
                standoff_R=fates['standoff_R'])


def fragments(mixed, separated, level):
    """Fragments escaping (kg/s): diffusive separation and the lower energetic share at the low end."""
    low = sum(ENERGETIC_SHARE[s][0] * separated[s]['fragments'] for s in MOLECULE_KG)
    high = sum(ENERGETIC_SHARE[s][1] * mixed[s]['fragments'] for s in MOLECULE_KG)
    return np.sqrt(low * high) if level == 'central' else (low if level == 'low' else high)


def ions_made(mixed, separated, level):
    low = sum(separated[s]['ions'] for s in MOLECULE_KG)
    high = sum(mixed[s]['ions'] for s in MOLECULE_KG)
    return np.sqrt(low * high) if level == 'central' else (low if level == 'low' else high)


def smallest(radii, values, target):
    return ab.smallest_radius(np.asarray(radii), np.asarray(values), target) if target > 0 else None


def assess(profile, activity, shadow_R, budget=None, september=SEPTEMBER_MOMENT, sweep=True, sweep_shadow_R=None,
           wake_R=None):
    """The sunlit exosphere's losses at one shadow radius, with and without a magnetosphere, and what reduces them;
    with wake_R, also without a magnetosphere with the wind absent within that radius of the Sun-Moon line.

    The moment sweep is taken at sweep_shadow_R (default: shadow_R)."""
    radii = np.array(RADII)
    at = np.r_[shadow_R, radii]
    mixed = sunlit(profile, activity, at)
    separated = sunlit(profile, activity, at, homopause_pa=HOMOPAUSE_PA)
    first = lambda d: {s: {k: float(v[0]) for k, v in d[s].items()} for s in d}
    exchange = {'mixed': charge_exchange(profile), 'separated': charge_exchange(profile, HOMOPAUSE_PA)}
    exo = {'mixed': Exosphere(profile, activity), 'separated': Exosphere(profile, activity, HOMOPAUSE_PA)}
    share = lambda level: 'separated' if level == 'low' else 'mixed'
    result = dict(shadow_radius_R=shadow_R, standoff_september_R=standoff_R(september),
                  sunlit=dict(ions_made_kg_s={lv: float(ions_made(first(mixed), first(separated), lv)) for lv in LEVELS},
                              fragments_escaping_kg_s={lv: float(fragments(first(mixed), first(separated), lv))
                                                       for lv in LEVELS},
                              charge_exchange_upper_bound_kg_s=exchange,
                              exobase_o2_share={'mixed': ab.O2_SHARE,
                                                'separated': exobase_shares(profile, HOMOPAUSE_PA)['O2']}))
    none_, sept = {}, {}
    for lv in LEVELS:
        frag = float(fragments(first(mixed), first(separated), lv))
        ions = float(ions_made(first(mixed), first(separated), lv))
        loss = no_magnetosphere(ions, exchange[share(lv)], lv)
        none_[lv] = dict(loss, fragments_kg_s=frag, total_kg_s=loss['plasma_kg_s'] + frag)
        mag = with_magnetosphere(exo[share(lv)].fates(shadow_R, september, lv), lv)
        sept[lv] = dict(mag, fragments_kg_s=frag, total_kg_s=mag['plasma_kg_s'] + frag)
    result.update(no_magnetosphere=none_, september_magnetosphere=sept)
    if wake_R is not None:
        sheltered = {'mixed': charge_exchange(profile, sheltered_R=wake_R),
                     'separated': charge_exchange(profile, HOMOPAUSE_PA, sheltered_R=wake_R)}
        wake = {}
        for lv in LEVELS:
            frag = float(fragments(first(mixed), first(separated), lv))
            loss = no_magnetosphere(float(ions_made(first(mixed), first(separated), lv)), sheltered[share(lv)], lv)
            wake[lv] = dict(loss, fragments_kg_s=frag, total_kg_s=loss['plasma_kg_s'] + frag)
        result['sunlit']['charge_exchange_in_wake_upper_bound_kg_s'] = sheltered
        result.update(wake_core_R=wake_R, no_magnetosphere_wake=wake)
    # the central loss against the shadow radius, without and with the September moment
    ions_r = np.sqrt(sum(separated[s]['ions'] for s in MOLECULE_KG) * sum(mixed[s]['ions'] for s in MOLECULE_KG))[1:]
    frag_low = sum(ENERGETIC_SHARE[s][0] * separated[s]['fragments'] for s in MOLECULE_KG)[1:]
    frag_high = sum(ENERGETIC_SHARE[s][1] * mixed[s]['fragments'] for s in MOLECULE_KG)[1:]
    frag_r = np.sqrt(frag_low * frag_high)
    none_r = np.array([no_magnetosphere(i, exchange['mixed'], 'central')['plasma_kg_s'] for i in ions_r]) + frag_r
    sept_r = np.array([with_magnetosphere(exo['mixed'].fates(x, september, 'central'), 'central')['plasma_kg_s']
                       for x in radii]) + frag_r
    result['against_radius'] = dict(radius_R=radii.tolist(), no_magnetosphere_kg_s=none_r.tolist(),
                                    september_magnetosphere_kg_s=sept_r.tolist(), fragments_kg_s=frag_r.tolist())
    if budget is not None:
        result['radius_for_a_tenth_of_budget_R'] = dict(
            no_magnetosphere=smallest(radii, none_r, 0.1 * budget),
            september_magnetosphere=smallest(radii, sept_r, 0.1 * budget))
    if sweep:
        # a stand-off inside the exobase leaves the exosphere in the wind: the loss without a magnetosphere
        x = shadow_R if sweep_shadow_R is None else sweep_shadow_R
        frag_x = float(np.interp(x, radii, frag_r))
        none_x = float(np.interp(x, radii, none_r))
        result['moment_sweep'] = dict(
            shadow_radius_R=x, moment_A_m2=list(MOMENTS), standoff_R=[standoff_R(m) for m in MOMENTS],
            loss_kg_s=[none_x if standoff_R(m) <= exo['mixed'].r_c else
                       with_magnetosphere(exo['mixed'].fates(x, m, 'central'), 'central')['plasma_kg_s'] + frag_x
                       for m in MOMENTS])
    return result


def central_losses(profile, activity, shadows_R):
    """The central loss (kg/s) of the sunlit exosphere outside each shadow: without a magnetosphere, with the ring
    fleet's wake sheltering the paths near the Sun-Moon line, and with the September moment."""
    mixed = sunlit(profile, activity, shadows_R)
    separated = sunlit(profile, activity, shadows_R, homopause_pa=HOMOPAUSE_PA)
    exchange = charge_exchange(profile)
    exo = Exosphere(profile, activity)
    out = dict(none=[], september=[], none_wake=[])
    for j, x in enumerate(shadows_R):
        at = lambda d: {s: {k: float(v[j]) for k, v in d[s].items()} for s in d}
        frag = float(fragments(at(mixed), at(separated), 'central'))
        ions = float(ions_made(at(mixed), at(separated), 'central'))
        out['none'].append(no_magnetosphere(ions, exchange, 'central')['plasma_kg_s'] + frag)
        sheltered = charge_exchange(profile, sheltered_R=wake_core_R(x))
        out['none_wake'].append(no_magnetosphere(ions, sheltered, 'central')['plasma_kg_s'] + frag)
        out['september'].append(with_magnetosphere(exo.fates(x, SEPTEMBER_MOMENT, 'central'), 'central')['plasma_kg_s']
                                + frag)
    return out


def allowed_with_exosphere(shield, treatment, activity, shadows=SHADOWS_R, budgets=lr.BUDGETS_KG_S):
    """Largest UV transmission per budget with the sunlit exosphere's central loss added, for each protected radius.

    The loss at each transmission is the loss response's ultraviolet-driven loss (with Earth's tide) at the traced
    state for that radius, plus the exosphere's, and without a magnetosphere also the solar wind's central proton
    sputtering. A protected radius holds only where it covers the thermosphere's heating, light beyond its aperture
    giving at most a tenth of the heat (traced.HEATING_SHARE); beyond that the transmission is not allowed."""
    protons = lr.solar_wind_losses()['central']['proton_sputtering_kg_s']
    out = {key: {} for key in MAGNETOSPHERE_KEYS}
    sweep = {}
    for x in shadows:
        table = []
        for row in lr.euv_sweep(shield, treatment, activity, glow=True, protected_R=x):
            entry = dict(transmission=row['leak_fraction'], uv_kg_s=lr.loss_estimate(row),
                         beyond_share=row.get('beyond_share', math.inf), exobase_radius_R=row.get('exobase_radius_R'),
                         exosphere_kg_s={})
            if row['status'] == 'thermal_column':
                try:
                    q = row['deposited_heat_w_m2']
                    _, profile, _ = solve_column(q, lr.state_config(shield, treatment, activity, x, row['leak_fraction'], q))
                    entry['exosphere_kg_s'] = {k: v[0] for k, v in central_losses(profile, activity, [x]).items()}
                except (ValueError, RuntimeError, FloatingPointError):
                    entry['exosphere_kg_s'] = {}
            table.append(entry)
        sweep[f'{x:g}'] = table
        f = [e['transmission'] for e in table]
        share = [e['beyond_share'] for e in table]
        for key in MAGNETOSPHERE_KEYS:
            loss = [e['uv_kg_s'] + e['exosphere_kg_s'][key] + (protons if key != 'september' else 0.0)
                    if key in e['exosphere_kg_s'] else math.inf for e in table]
            out[key][f'{x:g}'] = {f'{b:g}': _largest_within(f, loss, share, traced.HEATING_SHARE, b) for b in budgets}
    return dict(shield=shield, treatment=treatment, activity=activity, allowed=out, sweep=sweep)


def _largest_within(f, loss, share, limit, budget):
    """Largest transmission whose loss stays within the budget and whose heat from beyond the aperture within its
    share, interpolating the loss log-log and the share against log transmission; None if none."""
    if not (loss[0] <= budget and share[0] <= limit):
        return None
    for i in range(1, len(f)):
        over_loss, over_share = loss[i] > budget, share[i] > limit
        if not (over_loss or over_share):
            continue
        x0, x1 = math.log(f[i - 1]), math.log(f[i])
        ends = []
        if over_loss:
            y0, y1 = math.log(max(loss[i - 1], 1e-30)), math.log(loss[i])
            ends.append(x0 + (math.log(budget) - y0) * (x1 - x0) / (y1 - y0) if math.isfinite(y1) else x0)
        if over_share:
            h0, h1 = share[i - 1], share[i]
            ends.append(x0 + (limit - h0) * (x1 - x0) / (h1 - h0) if math.isfinite(h1) else x0)
        return float(math.exp(min(ends)))
    return float(f[-1])


def summarise(allowed, activities=('quiet', 'solar_maximum')):
    """Allowed transmissions per magnetosphere, shield, protected radius and budget: the range over the upper-air
    treatments and the given solar activities, and how many of those cases cannot meet the budget at any
    transmission."""
    out = {}
    for key in MAGNETOSPHERE_KEYS:
        out[key] = {}
        for shield in lr.SHIELDS:
            rows = [a for a in allowed if a['shield'] == shield and a['activity'] in activities]
            out[key][shield] = {}
            for x in SHADOWS_R:
                out[key][shield][f'{x:g}'] = {}
                for b in lr.BUDGETS_KG_S:
                    values = [a['allowed'][key][f'{x:g}'][f'{b:g}'] for a in rows]
                    held = [v for v in values if v is not None]
                    out[key][shield][f'{x:g}'][f'{b:g}'] = dict(
                        transmission=[min(held), max(held)] if held else None,
                        cases_without=len(values) - len(held), cases=len(values))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    absorption = json.loads((args.out / 'absorption_radius.json').read_text())
    if absorption['schema'] != ab.SCHEMA:
        raise SystemExit(f"unexpected absorption schema {absorption['schema']}")
    fmt = lambda v: f'{v:5.2f}' if v is not None else '  -  '
    rows = []
    shadow = lr.PROTECTED_R
    for c in absorption['cases']:
        _, profile, _ = ab.column(c['shield'], c['treatment'], c['activity'], c['allowed_leak_fraction'])
        row = dict(shield=c['shield'], treatment=c['treatment'], activity=c['activity'], budget_kg_s=c['budget_kg_s'],
                   allowed_leak_fraction=c['allowed_leak_fraction'], exobase_radius_R=c['exobase_radius_R'],
                   exobase_temperature_k=c['exobase_temperature_k'],
                   heating_radius_R=c['heating']['0.1']['radius_R'])
        row.update(assess(profile, c['activity'], shadow, c['budget_kg_s'], sweep=False, wake_R=wake_core_R(shadow)))
        rows.append(row)
        n, m = row['no_magnetosphere'], row['september_magnetosphere']
        tenth = row['radius_for_a_tenth_of_budget_R']
        print(f"{row['shield']:13s} {row['treatment']:11s} {row['activity']:13s} {row['budget_kg_s']:4g} kg/s: "
              f"exobase {row['exobase_radius_R']:.2f} R, shadow {shadow:.2f} R; made {row['sunlit']['ions_made_kg_s']['central']:6.2f}, "
              f"fragments {row['sunlit']['fragments_escaping_kg_s']['central']:6.2f}, exchange "
              f"{row['sunlit']['charge_exchange_upper_bound_kg_s']['mixed']:5.2f}; no field "
              f"{n['low']['total_kg_s']:6.2f}/{n['central']['total_kg_s']:6.2f}/{n['high']['total_kg_s']:6.2f}; "
              f"September {m['low']['total_kg_s']:6.2f}/{m['central']['total_kg_s']:6.2f}/{m['high']['total_kg_s']:6.2f}; "
              f"tenth at {fmt(tenth['no_magnetosphere'])} / {fmt(tenth['september_magnetosphere'])} R")
    allowed = [allowed_with_exosphere(shield, treatment, activity)
               for shield in lr.SHIELDS for treatment in lr.TREATMENTS for activity in lr.ACTIVITY]
    summary = summarise(allowed)
    stress = {name: summarise(allowed, (name,)) for name in lr.STRESS_CASES}
    for key, by_shield in summary.items():
        for shield, by_radius in by_shield.items():
            for x, by_budget in by_radius.items():
                cells = []
                for b, s in by_budget.items():
                    span = ('%.1e-%.1e' % tuple(s['transmission'])) if s['transmission'] else 'none'
                    cells.append(f"{b:>3s} kg/s {span:>15s} ({s['cases_without']} of {s['cases']} without)")
                print(f"{key:9s} {shield:13s} {x:>2s} R: " + '; '.join(cells))
    digest = lambda p: hashlib.sha256((lr.ROOT / p).read_bytes()).hexdigest()[:16]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={p: digest(p) for p in (
            'atmosphere/loss_response/exosphere.py', 'atmosphere/loss_response/absorption.py',
            'atmosphere/loss_response/model.py', 'atmosphere/loss_response/traced.py', 'atmosphere/thermal_column.py',
            'atmosphere/middle_atmosphere/results/limb_heat.json')}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        assumptions=dict(
            rates_s={a: rates(a) for a in lr.ACTIVITY}, heays_2017_quiet_1au_s=HEAYS_1AU_S,
            homopause_pa=HOMOPAUSE_PA, fragment_energy_ev=FRAGMENT_EV, energetic_share=ENERGETIC_SHARE,
            recombination_fragment_ev=RECOMBINATION_EV, charge_exchange_m2=CHARGE_EXCHANGE_M2,
            charge_exchange_share=CHARGE_EXCHANGE_SHARE, imf_t=IMF_T, solar_wind=lr.SOLAR_WIND,
            solar_wind_ranges=lr.SW_RANGES, magnetosphere=MAGNETOSPHERE, flaring=FLARING,
            electron_temperature_k=ELECTRON_K, september_moment_A_m2=SEPTEMBER_MOMENT,
            september_standoff_R=standoff_R(SEPTEMBER_MOMENT), outer_radius_R=OUTER_R,
            protected_radii_R=SHADOWS_R, ring_fleet_distance_km=RING_FLEET_DISTANCE_KM,
            wake_refill_screen_radii=WAKE_REFILL_SCREEN_RADII,
            wake_core_R={f'{x:g}': wake_core_R(x) for x in SHADOWS_R}, heating_share=traced.HEATING_SHARE),
        cases=rows, allowed_summary=summary, stress_summary=stress, allowed=allowed)
    (args.out / 'exosphere_loss.json').write_text(json.dumps(product, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
