"""A bulk mass-flux convection scheme for PlaSim on the Moon, one column at a time in Python, and its offline tests
against the CM1 boxes and the GCM's own columns (the record is massflux.md beside this file).

    climate/gcm/.venv/bin/python -m climate.gcm.massflux        # writes ../results/gcm/massflux_offline.json

The scheme follows Tiedtke (1989) in the form the IFS documentation gives it (Cy36r1 and Cy40r1, Part IV, chapter 6),
with the entrainment and detrainment of Bechtold and others (2008) and the CAPE closure with its boundary-layer term of
Bechtold and others (2014). One bulk updraught rises from a departure level through the column: dry to its
condensation level, then saturated, mixing in the surrounding air and detraining its own, raining part of its
condensate, freezing it between 273 and 250 K, and stopping where its kinetic energy runs out. A test parcel decides
the type: a cloud deeper than 200 hPa convects deep, a thinner one rooted at the ground shallow, and air lifted by the
large-scale ascent from a moist level convects from there (mid-level). A saturated downdraught starts where a mixture
of cloud and surrounding air sinks, with a third of the cloud-base flux, and evaporates the rain to stay saturated. The
cloud-base mass flux comes from the closure: deep convection removes its entraining, density-weighted CAPE (PCAPE),
less the part the boundary layer's own heating supplies (PCAPEbl), over a turnover time; shallow convection carries out
of the sub-cloud layer the moist static energy its other processes bring in; mid-level convection carries the
large-scale mass flux at its base. The tendencies are the convergence of the updraught's and downdraught's fluxes
of liquid-ice static energy and total water, against the compensating motion of the surrounding air taken upwind,
with the rain's formation, evaporation and melting as the only sources, so the column keeps its energy and water
exactly: the heat it gains is the latent heat of the rain and snow at the ground, and the water it loses is that
precipitation. Detrained condensate evaporates where it leaves the cloud, since PlaSim carries no cloud water. With
melt 'plasim' the falling water takes the phase PlaSim's mkrain gives it, and the heating goes over PlaSim's moist heat
capacity; massflux_fortran.py holds that form as Fortran for PlaSim, tested against this copy.

Lunar choices (each recorded in PARAMETERS with its source):
- Entrainment per metre follows the bulk rates measured on the CM1 sea box's clouds (box_entrainment.py) and the
  zero-buoyancy plume that reproduces its upper air (box_upper_air.py); per unit pressure it is eps times the local
  scale height R Tv / g, the per-e-folding count PlaSim's simplified Betts-Miller scheme uses (eps_bm).
- Depths, departure ranges and mixed layers are in pressure, which a column stretched by the ratio of the gravities
  keeps; lengths and Earth's constant rates per metre are divided by that ratio where no CM1 measurement stands for
  them, and each such value is marked as stretched.
- The surface test parcel is the IFS's in such a stretched column: the air of the IFS's lowest model level, 10 m up on
  Earth and 60 m on the Moon, carried down from PlaSim's lowest level by surface-layer similarity, with the IFS's
  excess from the surface fluxes at that height (the second round, massflux.md). It starts at the IFS's 1 m/s, or
  with start 'boundary_layer' at c_start times the boundary layer's convective velocity w*, the speed of a grid cell's
  strongest thermals (the third round).
- Times follow the column: the adjustment time is the cloud's own turnover time, H over its mean updraught speed, with
  no Earth bound on it; the rain forms at Earth's rates per second, as the gravity pair's microphysics does.

Arrays run top first, as PlaSim's levels do: full level k lies between half levels k (above) and k + 1 (below), half
level 0 is the model top and half level n the ground.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from climate.gcm import sbm

HERE = Path(__file__).resolve().parent

PLASIM = sbm.PLASIM
G_EARTH = 9.80665                  # standard gravity (shared/constants.json physics.standard_gravity_m_s2)

# The scheme's parameters, each with its source. Lengths per metre are lunar values; 'stretched' marks an Earth value
# divided by g_earth/g_moon = 6.04 because no CM1 measurement stands for it.
PARAMETERS = dict(
    # updraught entrainment (per metre): eps = eps_bg + eps_nb exp(-(z - z_base) / l_nb) above the cloud base, from
    # the CM1 boxes' bulk rates (massflux.md): 0.06-0.1e-3 /m through the middle and upper troposphere of both boxes'
    # deep clouds and 0.5-1e-4 /m in the zero-buoyancy plume that reproduces the sea box's upper air, 0.2-0.6e-3 /m
    # within about 2 km of the base wherever the base lies; shallow clouds entrain shallow_factor times as fast (the
    # IFS's f_eps = 2; the sea box's night cores 0.18-0.57e-3 /m)
    eps_bg=0.8e-4,
    eps_nb=2.0e-4,
    l_nb=2000.0,
    shallow_factor=2.0,
    rh_dependence=False,           # Bechtold and others (2008, 2014): eps scaled by (1.3 - RH) / (1.3 - rh_ref)
    rh_ref=0.85,
    # turbulent detrainment: delta = eps (1.6 - RH), the IFS's form for shallow clouds (Cy40r1 (6.9); Bechtold and
    # others 2014, (A.3)), taken for every type because CM1's deep clouds detrain nearly as fast as they entrain above
    # 12 km (massflux.md); delta_form 'ifs' gives deep clouds the IFS's delta_deep (1.6 - RH) instead (Cy40r1 (6.8),
    # 0.75e-4 /m stretched)
    delta_form='ratio',
    delta_deep=0.75e-4 * 1.6242 / G_EARTH,
    # the updraught's kinetic energy (IFS (6.10)): virtual mass, turbulence factor, drag
    gamma_vm=0.5, f_turb=2.0, beta_cd=1.875 * 0.506,
    w_base_deep=1.0,               # m/s at the departure level, deep and mid-level tests (IFS)
    dt_deep=0.2, dq_deep=1.0e-4,   # the deep test parcel's excess (IFS (6.21))
    deep_dp=20000.0,               # Pa: a cloud thicker than this convects deep (IFS)
    departure_dp=35000.0,          # Pa: deep test parcels depart up to this far above the ground (IFS)
    source_dp=0.0,                 # Pa: the lowest layer averaged for a surface-rooted deep parcel (the IFS's 50 hPa;
                                   # 0 takes the lowest level's own air, as CM1's updraughts carry it, massflux.md)
    sub_factor=0.4,                # test parcels mix in sub-cloud air at this share of the deep rate at the base
                                   # (the IFS Cy40r1 deep test's 0.4, (6.21))
    # the surface test parcel (massflux.md, the second round): the IFS's, in a column stretched by the ratio of the
    # gravities. Its air is that of the IFS's lowest model level, 10 m above Earth's ground and 60 m above the Moon's,
    # carried down from PlaSim's lowest level by surface-layer similarity (surface_layer_air), and its excess is the
    # IFS's from the surface fluxes (Cy40r1 (6.19)) with sigma_w at that height ((6.20))
    source_height_earth=10.0,      # m at Earth's gravity, stretched by g_earth / g (0: the lowest level's own air)
    ustar=0.1,                     # m/s, the IFS's friction velocity in sigma_w (6.20), also taken for the similarity
    excess_scale='source',         # sigma_w's height: 'source', the parcel's source height, or 'mixed_layer', the
                                   # first round's convective velocity w* of the sub-cloud layer
    stability='free_convection',   # the similarity's stability function in unstable air: 'free_convection',
                                   # phi_h = (1 - 16 zeta)^(-1/3); 'dyer_hicks', the IFS's (1 - 16 zeta)^(-1/2); 'neutral'
    # the surface deep test parcel's starting speed (massflux.md, the third round): 'ifs', the larger of w_base_deep
    # and the excess's sigma_w; 'boundary_layer', c_start times the boundary layer's convective velocity w* from the
    # surface buoyancy flux and the depth of its dry convective layer (boundary_layer_velocity), nil without a rising
    # surface buoyancy flux. c_start stands for the strongest thermals of a grid cell, the 99.9th percentile of the
    # mixed layer's vertical wind: 2.0-2.3 w* between 0.3 and 0.9 z_i in CM1's 2 km box of the sea's night, whose
    # spread peaks at Lenschow and others' (1980) 0.68 w*
    start='ifs',
    c_start=2.2,
    # microphysics (IFS (6.35)-(6.37)), at Earth's rates per second
    c0=1.4e-3, l_crit=5.0e-4, l_min=3.0e-4, w_max_precip=10.0, t_bf=268.16, t_ice=250.16,
    # downdraughts (IFS (6.13), (6.15)): eta of the cloud-base flux at the level of free sinking; their mixing
    # stretched from 2e-4 /m; they detrain over the lowest 60 hPa
    eta_down=0.3,
    eps_down=2.0e-4 * 1.6242 / G_EARTH,
    down_bottom_dp=6000.0,
    # closure
    tau_factor=5.0,                # tau = tau_factor H / w-bar (the IFS's alpha_nT): CM1's convection consumes the
                                   # scheme's PCAPE over 1.6-6.9 turnover times in its four quarters of deep convection,
                                   # median 5 (massflux.md); the IFS's 1 + 264/n at the lunar grid's Earth-equivalent
                                   # truncation, 77 for T21, gives 4.4
    tau_min=720.0,                 # s, the IFS's lower bound (a numerical floor, not an Earth time scale)
    t_star='cloud',                # T* of PCAPEbl (Bechtold and others 2014, (13)): g H / cp, H the cloud's depth
                                   # ('cloud'), or a number in K (the memo's 1 K casts the scale into alpha)
    elevated_ascent=True,          # convection rooted above the ground departs only under large-scale ascent
    u_bl_min=2.0,                  # m/s, the sub-cloud wind's floor in tau_bl over water (IFS)
    midlevel_rh=0.8,               # mid-level convection departs from levels at least this moist (IFS)
    midlevel_top=0.25,             # sigma: mid-level departure levels lie below this (Earth's 10 km, as pressure)
    cfl=1.0,                       # the mass flux's limit, in units of dp / (g dt) of the thinnest layer it crosses
    substeps=2,                    # steps a half layer is divided into for the plume
)


# ----------------------------------------------------------------------------------------------------- thermodynamics

def qsat(t: float, p: float, c=PLASIM) -> float:
    """PlaSim's saturation specific humidity (kg/kg) at t (K) and p (Pa), over water, as sbm.saturation."""
    q = c['rdbrv'] * c['ra1'] * math.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
    q = min(q, c['rdbrv'])
    return q / (1.0 - (1.0 / c['rdbrv'] - 1.0) * q)


def dqsat_dt(t: float, p: float, c=PLASIM) -> float:
    q0 = c['rdbrv'] * c['ra1'] * math.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
    q0 = min(q0, c['rdbrv'])
    cor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * q0)
    return q0 * cor * cor * c['ra2'] * (c['TMELT'] - c['ra4']) / (t - c['ra4']) ** 2


def ice_share(t: float, prm=PARAMETERS, c=PLASIM) -> float:
    """The frozen share of condensate at t (K): nil above freezing, all below t_ice, and the IFS's quadratic mixed
    phase between (IFS Cy36r1, (7.6))."""
    if t >= c['TMELT']:
        return 0.0
    if t <= prm['t_ice']:
        return 1.0
    a = (t - prm['t_ice']) / (c['TMELT'] - prm['t_ice'])
    return 1.0 - a * a


def latent(t: float, prm=PARAMETERS, c=PLASIM) -> float:
    """The latent heat (J/kg) of condensate at t: vaporization for its liquid share, sublimation for its ice."""
    return c['ALV'] + ice_share(t, prm, c) * (c['ALS'] - c['ALV'])


def adjust(sl: float, qt: float, z: float, p: float, g: float, prm=PARAMETERS, c=PLASIM):
    """Temperature, vapour and condensate (K, kg/kg, kg/kg) of air with liquid-ice static energy
    sl = cpd T + g z - L(T) l and total water qt at height z and pressure p: unsaturated air keeps its water as vapour,
    saturated air holds saturation and the rest as condensate, its ice share following its temperature."""
    cp = c['acpd']
    td = (sl - g * z) / cp
    if qt <= qsat(td, p, c):
        return td, qt, 0.0
    t = td
    for _ in range(12):                     # Newton on cp T + g z - L(T) (qt - qs(T)) - sl = 0
        qs = qsat(t, p, c)
        lv = latent(t, prm, c)
        dl = 0.0
        if prm['t_ice'] < t < c['TMELT']:
            a = (t - prm['t_ice']) / (c['TMELT'] - prm['t_ice'])
            dl = -2.0 * a / (c['TMELT'] - prm['t_ice']) * (c['ALS'] - c['ALV'])
        f = cp * t + g * z - lv * (qt - qs) - sl
        df = cp - dl * (qt - qs) + lv * dqsat_dt(t, p, c)
        step = f / df
        t -= step
        if abs(step) < 1.0e-7:
            break
    qs = qsat(t, p, c)
    return t, qs, max(qt - qs, 0.0)


def virtual(t: float, qv: float, c=PLASIM) -> float:
    return t * (1.0 + (1.0 / c['rdbrv'] - 1.0) * qv)


# ------------------------------------------------------------------------------------------------------------ geometry

def heights(t, q, p_full, p_half, g: float, c=PLASIM):
    """Heights (m above the ground) of the full and half levels, top first, from the virtual temperatures:
    each layer's thickness R Tv ln(p_below / p_above) / g, the full level placed at its pressure within the layer."""
    t, q, p_full, p_half = (np.asarray(a, float) for a in (t, q, p_full, p_half))
    n = t.size
    tv = t * (1.0 + (1.0 / c['rdbrv'] - 1.0) * q)
    zh = np.zeros(n + 1)
    zf = np.zeros(n)
    for k in range(n - 1, -1, -1):
        r = c['gascon'] * tv[k] / g
        zf[k] = zh[k + 1] + r * math.log(p_half[k + 1] / p_full[k])
        top = p_half[k] if p_half[k] > 0.0 else p_full[k] * 0.5
        zh[k] = zh[k + 1] + r * math.log(p_half[k + 1] / top)
    return zf, zh


class Column:
    """A column as the scheme reads it: full-level temperature, vapour, pressure and heights, half-level pressure and
    heights, the environment along the plume's path (full and half levels from the ground up) and each layer's mass."""

    def __init__(self, t, q, p_full, p_half, g: float, c=PLASIM, z_full=None, z_half=None):
        self.t, self.q = np.asarray(t, float), np.asarray(q, float)
        self.pf, self.ph = np.asarray(p_full, float), np.asarray(p_half, float)
        self.n = self.t.size
        self.g, self.c = g, c
        if z_full is None:
            self.zf, self.zh = heights(self.t, self.q, self.pf, self.ph, g, c)
        else:
            self.zf, self.zh = np.asarray(z_full, float), np.asarray(z_half, float)
        self.dp = np.diff(self.ph)                         # each layer's thickness (Pa), top first
        self.mass = self.dp / g                            # kg/m2
        self.ps = float(self.ph[-1])
        cp = c['acpd']
        self.s = cp * self.t + g * self.zf                 # dry static energy, J/kg
        self.tv = np.array([virtual(a, b, c) for a, b in zip(self.t, self.q)])
        self.qs = np.array([qsat(a, b, c) for a, b in zip(self.t, self.pf)])
        self.rh = self.q / self.qs
        # the environment at half levels 1..n-1, linear in log pressure between the full levels on each side
        self.th, self.qh = np.zeros(self.n + 1), np.zeros(self.n + 1)
        for i in range(1, self.n):
            w = math.log(self.ph[i] / self.pf[i - 1]) / math.log(self.pf[i] / self.pf[i - 1])
            self.th[i] = (1.0 - w) * self.t[i - 1] + w * self.t[i]
            self.qh[i] = (1.0 - w) * self.q[i - 1] + w * self.q[i]
        self.th[self.n], self.qh[self.n] = self.t[-1], self.q[-1]
        self.th[0], self.qh[0] = self.t[0], self.q[0]

    def path(self, k_start: int):
        """The plume's points from full level k_start up: (kind, index, z, p, environment T, q), kind 'f' for a full
        level and 'h' for a half level, ending at the highest full level."""
        pts = []
        for k in range(k_start, -1, -1):
            pts.append(('f', k, self.zf[k], self.pf[k], self.t[k], self.q[k]))
            if k > 0:
                pts.append(('h', k, self.zh[k], self.ph[k], self.th[k], self.qh[k]))
        return pts


# --------------------------------------------------------------------------------------------------------------- plume

def mix(phi: float, env_a: float, env_b: float, eps: float, dz: float) -> float:
    """A plume's value after rising dz while mixing in surrounding air at eps per metre, the surrounding value linear
    from env_a to env_b: the exact solution of d phi/dz = eps (phi_env - phi)."""
    if eps <= 0.0 or dz <= 0.0:
        return phi
    grad = (env_b - env_a) / dz
    decay = math.exp(-eps * dz)
    return env_b - grad / eps + (phi - env_a + grad / eps) * decay


def entrainment(kind: str, z: float, rh: float, z_base: float, prm=PARAMETERS) -> float:
    """The updraught's fractional entrainment (per metre) at height z, z_base being its cloud base and rh the
    surrounding relative humidity: eps_bg + eps_nb exp(-(z - z_base) / l_nb), scaled by (1.3 - rh) / (1.3 - rh_ref)
    when rh_dependence is set, and shallow_factor times as fast for shallow clouds."""
    eps = prm['eps_bg'] + prm['eps_nb'] * math.exp(-max(z - z_base, 0.0) / prm['l_nb'])
    if prm['rh_dependence']:
        eps *= max(1.3 - min(rh, 1.0), 0.3) / (1.3 - prm['rh_ref'])
    if kind == 'shallow':
        eps *= prm['shallow_factor']
    return eps


def detrainment(kind: str, eps: float, rh: float, prm=PARAMETERS) -> float:
    """The updraught's turbulent fractional detrainment (per metre)."""
    if kind == 'shallow' or prm['delta_form'] == 'ratio':
        return eps * (1.6 - min(rh, 1.0))
    return prm['delta_deep'] * (1.6 - min(rh, 1.0))


def conversion(l: float, t: float, w: float, prm=PARAMETERS, c=PLASIM) -> float:
    """The rate (per metre of ascent) at which condensate l (kg/kg) turns to rain or snow at temperature t in an
    updraught rising at w: c0 / (0.75 w) (1 - exp(-(l / l_crit)^2)), with the Bergeron-Findeisen enhancement below
    t_bf, nil below l_min (IFS (6.35)-(6.37))."""
    if l <= prm['l_min']:
        return 0.0
    cbf = 1.0
    if t < prm['t_bf']:
        cbf = 1.0 + 0.5 * math.sqrt(min(prm['t_bf'] - t, prm['t_bf'] - prm['t_ice']))
    lcrit = prm['l_crit'] / cbf
    w = min(max(w, 0.05), prm['w_max_precip'])
    return prm['c0'] * cbf / (0.75 * w) * (1.0 - math.exp(-(l / lcrit) ** 2))


def ascend(col: Column, pts: list, state: dict, kind: str, rain=True, entrain_below_base: bool = False,
           prm=PARAMETERS):
    """One plume along the path points pts (Column.path, the ground up) from the state at pts[0] (its liquid-ice
    static energy sl, total water qt, kinetic energy kin, and mass flux m, normally 1): below its condensation level
    it rises dry, undiluted unless entrain_below_base; above it, it mixes in the surrounding air, detrains its own and
    rains (rain True: by conversion; 'half': half the condensate it forms, the IFS's deep test; False: none), and it
    stops where its kinetic energy runs out. Returns its points (dicts of the plume's state) from the start to the
    top, the index of its first saturated point (None if it never saturates) and the top height (m)."""
    c, g = col.c, col.g
    cp = c['acpd']
    z0, p0, te0, qe0 = pts[0][2], pts[0][3], pts[0][4], pts[0][5]
    s_l, q_t, m, kin = state['sl'], state['qt'], state.get('m', 1.0), state['kin']
    t_a, qv_a, l_a = adjust(s_l, q_t, z0, p0, g, prm, c)
    tv0 = virtual(te0, qe0, c)
    out = [dict(kind=pts[0][0], index=pts[0][1], z=z0, p=p0, t=t_a, qv=qv_a, l=l_a, sl=s_l, qt=q_t, m=m, kin=kin,
                b=g * (virtual(t_a, qv_a, c) - l_a * t_a - tv0) / tv0, rain=0.0, snow=0.0, eps=0.0, delta=0.0)]
    base = 0 if l_a > 0.0 else None
    z_base = state.get('z_base', z0)
    top_z = None
    nsub = max(1, int(prm['substeps']))
    fac = 1.0 / (prm['f_turb'] * (1.0 + prm['gamma_vm']))
    for j in range(1, len(pts)):
        kind_b, idx_b, z_b, p_b, te_b, qe_b = pts[j]
        _, _, z_a, p_a, te_a, qe_a = pts[j - 1]
        stop = False
        for s in range(nsub):
            f0, f1 = s / nsub, (s + 1) / nsub
            za, zb = z_a + f0 * (z_b - z_a), z_a + f1 * (z_b - z_a)
            pa = math.exp(math.log(p_a) + f0 * math.log(p_b / p_a))
            pb = math.exp(math.log(p_a) + f1 * math.log(p_b / p_a))
            tea, teb = te_a + f0 * (te_b - te_a), te_a + f1 * (te_b - te_a)
            qea, qeb = qe_a + f0 * (qe_b - qe_a), qe_a + f1 * (qe_b - qe_a)
            dz = zb - za
            rhm = 0.5 * (qea + qeb) / qsat(0.5 * (tea + teb), math.sqrt(pa * pb), c)
            cloudy = base is not None
            if cloudy:
                eps = entrainment(kind, 0.5 * (za + zb), rhm, z_base, prm)
                delta = detrainment(kind, eps, rhm, prm)
            elif entrain_below_base:                  # a test parcel's thermal mixing in the sub-cloud layer
                eps = prm['sub_factor'] * (prm['eps_bg'] + prm['eps_nb'])
                delta = eps
            else:
                eps, delta = 0.0, 0.0
            l_prev = out[-1]['l']
            s_l = mix(s_l, cp * tea + g * za, cp * teb + g * zb, eps, dz)
            q_t = mix(q_t, qea, qeb, eps, dz)
            t_b, qv_b, l_b = adjust(s_l, q_t, zb, pb, g, prm, c)
            formed = 0.0
            if l_b > 0.0 and rain:
                if rain == 'half':
                    formed = 0.5 * max(l_b - l_prev * math.exp(-eps * dz), 0.0)
                else:
                    a = conversion(l_b, t_b, math.sqrt(2.0 * max(kin, 0.0)), prm, c)
                    if a > 0.0:
                        l0 = l_prev * math.exp(-eps * dz)
                        decay = math.exp(-a * dz)
                        kept = l0 * decay + (l_b - l0) / (a * dz) * (1.0 - decay)
                        formed = max(l_b - max(kept, 0.0), 0.0)
                if formed > 0.0:
                    q_t -= formed
                    s_l += latent(t_b, prm, c) * formed
                    t_b, qv_b, l_b = adjust(s_l, q_t, zb, pb, g, prm, c)
            tvb = virtual(teb, qeb, c)
            b_b = g * (virtual(t_b, qv_b, c) - l_b * t_b - tvb) / tvb
            # kinetic energy, with drag from mixing at the larger of the rates (IFS (6.10), (6.54))
            mu = max(eps, delta)
            bm = 0.5 * (out[-1]['b'] + b_b)
            if mu > 0.0:
                k2 = 2.0 * (1.0 + prm['beta_cd']) * mu
                dec = math.exp(-k2 * dz)
                kin_b = kin * dec + fac * bm * (1.0 - dec) / k2
            else:
                kin_b = kin + fac * bm * dz
            m_b = m * math.exp((eps - delta) * dz)
            if cloudy and b_b < 0.0 and kin_b > 0.0 and kin > 0.0:
                m_b = min(m_b, m * math.sqrt(min(kin_b / kin, 1.0)))     # organized detrainment (IFS (6.12))
            if kin_b <= 0.0:
                top_z = za + (kin / (kin - kin_b) if kin > kin_b else 0.0) * dz
                stop, m_b = True, 0.0
            share = ice_share(t_b, prm, c)
            mean_m = 0.5 * (m + m_b)
            out.append(dict(kind=kind_b if s == nsub - 1 else 's', index=idx_b, z=zb, p=pb, t=t_b, qv=qv_b, l=l_b,
                            sl=s_l, qt=q_t, m=m_b, kin=max(kin_b, 0.0), b=b_b, rain=mean_m * formed * (1.0 - share),
                            snow=mean_m * formed * share, eps=eps, delta=delta))
            m, kin = m_b, kin_b
            if base is None and l_b > 0.0:
                base = len(out) - 1
                z_base = zb
            if stop:
                break
        if stop:
            break
    if top_z is None:
        top_z = out[-1]['z']
    return out, base, top_z


# ------------------------------------------------------------------------------------------- test parcels and the type

def mixed_layer_depth(col: Column) -> float:
    """The height (m) of the lowest level's air's condensation level, lifted dry with its vapour kept: the depth of the
    sub-cloud layer its thermals cross."""
    c, g = col.c, col.g
    t0, q0 = col.t[-1], col.q[-1]
    s0 = c['acpd'] * t0 + g * col.zf[-1]
    for k in range(col.n - 1, -1, -1):
        for z, p in ((col.zf[k], col.pf[k]), (col.zh[k], col.ph[k])):
            if q0 >= qsat((s0 - g * z) / c['acpd'], p, c):
                return float(z)
    return float(col.zh[0])


def source_height(col: Column, prm=PARAMETERS) -> float:
    """The height (m) whose air the surface test parcel takes: the IFS's lowest model level, source_height_earth above
    Earth's ground, stretched by the ratio of Earth's gravity to the column's (60 m on the Moon); nil when it is 0."""
    return prm['source_height_earth'] * G_EARTH / col.g


def surface_excess(col: Column, shf: float, lhf: float, prm=PARAMETERS):
    """The surface test parcels' excess, temperature (K) and vapour (kg/kg), and the velocity scale (m/s) that starts
    the shallow test parcel, from the surface's sensible and latent heat fluxes shf and lhf (W/m2, upward): 1.5 times
    each flux over rho sigma_w (IFS Cy40r1 (6.19)). With excess_scale 'source', sigma_w is the IFS's (6.20), Holtslag
    and Moeng's (1991) profile of the vertical velocity's spread in a convective boundary layer at the height z of
    the parcel's source air (the source height, or the lowest level's if that is lower): (1.2 u*^3 + 1.5 kappa g z /
    T (w'theta_v')_s)^(1/3), with z / z_i cancelling against w*^3's z_i. With 'mixed_layer' it is the first round's
    reading, the sub-cloud layer's convective velocity w* = (g / T z_i (w'theta_v')_s)^(1/3), at least 0.1 m/s."""
    c = col.c
    t, q = col.t[-1], col.q[-1]
    rho = col.ph[-1] / (c['gascon'] * virtual(t, q, c))
    lv = c['ALV']
    flux = (shf / c['acpd'] + (1.0 / c['rdbrv'] - 1.0) * t * lhf / lv) / rho      # kinematic buoyancy flux, K m/s
    if prm['excess_scale'] == 'source':
        z = col.zf[-1]
        if prm['source_height_earth'] > 0.0:
            z = min(source_height(col, prm), z)
        ws = (1.2 * prm['ustar'] ** 3 + 1.5 * 0.4 * col.g / t * z * max(flux, 0.0)) ** (1.0 / 3.0)
    else:
        zi = mixed_layer_depth(col)
        ws = max(col.g / t * zi * max(flux, 0.0), 1.0e-6) ** (1.0 / 3.0)
        ws = max(ws, 0.1)
    return 1.5 * max(shf, 0.0) / (rho * c['acpd'] * ws), 1.5 * max(lhf, 0.0) / (rho * lv * ws), ws


def psi_heat(zeta: float, form: str = 'free_convection') -> float:
    """The integrated stability function for heat and vapour in unstable air (zeta < 0), psi(zeta) = the integral from
    0 to zeta of (1 - phi_h(x)) / x: for the free-convection form phi_h = (1 - 16 zeta)^(-1/3) (Carl and others 1973;
    Kader and Yaglom 1990), whose exponent free convection sets, as the lunar surface layer is (zeta near -10 at 0.3 km
    over the sea); for 'dyer_hicks' the IFS's (1 - 16 zeta)^(-1/2) (Paulson 1970's integral); nil for 'neutral'."""
    if form == 'neutral':
        return 0.0
    if form == 'dyer_hicks':
        return 2.0 * math.log((1.0 + math.sqrt(1.0 - 16.0 * zeta)) / 2.0)
    y = (1.0 - 16.0 * zeta) ** (1.0 / 3.0)
    return (1.5 * math.log((y * y + y + 1.0) / 3.0) - math.sqrt(3.0) * math.atan((2.0 * y + 1.0) / math.sqrt(3.0))
            + math.pi / math.sqrt(3.0))


def surface_layer_air(col: Column, shf: float, lhf: float, prm=PARAMETERS):
    """The temperature and vapour (K, kg/kg) by which the air at the source height exceeds the lowest level's in a
    convective surface layer, by Monin-Obukhov similarity from the surface fluxes (W/m2, upward) and the IFS's friction
    velocity u*: potential temperature and vapour differ between the lowest level's height z_n and the source height
    z_s by (w'theta')_s / u* (or (w'q')_s / u*) times (ln(z_n / z_s) - psi(z_n / L) + psi(z_s / L)) / kappa, L the
    Obukhov length; the potential temperature's difference is carried to the lowest level's pressure. Nil when the
    surface's virtual heat flux does not rise (no surface-rooted convection then; the excess is nil too), when the
    source height is 0, or when it is not below the lowest level."""
    z_s, z_n = source_height(col, prm), col.zf[-1]
    if z_s <= 0.0 or z_s >= z_n:
        return 0.0, 0.0
    c = col.c
    t, q = col.t[-1], col.q[-1]
    rho = col.ph[-1] / (c['gascon'] * virtual(t, q, c))
    ust = prm['ustar']
    wt = shf / (rho * c['acpd'])
    wq = lhf / (rho * c['ALV'])
    exner = (col.pf[-1] / col.ps) ** (c['gascon'] / c['acpd'])
    eps1 = 1.0 / c['rdbrv'] - 1.0
    theta = t / exner
    wtv = wt * (1.0 + eps1 * q) + eps1 * theta * wq                   # the surface's virtual heat flux, K m/s
    if wtv <= 0.0:
        return 0.0, 0.0
    obukhov = -ust ** 3 * theta * (1.0 + eps1 * q) / (0.4 * col.g * wtv)
    form = prm.get('stability', 'free_convection')
    f = (math.log(z_n / z_s) - psi_heat(z_n / obukhov, form) + psi_heat(z_s / obukhov, form)) / 0.4
    return wt / ust * f * exner, wq / ust * f


def boundary_layer_velocity(col: Column, shf: float, lhf: float, dt_l: float, dq_l: float):
    """The boundary layer's convective velocity w* (m/s; Deardorff's, as Lenschow and others 1980 define it, their
    eq. 12) and the depth z_i (m) it takes, from the surface fluxes (W/m2, upward) and the column:
    w* = (g / theta_v (w'theta_v')_s z_i)^(1/3), with the surface's
    kinematic virtual heat flux as surface_excess takes it and theta_v the surface parcel's source air (the lowest
    level's with the similarity's difference dt_l, dq_l). z_i is the depth of the dry convective layer that air rises
    through: the height where the surroundings' virtual potential temperature first exceeds its own (the parcel
    method), or where it saturates lifted dry with its vapour kept, whichever is lower, each linear in height between
    full levels. Nil, with no depth, without a rising surface buoyancy flux."""
    c = col.c
    t, q = col.t[-1], col.q[-1]
    rho = col.ph[-1] / (c['gascon'] * virtual(t, q, c))
    eps1 = 1.0 / c['rdbrv'] - 1.0
    flux = (shf / c['acpd'] + eps1 * t * lhf / c['ALV']) / rho
    if flux <= 0.0:
        return 0.0, None
    kap = c['gascon'] / c['acpd']
    n = col.n
    t0, q0 = t + dt_l, q + dq_l
    thv0 = virtual(t0, q0, c) * (col.ps / col.pf[-1]) ** kap
    thvb = virtual(t, q, c) * (col.ps / col.pf[-1]) ** kap
    if thvb >= thv0:
        z_t = col.zf[-1]
    else:
        z_t = col.zf[0]
        for k in range(n - 2, -1, -1):
            thva = virtual(col.t[k], col.q[k], c) * (col.ps / col.pf[k]) ** kap
            if thva > thv0:
                z_t = col.zf[k + 1] + (thv0 - thvb) / (thva - thvb) * (col.zf[k] - col.zf[k + 1])
                break
            thvb = thva
    db = qsat(t0, col.pf[-1], c) - q0
    if db <= 0.0:
        z_l = col.zf[-1]
    else:
        z_l = col.zf[0]
        for k in range(n - 2, -1, -1):
            da = qsat(t0 * (col.pf[k] / col.pf[-1]) ** kap, col.pf[k], c) - q0
            if da <= 0.0:
                z_l = col.zf[k + 1] + db / (db - da) * (col.zf[k] - col.zf[k + 1])
                break
            db = da
    z_i = min(z_t, z_l)
    return (col.g / thv0 * flux * z_i) ** (1.0 / 3.0), float(z_i)


def surface_start(col: Column, shf: float, lhf: float, ws: float, dt_l: float, dq_l: float, prm=PARAMETERS) -> float:
    """The surface deep test parcel's starting speed (m/s): with start 'ifs' the larger of w_base_deep and the
    excess's velocity scale ws (IFS (6.20)-(6.21)); with 'boundary_layer' c_start times the boundary layer's w*
    (boundary_layer_velocity), with no floor."""
    if prm.get('start', 'ifs') == 'boundary_layer':
        return prm['c_start'] * boundary_layer_velocity(col, shf, lhf, dt_l, dq_l)[0]
    return max(prm['w_base_deep'], ws)


def source_layer(col: Column, k: int, dp_mix: float):
    """Temperature and vapour of a departure level's air: the level's own, or for the lowest level the mean
    potential temperature and vapour by mass of the lowest dp_mix of the column carried to the lowest level's
    pressure (the IFS's 50 hPa source layer)."""
    c = col.c
    if k != col.n - 1 or dp_mix <= 0.0:
        return float(col.t[k]), float(col.q[k])
    kap = c['gascon'] / c['acpd']
    top = col.ps - dp_mix
    w_sum, th_sum, q_sum = 0.0, 0.0, 0.0
    for j in range(col.n - 1, -1, -1):
        lo, hi = col.ph[j + 1], max(col.ph[j], top)
        if hi >= lo:
            break
        w = lo - hi
        w_sum += w
        th_sum += w * col.t[j] * (col.ps / col.pf[j]) ** kap
        q_sum += w * col.q[j]
        if col.ph[j] <= top:
            break
    return th_sum / w_sum * (col.pf[-1] / col.ps) ** kap, q_sum / w_sum


def first_half_level(points: list, base: int):
    """The index of the first half level at or above point base, where the cloud's fluxes start."""
    for j in range(base, len(points)):
        if points[j]['kind'] == 'h':
            return j
    return None


def top_pressure(points: list, top_z: float) -> float:
    """The pressure at the plume's top height top_z, log pressure linear in height between the plume's points (which
    rise along the list), the last point's beyond them."""
    z = [x['z'] for x in points]
    return float(np.exp(np.interp(top_z, z, [math.log(x['p']) for x in points])))


def trial(col: Column, k: int, t0: float, q0: float, w0: float, kind: str, rain, prm=PARAMETERS) -> dict:
    """One test parcel from full level k: its points, condensation (cloud base) and top, the cloud's depth in
    pressure from its base to its top, and whether it saturates at a half level with kinetic energy to spare."""
    pts = col.path(k)
    state = dict(sl=col.c['acpd'] * t0 + col.g * col.zf[k], qt=q0, kin=0.5 * w0 * w0, m=1.0)
    points, base, top_z = ascend(col, pts, state, kind, rain=rain, entrain_below_base=True, prm=prm)
    out = dict(k=k, points=points, base=base, top_z=top_z, ok=False)
    if base is None:
        return out
    jb = first_half_level(points, base)
    if jb is None or points[jb]['kin'] <= 0.0:
        return out
    p_top = top_pressure(points, top_z)
    out.update(ok=True, base_point=jb, base_half=points[jb]['index'], p_base=points[jb]['p'], p_top=p_top,
               depth=points[jb]['p'] - p_top)
    return out


def classify(col: Column, shf: float, lhf: float, omega=None, prm=PARAMETERS) -> dict:
    """The type of convection the column takes, by the IFS's sequence: a surface parcel with the surface fluxes'
    excess gives a shallow candidate; deep test parcels depart from the lowest level (its lowest 50 hPa mixed) and each
    level above it up to departure_dp above the ground, and the first whose cloud is thicker than deep_dp makes the
    column deep; failing both, a moist level under large-scale ascent (omega < 0, Pa/s) convects from there."""
    c = col.c
    dt_s, dq_s, ws = surface_excess(col, shf, lhf, prm)
    dt_l, dq_l = surface_layer_air(col, shf, lhf, prm)
    t0, q0 = col.t[-1] + dt_l + max(prm['dt_deep'], dt_s), col.q[-1] + dq_l + max(prm['dq_deep'], dq_s)
    shallow = trial(col, col.n - 1, t0, q0, ws, 'shallow', rain=False, prm=prm)
    shallow['excess'] = (dt_s, dq_s, ws)
    shallow['surface_layer'] = (dt_l, dq_l)
    deep = None
    for k in range(col.n - 1, -1, -1):
        if col.ps - col.pf[k] > prm['departure_dp']:
            break
        if k < col.n - 1 and prm['elevated_ascent'] and (omega is None or omega[k] >= 0.0):
            continue
        ts, qs_ = source_layer(col, k, prm['source_dp'])
        if k == col.n - 1:                 # the surface parcel takes the larger of the two excesses (IFS (6.19), (6.21))
            ts, qs_ = ts + dt_l + max(prm['dt_deep'], dt_s), qs_ + dq_l + max(prm['dq_deep'], dq_s)
            w0 = surface_start(col, shf, lhf, ws, dt_l, dq_l, prm)
        else:
            ts, qs_, w0 = ts + prm['dt_deep'], qs_ + prm['dq_deep'], prm['w_base_deep']
        tr = trial(col, k, ts, qs_, w0, 'deep', rain=True, prm=prm)
        if tr['ok'] and tr['depth'] > prm['deep_dp']:
            deep = tr
            break
    if deep is not None:
        return dict(type='deep', trial=deep, shallow_trial=shallow, rooted=deep['k'] == col.n - 1)
    if shallow['ok'] and shallow['depth'] <= prm['deep_dp']:
        return dict(type='shallow', trial=shallow, shallow_trial=shallow, rooted=True)
    if omega is not None:
        for k in range(col.n - 1, -1, -1):
            if col.pf[k] < prm['midlevel_top'] * col.ps:
                break
            if omega[k] < 0.0 and col.rh[k] >= prm['midlevel_rh']:
                tr = trial(col, k, col.t[k] + prm['dt_deep'], col.q[k] + prm['dq_deep'], prm['w_base_deep'],
                           'mid', rain=True, prm=prm)
                if tr['ok']:
                    return dict(type='mid', trial=tr, shallow_trial=shallow, rooted=False)
    return dict(type=None, shallow_trial=shallow, rooted=False)


# -------------------------------------------------------------------------------------------- updraught and downdraught

def updraught(col: Column, decision: dict, prm=PARAMETERS) -> dict:
    """The full updraught from the test parcel's cloud base, at unit mass flux there: its points, and on half levels
    its mass flux, liquid-ice static energy and total water (nil below the base and above the top), the rain and snow
    it forms in each layer, its cloud base and top half levels and heights."""
    tr = decision['trial']
    kind = {'deep': 'deep', 'shallow': 'shallow', 'mid': 'deep'}[decision['type']]
    pts = col.path(tr['k'])
    jb = tr['base_point']
    # the path from the base half level up (Column.path lists full and half levels alternately)
    start = next(j for j, p in enumerate(pts) if p[0] == 'h' and p[1] == tr['base_half'])
    b = tr['points'][jb]
    z_cond = tr['points'][tr['base']]['z']                       # the parcel's condensation level
    state = dict(sl=b['sl'], qt=b['qt'], kin=b['kin'], m=1.0, z_base=z_cond)
    points, base, top_z = ascend(col, pts[start:], state, kind, rain=True, prm=prm)
    n = col.n
    m, sl, qt = np.zeros(n + 1), np.zeros(n + 1), np.zeros(n + 1)
    rain, snow = np.zeros(n), np.zeros(n)
    ib = tr['base_half']
    it = ib
    for p in points:
        if p['kind'] == 'h' and p['m'] > 0.0:
            m[p['index']], sl[p['index']], qt[p['index']] = p['m'], p['sl'], p['qt']
            it = p['index']
    for j in range(1, len(points)):
        layer = points[j]['index']
        rain[layer] += points[j]['rain']
        snow[layer] += points[j]['snow']
    cloudy = [p for p in points if p['z'] <= top_z + 1e-6]
    w_mean = float(np.mean([math.sqrt(2.0 * p['kin']) for p in cloudy])) if cloudy else 0.0
    return dict(points=points, m=m, sl=sl, qt=qt, rain=rain, snow=snow, base_half=ib, top_half=it,
                z_base=col.zh[ib], z_top=top_z, p_base=col.ph[ib], kind=kind, w_mean=w_mean)


def wet_bulb(t: float, q: float, z: float, p: float, g: float, prm=PARAMETERS, c=PLASIM):
    """Air brought to saturation by evaporating water into it at constant moist static energy: its temperature and
    vapour."""
    cp, lv = c['acpd'], c['ALV']
    h = cp * t + lv * q
    tw = t
    for _ in range(12):
        qs = qsat(tw, p, c)
        f = cp * tw + lv * qs - h
        step = f / (cp + lv * dqsat_dt(tw, p, c))
        tw -= step
        if abs(step) < 1e-7:
            break
    return tw, qsat(tw, p, c)


def downdraught(col: Column, up: dict, prm=PARAMETERS) -> dict:
    """The saturated downdraught at unit updraught base flux: from its level of free sinking (the highest half level
    above the cloud base, below the level of least moist static energy, where an even mixture of cloud air and
    surrounding air at its wet-bulb temperature is heavier than its surroundings) with -eta_down of the base flux,
    mixing at eps_down, kept saturated by evaporating rain (at most what falls through it), stopping where it turns
    buoyant and spreading out over the lowest down_bottom_dp. Returns its half-level mass flux, dry static energy and
    vapour, and the rain it evaporates in each layer (unit base flux)."""
    c, g, n = col.c, col.g, col.n
    cp, lv = c['acpd'], c['ALV']
    out = dict(m=np.zeros(n + 1), s=np.zeros(n + 1), q=np.zeros(n + 1), evap=np.zeros(n), lfs=None)
    precip = up['rain'] + up['snow']
    if precip.sum() <= 0.0 or up['kind'] == 'shallow':
        return out
    h_env = col.s + lv * col.q
    k_min = int(np.argmin(np.where(col.pf > col.ph[up['top_half']], h_env, np.inf)))
    lfs = None
    for i in range(up['top_half'] + 1, up['base_half']):           # half levels below the top, above the base
        if col.ph[i] < col.pf[k_min] or up['m'][i] <= 0.0:
            continue
        z, p = col.zh[i], col.ph[i]
        tu, qvu, lu = adjust(up['sl'][i], up['qt'][i], z, p, g, prm, c)
        tw, qw = wet_bulb(col.th[i], col.qh[i], z, p, g, prm, c)
        tm, qm, lm = adjust(0.5 * (up['sl'][i] + cp * tw + g * z), 0.5 * (up['qt'][i] + qw), z, p, g, prm, c)
        if virtual(tm, qm, c) - lm * tm < virtual(col.th[i], col.qh[i], c) and precip[:i].sum() > 0.0:
            lfs = i
            sd, qd = cp * tm + g * z, qm
            break
    if lfs is None:
        return out
    out['lfs'] = lfs
    md = -prm['eta_down']
    fall = precip[:lfs].sum()                                     # rain falling through the level of free sinking
    out['m'][lfs], out['s'][lfs], out['q'][lfs] = md, sd, qd
    bottom = col.ps - prm['down_bottom_dp']
    for i in range(lfs + 1, n + 1):
        k = i - 1                                                 # the layer between half levels i - 1 and i
        dz = col.zh[i - 1] - col.zh[i]
        z, p = col.zh[i], col.ph[i]
        te, qe = (col.th[i], col.qh[i]) if i < n else (col.t[-1], col.q[-1])
        sd = mix(sd, cp * col.th[i - 1] + g * col.zh[i - 1], cp * te + g * z, prm['eps_down'], dz)
        qd = mix(qd, col.qh[i - 1], qe, prm['eps_down'], dz)
        td = (sd - g * z) / cp
        tw, qw = wet_bulb(td, qd, z, p, g, prm, c)
        need = max(qw - qd, 0.0) * -md                            # kg/m2/s to saturate the downdraught
        fall += precip[k]
        e = min(need, fall)
        if need > 0.0:
            f = e / need
            td, qd = td + f * (tw - td), qd + f * (qw - qd)
            sd = cp * td + g * z
        fall -= e
        out['evap'][k] += e
        if i == n:
            break
        if virtual(td, qd, c) > virtual(te, qe, c):              # buoyant: it detrains in this layer
            break
        if p > bottom:                                            # spreading out over the lowest down_bottom_dp
            md = -prm['eta_down'] * max(col.ps - p, 0.0) / prm['down_bottom_dp']
        out['m'][i], out['s'][i], out['q'][i] = md, sd, qd
    return out


# --------------------------------------------------------------------------------------------- fluxes and tendencies

def tendencies(col: Column, up: dict, down: dict, scale: float, prm=PARAMETERS, melt=True, t_phase=None,
               t_ground=None) -> dict:
    """The column's tendencies (K/s, kg/kg/s) and precipitation for the updraught and downdraught scaled to a cloud
    base flux of `scale` (kg/m2/s): the convergence of their fluxes of liquid-ice static energy and total water
    against the compensating motion of the surrounding air taken upwind, the fluxes in the sub-cloud layer falling
    linearly in pressure to nil at the ground (IFS 6.5), the rain and snow formed as sources and the downdraught's
    evaporation. The falling water's phase follows `melt`:
    - True: snow melts in the first layer above freezing that it falls into, the scheme taking the heat;
    - False: the falling water keeps the phase it formed in;
    - 'plasim': it takes the phase PlaSim's mkrain gives it (rainmod.f90), all of it rain in a layer whose
      temperature t_phase (PlaSim's at the step's start, the column's by default) is above freezing and snow
      otherwise, and the downdraught evaporates that phase; mkrain itself takes the heat of melting and freezing, so
      the rain and snow handed on stay non-negative in mkrain's running totals however they change phase.
    The heating is converted to temperature with PlaSim's moist heat capacity cp (1 + ADV q). Returns the
    tendencies, the half-level fluxes and mass fluxes, the precipitation flux at each half level and at the ground,
    and the rain and snow each layer hands on net of the downdraught's evaporation."""
    c, g, n = col.c, col.g, col.n
    cp, lv, ls = c['acpd'], c['ALV'], c['ALS']
    mu, md = up['m'] * scale, down['m'] * scale
    s_env, q_env = col.s, col.q
    fs, fq = np.zeros(n + 1), np.zeros(n + 1)
    ib = up['base_half']

    def upwind(i, mnet):
        """The surrounding air crossing half level i: from the level above when the net cloud flux rises (the
        surroundings sink), from the level below when it falls."""
        return (s_env[i - 1], q_env[i - 1]) if mnet > 0.0 else (s_env[i], q_env[i])
    sh_b, qh_b = upwind(ib, mu[ib] + md[ib])
    fu_s, fu_q = mu[ib] * (up['sl'][ib] - sh_b), mu[ib] * (up['qt'][ib] - qh_b)    # the updraught's at the base
    for i in range(1, n):
        if i <= ib:
            mnet = mu[i] + md[i]
            sh, qh = upwind(i, mnet)
            fs[i] = mu[i] * (up['sl'][i] - sh) + md[i] * (down['s'][i] - sh)
            fq[i] = mu[i] * (up['qt'][i] - qh) + md[i] * (down['q'][i] - qh)
        else:                       # below the base the updraught's fluxes fall linearly in pressure to the ground
            zf = (col.ps - col.ph[i]) / (col.ps - col.ph[ib])
            mu[i] = mu[ib] * zf
            sh, qh = upwind(i, mu[i] + md[i])
            fs[i] = zf * fu_s + md[i] * (down['s'][i] - sh)
            fq[i] = zf * fu_q + md[i] * (down['q'][i] - qh)
    rain, snow = up['rain'] * scale, up['snow'] * scale
    evap = down['evap'] * scale
    t_phase = col.t if t_phase is None else np.asarray(t_phase, float)
    dsdt = np.zeros(n)
    dqdt = np.zeros(n)
    p_rain, p_snow = np.zeros(n + 1), np.zeros(n + 1)
    melted = np.zeros(n)
    net_rain, net_snow = np.zeros(n), np.zeros(n)
    for k in range(n):
        heat = lv * rain[k] + ls * snow[k]
        pr, ps_ = p_rain[k] + rain[k], p_snow[k] + snow[k]
        if melt == 'plasim':
            if t_phase[k] > c['TMELT']:
                pr, ps_ = pr + ps_, 0.0
            else:
                pr, ps_ = 0.0, pr + ps_
            e = min(evap[k], pr + ps_)
            e_rain, e_snow = (e, 0.0) if t_phase[k] > c['TMELT'] else (0.0, e)
        else:
            # rain and snow falling into the layer, less what the downdraught evaporates here (snow first, as rain)
            e = min(evap[k], pr + ps_)
            e_snow = min(e, ps_) if col.t[k] < c['TMELT'] else 0.0
            e_rain = e - e_snow
            if e_rain > pr:
                e_snow += e_rain - pr
                e_rain = pr
        pr -= e_rain
        ps_ -= e_snow
        heat -= lv * e_rain + ls * e_snow
        net_rain[k], net_snow[k] = rain[k] - e_rain, snow[k] - e_snow
        if melt is True and ps_ > 0.0 and col.t[k] > c['TMELT']:
            melted[k] = ps_
            heat -= (ls - lv) * ps_
            net_rain[k] += ps_
            net_snow[k] -= ps_
            pr += ps_
            ps_ = 0.0
        p_rain[k + 1], p_snow[k + 1] = pr, ps_
        dsdt[k] = g * (fs[k + 1] - fs[k]) / col.dp[k] + g * heat / col.dp[k]
        dqdt[k] = g * (fq[k + 1] - fq[k]) / col.dp[k] - g * (rain[k] + snow[k] - e_rain - e_snow) / col.dp[k]
    return dict(dtdt=dsdt / (cp * (1.0 + c['ADV'] * col.q)), dqdt=dqdt, fs=fs, fq=fq, mu=mu, md=md, p_rain=p_rain,
                p_snow=p_snow, rain_formed=rain, snow_formed=snow, evap_down=evap, melt=melted, net_rain=net_rain,
                net_snow=net_snow, rain_surface=float(p_rain[n]), snow_surface=float(p_snow[n]))


# --------------------------------------------------------------------------------------------------------- closures

def pcape(col: Column, up: dict) -> float:
    """The updraught's density-weighted CAPE (Pa): its buoyancy (T_v,u - T_v)/T_v less its condensate, counted where
    positive, by pressure from its base to its top (Bechtold and others 2014, (2))."""
    total = 0.0
    pts = up['points']
    for j in range(1, len(pts)):
        a, b = pts[j - 1], pts[j]
        ba = a['b'] / col.g
        bb = b['b'] / col.g
        total += 0.5 * (max(ba, 0.0) + max(bb, 0.0)) * (a['p'] - b['p'])
    return total


def pcape_rate(col: Column, up: dict, tend: dict) -> float:
    """The rate (Pa/s) at which tendencies reduce the updraught's PCAPE, the updraught held as it is: the
    environment's virtual-temperature tendency over its virtual temperature, by pressure, over the layers where the
    updraught is buoyant (Bechtold and others 2014, (7), with the scheme's own tendencies)."""
    c = col.c
    eps1 = 1.0 / c['rdbrv'] - 1.0
    buoyant = np.zeros(col.n, bool)
    for p in up['points']:
        if p['kind'] == 'f' and p['b'] > 0.0:
            buoyant[p['index']] = True
    dtv = tend['dtdt'] * (1.0 + eps1 * col.q) + eps1 * col.t * tend['dqdt']
    return float(np.sum(np.where(buoyant, dtv / col.tv * col.dp, 0.0)))


def boundary_pcape(col: Column, up: dict, bl_dtdt, bl_dqdt, tau_bl: float, prm=PARAMETERS) -> float:
    """PCAPEbl (Pa): tau_bl / T* times the sub-cloud layer's virtual-temperature tendency from processes other than
    convection, by pressure from the ground to the cloud base (Bechtold and others 2014, (17))."""
    c = col.c
    eps1 = 1.0 / c['rdbrv'] - 1.0
    total = 0.0
    for k in range(col.n - 1, -1, -1):
        if col.ph[k + 1] <= up['p_base']:
            break
        dp = col.ph[k + 1] - max(col.ph[k], up['p_base'])
        dtv = bl_dtdt[k] * (1.0 + eps1 * col.q[k]) + eps1 * col.t[k] * bl_dqdt[k]
        total += dtv * dp
    t_star = prm['t_star']
    if t_star == 'cloud':
        t_star = col.g * max(up['z_top'] - up['z_base'], 1.0) / col.c['acpd']
    return tau_bl / t_star * total


def sub_cloud_wind(col: Column, up: dict, u_bl) -> float:
    """The sub-cloud layer's wind speed (m/s): u_bl itself, or for wind speeds by level their mean by pressure from the
    ground to the cloud base."""
    if np.ndim(u_bl) == 0:
        return float(u_bl)
    total, mass = 0.0, 0.0
    for k in range(col.n - 1, -1, -1):
        if col.ph[k + 1] <= up['p_base']:
            break
        dp = col.ph[k + 1] - max(col.ph[k], up['p_base'])
        total += u_bl[k] * dp
        mass += dp
    return total / mass if mass > 0.0 else float(u_bl[-1])


def boundary_mse_supply(col: Column, up: dict, bl_dtdt, bl_dqdt) -> float:
    """The moist static energy (W/m2) that processes other than convection bring into the sub-cloud layer, from the
    ground to the cloud base."""
    c = col.c
    total = 0.0
    for k in range(col.n - 1, -1, -1):
        if col.ph[k + 1] <= up['p_base']:
            break
        dp = col.ph[k + 1] - max(col.ph[k], up['p_base'])
        total += (c['acpd'] * bl_dtdt[k] + c['ALV'] * bl_dqdt[k]) * dp / col.g
    return total


# ----------------------------------------------------------------------------------------------------------- column

def column(t, q, p_full, p_half, g: float, dt: float = 1800.0, shf: float = 0.0, lhf: float = 0.0, bl_dtdt=None,
           bl_dqdt=None, omega=None, u_bl=0.0, land: bool = False, prm=PARAMETERS, c=PLASIM, z_full=None,
           z_half=None, melt=True, t_phase=None) -> dict:
    """The scheme on one column (arrays top first): temperature (K), specific humidity, full and half level pressures
    (Pa; half level 0 the top, n the ground), gravity, PlaSim's time step dt (s), the surface's sensible and latent heat
    fluxes (W/m2, upward), the tendencies of temperature and vapour from processes other than convection (K/s,
    kg/kg/s by level: the boundary layer's turbulence, radiation and advection), the large-scale vertical velocity
    (Pa/s by level, for mid-level convection), the sub-cloud wind speed (m/s; or the wind speed by level, whose mean by
    mass from the ground to the cloud base is taken) and whether the ground is land; the falling water's phase follows
    melt and t_phase (tendencies). Returns
    the type, the cloud's base and top, its base mass flux and closure terms, the tendencies (K/s, kg/kg/s), the
    half-level mass fluxes (kg/m2/s), the rain and snow formed and reaching the ground before PlaSim's evaporation of
    falling rain (kg/m2/s), and the test parcels' results."""
    col = Column(t, q, p_full, p_half, g, c, z_full, z_half)
    n = col.n
    bl_dtdt = np.zeros(n) if bl_dtdt is None else np.asarray(bl_dtdt, float)
    bl_dqdt = np.zeros(n) if bl_dqdt is None else np.asarray(bl_dqdt, float)
    decision = classify(col, shf, lhf, omega, prm)
    out = dict(type=decision['type'], rooted=decision['rooted'], dtdt=np.zeros(n), dqdt=np.zeros(n),
               mu=np.zeros(n + 1), md=np.zeros(n + 1), rain_formed=np.zeros(n), snow_formed=np.zeros(n),
               net_rain=np.zeros(n), net_snow=np.zeros(n),
               evap_down=np.zeros(n), rain_surface=0.0, snow_surface=0.0, m_base=0.0, z_base=None, z_top=None,
               p_base=None, p_top=None, pcape=0.0, pcape_bl=0.0, tau=None, tau_bl=None, supply=None, cfl_limited=False,
               shallow_trial=summary(decision['shallow_trial']), trial=summary(decision.get('trial')), column=col)
    if decision['type'] is None:
        return out
    up = updraught(col, decision, prm)
    down = downdraught(col, up, prm)
    unit = tendencies(col, up, down, 1.0, prm, melt, t_phase)
    kind = decision['type']
    pc = pcape(col, up)
    out.update(z_base=up['z_base'], z_top=up['z_top'], p_base=up['p_base'], pcape=pc, up=up, down=down,
               p_top=top_pressure(up['points'], up['z_top']))
    if kind == 'deep':
        tau_c = max(up['z_top'] - up['z_base'], 1.0) / max(up['w_mean'], 0.1)
        tau = max(prm['tau_factor'] * tau_c, prm['tau_min'], dt)
        if decision['rooted']:
            tau_bl = tau_c if land else max(up['z_base'], 1.0) / max(sub_cloud_wind(col, up, u_bl), prm['u_bl_min'])
            pbl = boundary_pcape(col, up, bl_dtdt, bl_dqdt, tau_bl, prm)
        else:
            tau_bl, pbl = None, 0.0
        rate = pcape_rate(col, up, unit)
        m_base = max(0.0, (pc - pbl) / tau / rate) if rate > 0.0 else 0.0
        out.update(tau=tau, tau_c=tau_c, tau_bl=tau_bl, pcape_bl=pbl, pcape_rate_unit=rate)
    elif kind == 'shallow':
        supply = boundary_mse_supply(col, up, bl_dtdt, bl_dqdt)
        ib = up['base_half']
        h_up = up['sl'][ib] + c['ALV'] * up['qt'][ib]
        h_env = c['acpd'] * col.th[ib] + g * col.zh[ib] + c['ALV'] * col.qh[ib]
        excess = h_up - h_env
        m_base = max(0.0, supply / excess) if excess > 0.0 else 0.0
        out.update(supply=supply, base_excess=excess)
    else:                                                         # mid-level: the large-scale mass flux at the base
        k = decision['trial']['k']
        m_base = max(0.0, -omega[k] / g)
    # the mass flux's limit: its net flux may cross no layer faster than cfl times that layer's mass a step
    mnet = np.abs(unit['mu'] + unit['md'])
    limit = np.inf
    for i in range(1, n):
        if mnet[i] > 0.0:
            limit = min(limit, prm['cfl'] * min(col.mass[i - 1], col.mass[i]) / dt / mnet[i])
    out['cfl_limited'] = m_base > limit
    m_base = min(m_base, limit)
    tend = tendencies(col, up, down, m_base, prm, melt, t_phase)
    out.update(m_base=m_base, dtdt=tend['dtdt'], dqdt=tend['dqdt'], mu=tend['mu'], md=tend['md'],
               rain_formed=tend['rain_formed'], snow_formed=tend['snow_formed'], evap_down=tend['evap_down'],
               rain_surface=tend['rain_surface'], snow_surface=tend['snow_surface'], fs=tend['fs'], fq=tend['fq'],
               p_rain=tend['p_rain'], p_snow=tend['p_snow'], melt=tend['melt'], net_rain=tend['net_rain'],
               net_snow=tend['net_snow'])
    return out


def summary(tr):
    """A test parcel's result without its points."""
    if tr is None:
        return None
    keep = ('k', 'ok', 'top_z', 'p_base', 'p_top', 'depth', 'base_half', 'excess')
    return {k: tr[k] for k in keep if k in tr}


# ================================================================================================= offline tests

DAY = 86400.0
STEP_S = 1800.0                   # PlaSim's time step at the runs' T21 (deltsec)
QUARTERS = (('night', -180.0, -90.0), ('morning', -90.0, 0.0), ('afternoon', 0.0, 90.0), ('evening', 90.0, 180.0))
BOXES = (('sea', 'box_near_side_sea_day3', (59.0, 88.5)), ('land', 'box_0e_design', (29.5, 59.0)))
BANDS_M = ((0.0, 6000.0), (6000.0, 10000.0), (10000.0, 16000.0), (16000.0, 35000.0), (35000.0, 48000.0),
           (48000.0, 63000.0), (63000.0, 150000.0))
GCM_RUN, GCM_YEARS = 'A28R_301k_cold_sbm_rh07', (10, 14)


def quarter_name(hour_deg: float) -> str:
    h = (float(hour_deg) + 180.0) % 360.0 - 180.0
    return next(name for name, a, b in QUARTERS if a <= h < b)


def omega_equivalent(heating, t, zh, rho, g) -> np.ndarray:
    """The large-scale vertical velocity (Pa/s, ground up) a CM1 box's forcing stands for at each level: the ascent
    whose adiabatic cooling gives the forcing's heating (the mean ascent acting on the box's air and the day-night
    table, which holds the GCM's circulation by local time), w = -heating / (dT/dz + g / cp), the stability floored at
    a tenth of g / cp."""
    from climate.crm import cm1_output as co
    stab = np.maximum(np.gradient(np.asarray(t, float), zh) + g / co.CP, 0.1 * g / co.CP)
    return np.asarray(rho, float) * g * np.asarray(heating, float) / stab


def cm1_outputs(case: str, span) -> dict:
    """A CM1 box's domain-mean columns at each diagnostic output over the span, as the scheme takes them (top first):
    temperature, specific humidity, full and half level pressures and heights; the tendencies of processes other than
    convection (radiation, the prescribed day-night table and the mean ascent acting on the box's air, and the surface
    fluxes laid into the lowest level) and CM1's own convective tendencies (microphysics, resolved transport and the
    boundary-layer scheme with its dissipation), K/s and kg/kg/s; the surface fluxes (W/m2), the mean ascent as omega
    (Pa/s), the sub-cloud wind speed, CM1's rain (mm/day), CAPE and inhibition, and each output's hour angle."""
    import json as _json
    from climate.crm import box_budget as bb
    from climate.crm import cm1_output as co
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    from climate.gcm import sbm_offline as so
    folder = ra.RUNS / case
    record = _json.loads((folder / 'case.json').read_text())
    cfg = record['configuration']
    geo = ra.case_geometry(folder)
    zh, zw = geo['zh'], geo['zw']
    dz = np.diff(zw)
    layout = ra.grads_layout(folder / 'cm1out_diag_s.ctl')[1]
    text = (folder / 'namelist.input').read_text()
    var18, var19 = bb.namelist_value(text, 'var18'), bb.namelist_value(text, 'var19')
    z_tab, th_tab = bb.read_day_night(folder / 'terluna_lsadv.txt')
    _, qv_tab = bb.read_day_night(folder / 'terluna_lsadv.txt', 'qv')
    wls = np.loadtxt(folder / 'terluna_wls.txt', skiprows=1)
    w_mean = np.interp(zh, wls[:, 0], wls[:, 1])
    g = record['build']['gravity_m_s2']
    rows = []
    for path in sorted(folder.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(folder, path, layout)
        if not span[0] * DAY + 1.0 < d['mtime'] <= span[1] * DAY + 1.0:
            continue
        ex = (d['prs'] / co.P00) ** (co.RD / co.CP)
        rho = d['rho']
        group = lambda names: sum(d[n] for n in names if n in d)
        heat_forcing = (d['ptb_lsw'] + bb.day_night_heating(z_tab, th_tab, zh, d['mtime'], var18, var19)) * ex
        heat = d['ptb_rad'] * ex + heat_forcing
        moist = d['qvb_lsw'] + bb.day_night_heating(z_tab, qv_tab, zh, d['mtime'], var18, var19)
        heat = heat.copy()
        moist = moist.copy()
        heat[0] += d['hfx'] / (co.CP * rho[0] * dz[0])
        moist[0] += d['qfx'] / (rho[0] * dz[0])
        table_t = bb.day_night_heating(z_tab, th_tab, zh, d['mtime'], var18, var19) * ex
        table_q = bb.day_night_heating(z_tab, qv_tab, zh, d['mtime'], var18, var19)
        conv_t = (d['ptb_mp'] + group(bb.GROUPS['resolved_transport']) + group(bb.GROUPS['turbulence'])) * ex
        conv_q = d['qvb_mp'] + group(bb.WATER['own_transport']) + group(bb.WATER['boundary_layer'])
        pw = so.w_level_pressures(d['prs'], d['psfc'], zh, zw)
        hour = local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], d['mtime'], solar_day_s())
        below = zh < max(d['hpbl'], 1000.0)
        rows.append(dict(
            day=d['mtime'] / DAY, hour=hour, quarter=quarter_name(hour), t=d['t'][::-1], q=d['shv'][::-1],
            p_full=d['prs'][::-1], p_half=pw[::-1], z_full=zh[::-1], z_half=zw[::-1], bl_dtdt=heat[::-1],
            bl_dqdt=moist[::-1], conv_dtdt=conv_t[::-1], conv_dqdt=conv_q[::-1], rad_dtdt=(d['ptb_rad'] * ex)[::-1],
            omega=omega_equivalent(heat_forcing, d['t'], zh, rho, g)[::-1], shf=float(d['hfx']),
            lhf=float(d['qfx']) * bb.XLV,
            u_bl=float(np.mean(d['wsp'][below])), rain=float(d['prate']) * DAY, cape=float(d['cape']),
            table_dtdt=table_t[::-1], table_dqdt=table_q[::-1], w_mean=w_mean[::-1],
            cin=float(d['cin']), rho=rho[::-1], g=g, mass=(rho * dz)[::-1], ust=float(d['ust'])))
    return dict(case=case, span=list(span), g=g, rows=rows, zh=zh, zw=zw, land=cfg['site']['surface'] == 'land',
                inputs={'runner': record.get('runner'), 'build': record['build']['label']})


def cm1_microphysics(case: str, span) -> dict:
    """CM1's running water totals (cm1out_stats) by quarter: the water condensed, evaporated back from cloud,
    evaporated from rain, snow and graupel on the way down, and reaching the ground (mm/day), and the falling water
    formed, the last two together."""
    import json as _json
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    folder = ra.RUNS / case
    ctl = (folder / 'cm1out_stats.ctl').read_text().splitlines()
    i = next(k for k, line in enumerate(ctl) if line.lower().startswith('vars'))
    names = [line.split()[0] for line in ctl[i + 1:i + 1 + int(ctl[i].split()[1])]]
    tot = dict(zip(names, np.fromfile(folder / 'cm1out_stats.dat', dtype='<f4').reshape(-1, len(names)).T))
    record = _json.loads((folder / 'case.json').read_text())
    grid, cfg = record['grid'], record['configuration']
    area = grid['nx'] * grid.get('ny', 1) * grid['dx_m'] ** 2
    t = tot['mtime'].astype(float)
    keep = (t > span[0] * DAY) & (t <= span[1] * DAY + 1.0)
    idx = np.flatnonzero(keep)
    out = {q[0]: dict(seconds=0.0, condensed=0.0, from_cloud=0.0, on_the_way_down=0.0, ground=0.0) for q in QUARTERS}
    for a, b in zip(idx[:-1], idx[1:]):
        mid = 0.5 * (t[a] + t[b])
        name = quarter_name(local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], mid, solar_day_s()))
        o = out[name]
        o['seconds'] += t[b] - t[a]
        for key, var in (('condensed', 'tcond'), ('from_cloud', 'tevac'), ('on_the_way_down', 'tevar'),
                         ('ground', 'train')):
            o[key] += float(tot[var][b] - tot[var][a]) / area
    res = {}
    for name, o in out.items():
        s = o.pop('seconds')
        res[name] = {f'{k}_mm_day': v / s * DAY for k, v in o.items()}
        res[name]['falling_formed_mm_day'] = res[name]['on_the_way_down_mm_day'] + res[name]['ground_mm_day']
    return res


def evaporate(col: Column, rain, snow, dtdt, dqdt, gamma: float, step2: float = 2.0 * STEP_S, t_old=None,
              t_ground=None) -> dict:
    """PlaSim's mkrain (rainmod.f90) on the convective rain and snow the scheme hands it by layer (kg/m2/s, top first,
    net of its downdraught's evaporation), as the mass-flux patch leaves it for nprc 3 (massflux_fortran.py), in
    mkrain's own arithmetic and units (each layer's water a step times its share of sigma). Layer by layer from the
    top the handed water joins the running rain and snow; all of it turns to snow where the layer's temperature at the
    step's start (t_old, the column's by default) is at or below freezing and to rain above, the latent heat of the
    change going into that layer, and in the lowest layer the ground's temperature t_ground converts it again. Each
    unsaturated layer, its state the column's plus the tendencies dtdt and dqdt over the step (step2, PlaSim's
    deltsec2), then evaporates from the falling rain and snow a share gamma of its saturation deficit over one plus
    L/cp dqs/dT, cp PlaSim's moist heat capacity, at most what falls; and under nprc 3 at most what the water handed on
    below still needs, the floor: the least running total of the handed water from that layer down, less what mkrain
    has evaporated above. So the downdraught's evaporation, already taken from the handed water, never leaves the
    running totals negative. Returns the rain and snow reaching the ground (kg/m2/s) and mkrain's heating and
    moistening (K/s, kg/kg/s)."""
    c = col.c
    n = col.n
    ps = col.ps
    dsig = col.dp / ps
    t_old = col.t if t_old is None else np.asarray(t_old, float)
    to_units = col.g / ps                        # kg/m2/s to mkrain's dq/dt dsigma
    prcl = (np.asarray(rain, float) + np.asarray(snow, float)) * to_units
    prscl = np.asarray(snow, float) * to_units
    running = np.cumsum(prcl)
    floor = np.minimum.accumulate(running[::-1])[::-1]
    taken = 0.0
    zprc, zprsc = 0.0, 0.0
    heat, wet = np.zeros(n), np.zeros(n)
    for k in range(n):
        zprc += prcl[k] - prscl[k]
        zprsc += prscl[k]
        zlcp = (c['ALS'] - c['ALV']) / (c['acpd'] * (1.0 + c['ADV'] * col.q[k]))
        for temp, add in ((t_old[k], False),) + (((t_ground, True),) if k == n - 1 and t_ground is not None else ()):
            if zprsc > 0.0 and temp > c['TMELT']:
                heat[k] = (heat[k] if add else 0.0) - zprsc / dsig[k] * zlcp
                zprc, zprsc = zprsc + zprc, 0.0
            elif zprc > 0.0 and temp <= c['TMELT']:
                heat[k] = (heat[k] if add else 0.0) + zprc / dsig[k] * zlcp
                zprsc, zprc = zprc + zprsc, 0.0
        zq = col.q[k] + dqdt[k] * step2
        zt = t_old[k] + dtdt[k] * step2
        zqsat = c['rdbrv'] * c['ra1'] * math.exp(c['ra2'] * (zt - c['TMELT']) / (zt - c['ra4'])) / col.pf[k]
        zqsat = min(zqsat, c['rdbrv'])
        zcor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * zqsat)
        zqsat = max(zqsat * zcor, 0.0)
        remaining = max(floor[k] - taken, 0.0)
        for phase in ('rain', 'snow'):
            avail = zprc if phase == 'rain' else zprsc
            if not (zq < zqsat and avail > 0.0):
                continue
            lat = c['ALV'] if phase == 'rain' else c['ALS']
            zlcpe = lat / (c['acpd'] * (1.0 + c['ADV'] * col.q[k]))
            e = min(gamma * (zqsat - zq) * dsig[k] / step2 / (1.0 + zlcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat
                                                               * zcor / (zt - c['ra4']) ** 2), avail, remaining)
            remaining -= e
            taken += e
            if phase == 'rain':
                zprc -= e
            else:
                zprsc -= e
            heat[k] -= e / dsig[k] * zlcpe
            wet[k] += e / dsig[k]
    return dict(rain=zprc / to_units, snow=zprsc / to_units, dtdt=heat, dqdt=wet)


def cm1_samples(case: str, span) -> dict:
    """The CM1 box's clouds by quarter from its 3-D output (box_entrainment's samples: cloud, cloudy updraught and
    core), to the column's top: each sample's mean mass flux (kg/m2/s over the domain) by height, and the tops of the
    columns' lowest unbroken cloud whose base lies below 16 km (m)."""
    from climate.crm import box_entrainment as be
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    info = be.case_info(case)
    cfg = info['record']['configuration']
    zh = info['zh']
    acc = {q[0]: dict(n=0, columns=0, updraft=np.zeros(zh.size), core=np.zeros(zh.size), tops=[]) for q in QUARTERS}
    for n in sorted(int(f.name[8:14]) for f in info['case'].glob('cm1out_t*_s.dat')):
        t = (n - 1) * info['tap']
        if not span[0] * DAY + 1.0 < t <= span[1] * DAY + 1.0:
            continue
        x = be.thermo(ra.read_snapshot(info['case'], n), zh, info['g'])
        m = be.masks(x)
        a = acc[quarter_name(local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], t, solar_day_s()))]
        cols = x['t'].shape[1]
        a['n'] += 1
        a['columns'] += cols
        for key in ('updraft', 'core'):
            a[key] += (x['rho'] * x['w'] * m[key]).sum(axis=1) / cols
        a['tops'].append(be.cloud_tops(x['qcl'], zh, limit_m=be.LOWER_M))
    out = {}
    for name, a in acc.items():
        tops = np.concatenate(a['tops']) if a['tops'] else np.zeros(0)
        out[name] = dict(snapshots=a['n'], zh=zh, updraft=a['updraft'] / max(a['n'], 1), core=a['core'] / max(a['n'], 1),
                         cloudy_share=float(tops.size / max(a['columns'], 1)),
                         top_median_m=float(np.median(tops)) if tops.size else None,
                         top_p90_m=float(np.percentile(tops, 90)) if tops.size else None,
                         top_above_16km_share=float(np.mean(tops > 16000.0)) if tops.size else None)
    return out


FIFTEEN_LOWEST_KM = 0.3          # the provisional fifteen levels' lowest full level, until levels.py holds a set


def grid_sets() -> dict:
    """The half levels (sigma, the top layer's bottom first, ending at 1) of PlaSim's ten levels, the fourteen of
    levels.py's l14, and fifteen: levels.py's fifteen-level set when it holds one, otherwise a provisional set, l14
    with its lowest layer split by a half level (sigma to four places) that puts the lowest full level near
    FIFTEEN_LOWEST_KM at levels.py's scale height."""
    from climate.gcm import levels
    out = {'10': np.asarray(levels.half_levels('10'), float), '14': np.asarray(levels.half_levels('l14'), float)}
    fifteen = sorted(k for k, v in levels.LEVEL_SETS.items() if len(v) == 15)
    if fifteen:
        out['15'] = np.asarray(levels.half_levels(fifteen[0]), float)
        out['source'] = f'climate/gcm/levels.py (sets 10, l14 and {fifteen[0]})'
    else:
        split = round(2.0 * math.exp(-FIFTEEN_LOWEST_KM / levels.SCALE_HEIGHT_KM) - 1.0, 4)
        out['15'] = np.r_[out['14'][:-1], split, 1.0]
        out['source'] = ('climate/gcm/levels.py (sets 10 and l14); fifteen provisional, l14 with a half level at '
                         f'sigma {split}')
    return out


def to_grid(x: dict, half_sigma) -> dict:
    """A CM1 output laid on a PlaSim grid at its surface pressure: each layer takes the mean by mass of CM1's
    temperature, vapour, forcing and convective tendencies, and large-scale vertical velocity over its pressure range,
    the full levels midway in sigma between the half levels; heights from the layers' virtual temperatures."""
    ph_cm1 = np.asarray(x['p_half'], float)                    # top first: CM1's top ... the ground
    ps = float(ph_cm1[-1])
    half = np.concatenate([[0.0], np.asarray(half_sigma, float)])
    target = half * ps
    full = 0.5 * (half[:-1] + half[1:]) * ps

    def mean(values):
        v = np.asarray(values, float)
        held = np.concatenate([[0.0], np.cumsum(v * np.diff(ph_cm1))])
        f = np.interp(target, ph_cm1, held)
        covered = target[1:] - np.maximum(target[:-1], ph_cm1[0])
        out = (f[1:] - f[:-1]) / np.where(covered > 0.0, covered, 1.0)
        top = covered <= 0.0
        out[top] = v[0]                                         # a layer wholly above CM1's top takes its top value
        return out
    out = {k: mean(x[k]) for k in ('t', 'q', 'bl_dtdt', 'bl_dqdt', 'conv_dtdt', 'conv_dqdt', 'omega', 'rad_dtdt')}
    out.update(p_full=full, p_half=np.concatenate([[0.0], target[1:]]), z_full=None, z_half=None, g=x['g'],
               shf=x['shf'], lhf=x['lhf'], u_bl=x['u_bl'])
    return out


def band_sums(dtdt, dqdt, ph, z_of_p) -> dict:
    """Heating (W/m2) and moistening (mm/day) of the bands BANDS_M from tendencies held as layer values between the
    half-level pressures ph (top first), the bands' bounds placed in pressure by z_of_p (heights to pressure)."""
    c = PLASIM
    out = {}
    for a, b in BANDS_M:
        pa, pb = float(z_of_p(a)), float(z_of_p(b))            # pressure at the band's bottom and top
        heat, wet = 0.0, 0.0
        for k in range(len(dtdt)):
            lo, hi = max(ph[k], pb), min(ph[k + 1], pa)
            if hi > lo:
                heat += dtdt[k] * (hi - lo)
                wet += dqdt[k] * (hi - lo)
        key = f'{a / 1000:g}-{b / 1000:g} km'
        out[key] = dict(heat_w_m2=heat * c['acpd'] / PLASIM_G, moist_mm_day=wet / PLASIM_G * DAY)
    return out


PLASIM_G = 1.6242                 # the bands' masses use CM1's gravity, the box's own (case.json build)
PROFILE_Z_M = (1000.0, 2000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 12000.0, 14000.0, 16000.0,
               20000.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, 55000.0, 60000.0, 65000.0, 70000.0)
GAMMA = dict(sea=0.027, land=0.046, plasim=0.01)   # CM1 boxes' shares on PlaSim's pooled runs (rain_evaporation.json)


def one_output(x: dict, land: bool, gamma: dict, prm=PARAMETERS, c=PLASIM) -> dict:
    """The scheme on one column with PlaSim's evaporation of its falling rain: the type, the cloud, the closure's
    terms, the mass fluxes and the tendencies of convection and evaporation together, the precipitation formed, at
    the ground before evaporation and after it at each gamma, all on the column's own levels."""
    r = column(x['t'], x['q'], x['p_full'], x['p_half'], x['g'], STEP_S, x['shf'], x['lhf'], x['bl_dtdt'],
               x['bl_dqdt'], x['omega'], x['u_bl'], land, prm=prm, c=c, z_full=x['z_full'], z_half=x['z_half'],
               melt='plasim')
    col = r['column']
    ev = {k: evaporate(col, r['net_rain'], r['net_snow'], r['dtdt'], r['dqdt'], g_) for k, g_ in gamma.items()}
    main = ev['site']
    return dict(type=r['type'], rooted=r['rooted'], z_base=r['z_base'], z_top=r['z_top'], m_base=r['m_base'],
                pcape=r['pcape'], pcape_bl=r['pcape_bl'], tau=r['tau'], supply=r['supply'],
                zf=col.zf, zh=col.zh, ph=col.ph, mu=r['mu'], md=r['md'],
                dtdt=r['dtdt'] + main['dtdt'], dqdt=r['dqdt'] + main['dqdt'],
                formed=float((r['rain_formed'] + r['snow_formed']).sum()) * DAY,
                handed=float((r['net_rain'] + r['net_snow']).sum()) * DAY,
                ground={k: (e['rain'] + e['snow']) * DAY for k, e in ev.items()})


def quarter_summary(results: list, rows: list, zh_cm1, mass_cm1) -> dict:
    """Means and shares over one quarter's outputs of the scheme's results, beside CM1's at the same outputs."""
    n = len(results)
    conv = [r for r in results if r['type'] in ('deep', 'shallow', 'mid')]
    deep = [r for r in results if r['type'] == 'deep']
    shallow = [r for r in results if r['type'] == 'shallow']
    prof = lambda key, z_key: np.mean([np.interp(PROFILE_Z_M, r[z_key][::-1], r[key][::-1]) for r in results], axis=0)
    out = dict(outputs=n,
               type_share={k: sum(1 for r in results if r['type'] == k) / n for k in ('deep', 'shallow', 'mid', None)},
               deep_rooted_share=float(np.mean([r['rooted'] for r in deep])) if deep else None,
               cloud={kind: dict(base_median_m=float(np.median([r['z_base'] for r in lst])),
                                 top_median_m=float(np.median([r['z_top'] for r in lst])),
                                 top_p90_m=float(np.percentile([r['z_top'] for r in lst], 90)),
                                 base_mass_flux_median=float(np.median([r['m_base'] for r in lst])))
                      for kind, lst in (('deep', deep), ('shallow', shallow)) if lst},
               base_mass_flux_mean=float(np.mean([r['m_base'] for r in results])),
               updraught_mass_flux=prof('mu', 'zh').tolist(), downdraught_mass_flux=prof('md', 'zh').tolist(),
               heating_k_day=(prof('dtdt', 'zf') * DAY).tolist(), moistening_g_kg_day=(prof('dqdt', 'zf') * DAY * 1e3).tolist(),
               precipitation_mm_day=dict(formed=float(np.mean([r['formed'] for r in results])),
                                         handed_on=float(np.mean([r['handed'] for r in results])),
                                         **{f'ground_gamma_{k}': float(np.mean([r['ground'][k] for r in results]))
                                            for k in results[0]['ground']}),
               pcape_median_pa=float(np.median([r['pcape'] for r in deep])) if deep else None,
               tau_median_h=float(np.median([r['tau'] for r in deep])) / 3600.0 if deep else None)
    # the bands, the scheme's on its own layers and CM1's on its own, both placed in pressure by CM1's column
    bands_s, bands_c = [], []
    for r, x in zip(results, rows):
        zf_cm1 = x['z_full'][::-1]
        lnp = np.log(x['p_full'][::-1])
        z_of_p = lambda z: math.exp(float(np.interp(z, zf_cm1, lnp))) if z > zf_cm1[0] else float(x['p_half'][-1])
        bands_s.append(band_sums(r['dtdt'], r['dqdt'], r['ph'], z_of_p))
        bands_c.append(band_sums(x['conv_dtdt'], x['conv_dqdt'], x['p_half'], z_of_p))
    keys = list(bands_s[0])
    out['bands'] = {k: dict(scheme_heat_w_m2=float(np.mean([b[k]['heat_w_m2'] for b in bands_s])),
                            cm1_heat_w_m2=float(np.mean([b[k]['heat_w_m2'] for b in bands_c])),
                            scheme_moist_mm_day=float(np.mean([b[k]['moist_mm_day'] for b in bands_s])),
                            cm1_moist_mm_day=float(np.mean([b[k]['moist_mm_day'] for b in bands_c]))) for k in keys}
    cm1_prof = lambda key: np.mean([np.interp(PROFILE_Z_M, x['z_full'][::-1], x[key][::-1]) for x in rows], axis=0)
    out['cm1'] = dict(rain_mm_day=float(np.mean([x['rain'] for x in rows])),
                      heating_k_day=(cm1_prof('conv_dtdt') * DAY).tolist(),
                      moistening_g_kg_day=(cm1_prof('conv_dqdt') * DAY * 1e3).tolist(),
                      cape_j_kg=float(np.mean([x['cape'] for x in rows])), cin_j_kg=float(np.mean([x['cin'] for x in rows])))
    return out


def run_box(tag: str, box: dict, grid=None, prm=PARAMETERS, c=PLASIM) -> dict:
    """The scheme on every output of a CM1 box, on CM1's own levels (grid None) or laid on a PlaSim grid (half levels
    in sigma), summarised by quarter."""
    rows = box['outputs']['rows']
    gamma = dict(site=GAMMA[tag], plasim=GAMMA['plasim'])
    by_q = {q[0]: ([], []) for q in QUARTERS}
    for x in rows:
        xg = x if grid is None else to_grid(x, grid)
        r = one_output(xg, box['outputs']['land'], gamma, prm, c)
        by_q[x['quarter']][0].append(r)
        by_q[x['quarter']][1].append(x)
    return {name: quarter_summary(res, rws, box['outputs']['zh'], None) for name, (res, rws) in by_q.items()}


def gcm_site_columns(run: str, years, box: str) -> dict:
    """The GCM's 3-day mean columns in the cells of a CM1 box's site (its case's gcm_site cells), as the scheme takes
    them: temperature, vapour and surface pressure; the surface's sensible and latent heat (W/m2, upward); the
    tendencies of processes other than convection, the surface fluxes laid into the lowest layer with the radiation's
    and the dynamics' moistening at the mean's last step (snapshots; the dynamics' heating is not written); omega from
    PlaSim's continuity on the means' winds (gcm_site.omega); the 10 m wind; the GCM's convective and total rain; the
    site's hour angle; cos-latitude weights."""
    import json as _json
    import netCDF4
    from climate.crm import gcm_site as gs
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import planet
    from climate.gcm import sbm_offline as so
    record = _json.loads((ra.RUNS / box / 'case.json').read_text())
    cells = [tuple(cc) for cc in record['gcm_site']['cells']]
    lon_site = record['configuration']['site']['lon_deg']
    c = so.run_constants(run)
    pl = planet()
    out = {k: [] for k in ('t', 'q', 'ps', 'shf', 'lhf', 'heat', 'moist', 'omega', 'wind', 'prc', 'pr', 'hour', 'w',
                           'year')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            f = {k: np.asarray(d[k][:], float) for k in ('ta', 'hus', 'ps', 'hfss', 'hfls', 'prc', 'pr', 'dt_lw',
                                                         'dt_sw', 'dq_dyn', 'ua', 'va', 'uas', 'vas', 'czen')}
        edges = gs.plasim_edges(sigma)
        sp, wind = gs.Spectral(lat.size, lon.size), gs.Spectral(lat.size, lon.size, gs.TRUNCATION + 1)
        om = gs.omega(f['ua'], f['va'], f['ps'] * 100.0, sigma, pl['radius_m'], sp, wind)
        idx = [(int(np.argmin(abs(lat - a))), int(np.argmin(abs(lon - b)))) for a, b in cells]
        hour = (lon_site - gs.subsolar_longitude(f['czen'], lat, lon) + 180.0) % 360.0 - 180.0
        dsig = np.diff(edges)
        for n in range(f['ta'].shape[0]):
            for j, i in idx:
                ps = f['ps'][n, j, i] * 100.0
                mass = ps * dsig / c['ga']
                heat = f['dt_lw'][n, :, j, i] + f['dt_sw'][n, :, j, i]
                moist = f['dq_dyn'][n, :, j, i] * 1000.0 / mass             # m/s of water per layer to kg/kg/s
                shf, lhf = -f['hfss'][n, j, i], -f['hfls'][n, j, i]
                heat = heat.copy()
                moist = moist.copy()
                heat[-1] += shf / (c['acpd'] * mass[-1])
                moist[-1] += lhf / c['ALV'] / mass[-1]
                out['t'].append(f['ta'][n, :, j, i])
                out['q'].append(f['hus'][n, :, j, i])
                out['ps'].append(ps)
                out['shf'].append(shf)
                out['lhf'].append(lhf)
                out['heat'].append(heat)
                out['moist'].append(moist)
                out['omega'].append(om[n, :, j, i])
                out['wind'].append(math.hypot(f['uas'][n, j, i], f['vas'][n, j, i]))
                out['prc'].append(f['prc'][n, j, i] * 1000.0 * DAY)
                out['pr'].append(f['pr'][n, j, i] * 1000.0 * DAY)
                out['hour'].append(hour[n])
                out['w'].append(math.cos(math.radians(lat[j])))
                out['year'].append(year)
    res = {k: np.asarray(v) for k, v in out.items()}
    res.update(sigma=sigma, edges=edges, constants=c, cells=cells)
    return res


def run_gcm(site: str, g: dict, prm=PARAMETERS) -> dict:
    """The scheme on the GCM's site columns by quarter of the lunar day, beside the GCM's own rain, with PlaSim's
    evaporation of falling rain at the site's CM1 constant and at PlaSim's 0.01."""
    c = g['constants']
    sigma, edges = g['sigma'], g['edges']
    land = site == 'land'
    acc = {q[0]: dict(w=[], formed=[], ground_site=[], ground_plasim=[], prc=[], pr=[], type=[], top=[], mb=[])
           for q in QUARTERS}
    for k in range(len(g['ps'])):
        ps = g['ps'][k]
        x = dict(t=g['t'][k], q=g['q'][k], p_full=sigma * ps, p_half=edges * ps, z_full=None, z_half=None,
                 g=c['ga'], shf=g['shf'][k], lhf=g['lhf'][k], bl_dtdt=g['heat'][k], bl_dqdt=g['moist'][k],
                 omega=g['omega'][k], u_bl=g['wind'][k])
        r = one_output(x, land, dict(site=GAMMA[site], plasim=GAMMA['plasim']), prm, c)
        a = acc[quarter_name(g['hour'][k])]
        a['w'].append(g['w'][k])
        a['formed'].append(r['formed'])
        a['ground_site'].append(r['ground']['site'])
        a['ground_plasim'].append(r['ground']['plasim'])
        a['prc'].append(g['prc'][k])
        a['pr'].append(g['pr'][k])
        a['type'].append(r['type'])
        a['top'].append(r['z_top'] if r['z_top'] is not None else np.nan)
        a['mb'].append(r['m_base'])
    out = {}
    for name, a in acc.items():
        w = np.asarray(a['w'])
        mean = lambda v: float(np.sum(np.asarray(v, float) * w) / np.sum(w))
        types = np.asarray([str(t) for t in a['type']])
        tops = np.asarray(a['top'], float)
        deep = types == 'deep'
        out[name] = dict(columns=int(w.size), gcm_convective_rain_mm_day=mean(a['prc']), gcm_rain_mm_day=mean(a['pr']),
                         formed_mm_day=mean(a['formed']), ground_gamma_site_mm_day=mean(a['ground_site']),
                         ground_gamma_plasim_mm_day=mean(a['ground_plasim']),
                         type_share={t: float(np.sum(w[types == t]) / np.sum(w)) for t in ('deep', 'shallow', 'mid', 'None')},
                         deep_top_median_m=float(np.nanmedian(tops[deep])) if deep.any() else None,
                         base_mass_flux_mean=mean(a['mb']))
    return out


# --------------------------------------------------------------------------------------------- single column

def mix_boundary_layer(t, q, pf, ph, c, excess: float = 0.2):
    """A mixed boundary layer: the layers from the ground up to the last whose virtual potential temperature lies
    within `excess` K of the lowest layer's take the mean potential temperature and vapour of them all by mass (a
    stand-in for PlaSim's vertical diffusion and dry adjustment, enthalpy and water kept)."""
    kap = c['gascon'] / c['acpd']
    n = len(t)
    th = t * (1.0e5 / pf) ** kap
    thv = th * (1.0 + (1.0 / c['rdbrv'] - 1.0) * q)
    top = n - 1
    for k in range(n - 2, -1, -1):
        if thv[k] <= thv[-1] + excess:
            top = k
        else:
            break
    if top == n - 1:
        return t, q
    m = np.diff(ph)[top:]
    ex = (pf[top:] / 1.0e5) ** kap
    # keep the layers' enthalpy: the mean potential temperature weighted by mass and Exner
    th_m = np.sum(t[top:] * m) / np.sum(ex * m)
    q_m = np.sum(q[top:] * m) / np.sum(m)
    t, q = t.copy(), q.copy()
    t[top:] = th_m * ex
    q[top:] = q_m
    return t, q


def condense_large_scale(t, q, pf, mass, c):
    """Saturation adjustment of each supersaturated layer, the condensate raining out at once: the new temperature
    and vapour, and the rain (kg/m2/s over a unit second, i.e. kg/m2)."""
    rain = np.zeros(len(t))
    t, q = t.copy(), q.copy()
    for k in range(len(t)):
        for _ in range(3):
            qs = qsat(t[k], pf[k], c)
            if q[k] <= qs:
                break
            lv = c['ALS'] if t[k] < c['TMELT'] else c['ALV']
            dq = (q[k] - qs) / (1.0 + lv / c['acpd'] * dqsat_dt(t[k], pf[k], c))
            t[k] += lv / c['acpd'] * dq
            q[k] -= dq
            rain[k] += dq * mass[k]
    return t, q, rain


def scm(box: dict, tag: str, grid=None, days: float | None = None, prm=PARAMETERS, c=PLASIM, dt: float = STEP_S,
        start=None, cycles: int = 1):
    """A single column driven by a CM1 box's own radiation, day-night table, mean ascent (acting on the column's own
    air) and surface fluxes, interpolated in time between its outputs, with the scheme, PlaSim's evaporation of falling
    rain at the site's constant, a saturation adjustment and a mixed boundary layer; on CM1's levels or a PlaSim grid.
    It starts from the box's first output, or from `start` (temperature and vapour on the column's levels), and runs
    to its last (or `days`) `cycles` times over, the forcing repeating each lunar day. Returns the rain by quarter and
    the column's temperature and vapour of the last cycle against CM1's means over the span, and each cycle's."""
    rows = box['outputs']['rows']
    land = box['outputs']['land']
    gridded = [r if grid is None else to_grid(r, grid) for r in rows]
    for gr, r in zip(gridded, rows):
        if grid is not None:
            for key in ('table_dtdt', 'table_dqdt', 'w_mean'):
                gr[key] = to_grid(dict(r, t=r[key], q=r[key]), grid)['t']
        gr['day'], gr['hour'], gr['quarter'] = r['day'], r['hour'], r['quarter']
    first = gridded[0]
    t, q = (first['t'].copy(), first['q'].copy()) if start is None else (np.array(start['t'], float),
                                                                         np.array(start['q'], float))
    t_start = t.copy()
    pf, ph, g = first['p_full'], first['p_half'], first['g']
    mass = np.diff(ph) / g
    times = np.array([r['day'] for r in gridded]) * DAY
    end = times[-1] if days is None else times[0] + days * DAY
    history = []
    for _ in range(cycles):
        rain = {qn[0]: [0.0, 0.0] for qn in QUARTERS}
        acc_t, acc_q, steps = np.zeros_like(t), np.zeros_like(q), 0
        now = times[0]
        while now < end:
            t, q, fell, hour = scm_step(t, q, pf, ph, mass, g, gridded, times, now, land, tag, prm, c, dt)
            name = quarter_name(hour)
            rain[name][0] += fell * dt
            rain[name][1] += dt
            acc_t += t
            acc_q += q
            steps += 1
            now += dt
        history.append(dict(rain_mm_day={k: v[0] / v[1] * DAY if v[1] else None for k, v in rain.items()},
                            t_mean=acc_t / steps, q_mean=acc_q / steps))
    cm1_t = np.mean([r['t'] for r in gridded], axis=0)
    cm1_q = np.mean([r['q'] for r in gridded], axis=0)
    z_full, _ = heights(history[-1]['t_mean'], history[-1]['q_mean'], pf, ph, g, c)
    return dict(rain_mm_day=history[-1]['rain_mm_day'], z_m=z_full, t_mean=history[-1]['t_mean'], t_cm1=cm1_t,
                q_mean=history[-1]['q_mean'], q_cm1=cm1_q, t_end=t, t_start=t_start, steps=steps, cycles=history)


def scm_step(t, q, pf, ph, mass, g, gridded, times, now, land, tag, prm, c, dt):
    """One step of the single column: the forcing at time now, the scheme and the evaporation of its rain, a
    saturation adjustment and the mixed boundary layer. Returns the new temperature and vapour, the rain and snow
    reaching the ground (kg/m2/s) and the hour angle."""
    kap = c['gascon'] / c['acpd']
    j = int(np.clip(np.searchsorted(times, now) - 1, 0, len(times) - 2))
    f = (now - times[j]) / (times[j + 1] - times[j])
    lerp = lambda key: (1.0 - f) * gridded[j][key] + f * gridded[j + 1][key]
    hour = gridded[j]['hour'] + f * ((gridded[j + 1]['hour'] - gridded[j]['hour'] + 180.0) % 360.0 - 180.0)
    z_full, z_half = heights(t, q, pf, ph, g, c)
    w = lerp('w_mean')
    th = t * (1.0e5 / pf) ** kap
    heat = lerp('rad_dtdt') + lerp('table_dtdt') - w * np.gradient(th, z_full) * (pf / 1.0e5) ** kap
    moist = lerp('table_dqdt') - w * np.gradient(q, z_full)
    shf, lhf = lerp('shf'), lerp('lhf')
    heat[-1] += shf / (c['acpd'] * mass[-1])
    moist[-1] += lhf / c['ALV'] / mass[-1]
    omega = -(pf / (c['gascon'] * t)) * g * w
    r = column(t, q, pf, ph, g, dt, shf, lhf, heat, moist, omega, lerp('u_bl'), land, prm=prm, c=c, z_full=z_full,
               z_half=z_half, melt='plasim')
    ev = evaporate(r['column'], r['net_rain'], r['net_snow'], r['dtdt'], r['dqdt'], GAMMA[tag])
    t = t + dt * (heat + r['dtdt'] + ev['dtdt'])
    q = np.maximum(q + dt * (moist + r['dqdt'] + ev['dqdt']), 1.0e-9)
    t, q, ls_rain = condense_large_scale(t, q, pf, mass, c)
    t, q = mix_boundary_layer(t, q, pf, ph, c)
    return t, q, ev['rain'] + ev['snow'] + ls_rain.sum() / dt, hour


# ------------------------------------------------------------------------------------------------ the record

RESULTS = HERE.parent / 'results' / 'gcm'
OUT = RESULTS / 'massflux_offline.json'
SCHEMA = 'terluna.climate.gcm-massflux-offline/1'
SCM_HEIGHTS_M = (2000.0, 4000.0, 6000.0, 8000.0, 10000.0, 13500.0, 19000.0, 29000.0, 42000.0, 56000.0, 70000.0)
SCM_CYCLES = 3
VARIANTS = {'undiluted trigger': dict(sub_factor=0.0), 'no boundary-layer term': dict(t_star=1.0e12),
            "the first round's surface parcel": dict(source_height_earth=0.0, excess_scale='mixed_layer')}
# On fifteen levels: the first round's surface parcel (the lowest level's air with the excess at the sub-cloud layer's
# w*), the IFS's excess at the lowest level's height without the similarity, and the adopted parcel undiluted below
# its condensation level.
VARIANTS_15 = {"the first round's surface parcel": dict(source_height_earth=0.0, excess_scale='mixed_layer'),
               "the excess at the lowest level, no similarity": dict(source_height_earth=0.0),
               'undiluted below the condensation level': dict(sub_factor=0.0)}
LITERATURE = (
    ('IFS Documentation Cy36r1 (2010), Part IV, chapter 6', 'read', 'session-tools/agents/climate-explain/'
     'ifs_cy36r1_part4.txt on the research drive, sections 6.1-6.11'),
    ('IFS Documentation Cy40r1 (2014), Part IV, chapter 6', 'read', 'ecmwf.int elibrary 9204, sha256 f9fec6f0246be4f3; '
     'sections 6.3-6.4: one entrainment, the RH-dependent detrainment, the deep test parcel, PCAPE and its '
     'boundary-layer term, the shallow closure'),
    ('Bechtold and others (2014), J. Atmos. Sci. 71, 734-753', 'read as ECMWF Technical Memorandum 705 (2013), the '
     'version submitted to the journal', 'ecmwf.int elibrary 8018, sha256 2e057dd36bfae7bd; sections 2 and 5, appendix'),
    ('Bechtold and others (2008), QJRMS 134, 1337-1351', 'read as ECMWF Technical Memorandum 556 (2008)',
     'ecmwf.int elibrary 8011, sha256 431b13baaf05e12e; section 2.1'),
    ('Tiedtke (1989), Mon. Wea. Rev. 117, 1779-1800', 'not read', 'the AMS site refused the download; taken as the '
     'IFS documentation describes it'),
)


def tau_measure(box: dict, prm=PARAMETERS, c=PLASIM) -> dict:
    """CM1's adjustment time for the scheme's PCAPE, by quarter: on each output where the scheme convects deep, the
    rate at which CM1's own convective tendencies (microphysics, resolved transport, boundary layer) would remove the
    scheme's PCAPE with its updraught held, the PCAPE (less PCAPEbl) over that rate, and its ratio to the cloud's
    turnover time; medians over the quarter's deep outputs with the share of them on which CM1 removes PCAPE."""
    rows = box['outputs']['rows']
    land = box['outputs']['land']
    acc = {q[0]: [] for q in QUARTERS}
    for x in rows:
        r = column(x['t'], x['q'], x['p_full'], x['p_half'], x['g'], STEP_S, x['shf'], x['lhf'], x['bl_dtdt'],
                   x['bl_dqdt'], x['omega'], x['u_bl'], land, prm=prm, c=c, z_full=x['z_full'], z_half=x['z_half'])
        if r['type'] != 'deep':
            continue
        rate = pcape_rate(r['column'], r['up'], dict(dtdt=x['conv_dtdt'], dqdt=x['conv_dqdt']))
        acc[x['quarter']].append((r['pcape'], r['pcape_bl'], r['tau_c'], rate))
    out = {}
    for name, v in acc.items():
        if not v:
            out[name] = dict(deep_outputs=0)
            continue
        a = np.array(v)
        pos = a[:, 3] > 0.0
        tau = np.where(pos, a[:, 0] / np.where(pos, a[:, 3], 1.0), np.nan)
        out[name] = dict(deep_outputs=int(len(v)), pcape_median_pa=float(np.median(a[:, 0])),
                         turnover_median_h=float(np.median(a[:, 2])) / 3600.0,
                         cm1_removal_median_pa_s=float(np.median(a[:, 3])), cm1_removes_share=float(pos.mean()),
                         tau_cm1_median_h=float(np.nanmedian(tau)) / 3600.0 if pos.any() else None,
                         ratio_to_turnover_median=float(np.nanmedian(tau / a[:, 2])) if pos.any() else None)
    return out


def surface_layer_excess(box: dict, half_sigma, prm=PARAMETERS, c=PLASIM) -> dict:
    """By quarter, how much richer in moist static energy CM1's lowest level's air (50 m) is (kJ/kg, median and 10th
    and 90th percentiles) than the mean by mass of PlaSim's lowest layer over the same column, the air the grid's
    lowest level would hold; and what the surface test parcel adds to that layer's air: the surface-layer similarity
    (surface_layer_air, with the IFS's friction velocity and with CM1's own, and with CM1's in the Dyer-Hicks and
    neutral forms), the flux excess at the source height, and CM1's 50 m air less the parcel."""
    acc = {q[0]: dict(layer=[], similarity=[], similarity_cm1_ustar=[], dyer_hicks=[], neutral=[], excess=[], parcel=[])
           for q in QUARTERS}
    pct = lambda v: dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))
    for x in box['outputs']['rows']:
        xg = to_grid(x, half_sigma)
        g = x['g']
        h0 = c['acpd'] * x['t'][-1] + g * x['z_full'][-1] + c['ALV'] * x['q'][-1]
        col = Column(xg['t'], xg['q'], xg['p_full'], xg['p_half'], g, c)
        hl = c['acpd'] * col.t[-1] + g * col.zf[-1] + c['ALV'] * col.q[-1]
        dt_l, dq_l = surface_layer_air(col, x['shf'], x['lhf'], prm)
        ust = max(x.get('ust', prm['ustar']), 0.01)
        dt_u, dq_u = surface_layer_air(col, x['shf'], x['lhf'], dict(prm, ustar=ust))
        dt_d, dq_d = surface_layer_air(col, x['shf'], x['lhf'], dict(prm, ustar=ust, stability='dyer_hicks'))
        dt_0, dq_0 = surface_layer_air(col, x['shf'], x['lhf'], dict(prm, ustar=ust, stability='neutral'))
        dt_s, dq_s, _ = surface_excess(col, x['shf'], x['lhf'], prm)
        excess = c['acpd'] * max(prm['dt_deep'], dt_s) + c['ALV'] * max(prm['dq_deep'], dq_s)
        a = acc[x['quarter']]
        a['layer'].append((h0 - hl) / 1e3)
        a['similarity'].append((c['acpd'] * dt_l + c['ALV'] * dq_l) / 1e3)
        a['similarity_cm1_ustar'].append((c['acpd'] * dt_u + c['ALV'] * dq_u) / 1e3)
        a['dyer_hicks'].append((c['acpd'] * dt_d + c['ALV'] * dq_d) / 1e3)
        a['neutral'].append((c['acpd'] * dt_0 + c['ALV'] * dq_0) / 1e3)
        a['excess'].append(excess / 1e3)
        a['parcel'].append((h0 - hl - c['acpd'] * dt_l - c['ALV'] * dq_l - excess) / 1e3)
    return {k: dict(h50_minus_layer_kj_kg=pct(v['layer']), similarity_kj_kg=pct(v['similarity']),
                    similarity_cm1_ustar_kj_kg=pct(v['similarity_cm1_ustar']),
                    dyer_hicks_cm1_ustar_kj_kg=pct(v['dyer_hicks']), neutral_cm1_ustar_kj_kg=pct(v['neutral']),
                    excess_kj_kg=pct(v['excess']),
                    h50_minus_parcel_kj_kg=pct(v['parcel']),
                    correlation_layer_similarity=float(np.corrcoef(v['layer'], v['similarity'])[0, 1])
                    if np.std(v['similarity']) > 0.0 else None)
            for k, v in acc.items()}


def trigger_record(box: dict, half_sigma, prm=PARAMETERS, c=PLASIM) -> dict:
    """By quarter, the surface test parcel of the deep type on a box's columns laid on a grid: the share whose kinetic
    energy lasts to the first half level above its condensation level, the inhibition it meets below its level of free
    convection (J/kg, negative buoyancy integrated by height, median), the energy it starts with (median), and CM1's own
    CIN (median) for its 50 m air."""
    acc = {q[0]: dict(ok=[], cin=[], cm1=[], w0=[]) for q in QUARTERS}
    for x in box['outputs']['rows']:
        xg = to_grid(x, half_sigma)
        col = Column(xg['t'], xg['q'], xg['p_full'], xg['p_half'], x['g'], c)
        dt_s, dq_s, ws = surface_excess(col, x['shf'], x['lhf'], prm)
        dt_l, dq_l = surface_layer_air(col, x['shf'], x['lhf'], prm)
        w0 = surface_start(col, x['shf'], x['lhf'], ws, dt_l, dq_l, prm)
        acc[x['quarter']]['w0'].append(w0)
        tr = trial(col, col.n - 1, col.t[-1] + dt_l + max(prm['dt_deep'], dt_s),
                   col.q[-1] + dq_l + max(prm['dq_deep'], dq_s), w0, 'deep', rain=True, prm=prm)
        pts = tr['points']
        cin = 0.0
        for j in range(1, len(pts)):
            if pts[j]['b'] > 0.0 and pts[j]['l'] > 0.0:
                break
            cin -= 0.5 * (min(pts[j - 1]['b'], 0.0) + min(pts[j]['b'], 0.0)) * (pts[j]['z'] - pts[j - 1]['z'])
        a = acc[x['quarter']]
        a['ok'].append(bool(tr['ok']))
        a['cin'].append(cin)
        a['cm1'].append(x['cin'])
    return {k: dict(reaches_cloud_share=float(np.mean(v['ok'])), inhibition_median_j_kg=float(np.median(v['cin'])),
                    starting_energy_j_kg=float(np.median(0.5 * np.square(v['w0']))),
                    buoyancy_share_in_kinetic_energy=1.0 / (prm['f_turb'] * (1.0 + prm['gamma_vm'])),
                    cm1_cin_median_j_kg=float(np.median(v['cm1'])))
            for k, v in acc.items()}


COLD_POOL_K = 1.0                 # K: a column whose density potential temperature at 50 m lies this far below the
                                  # level's mean is in a cold pool (and, as a check, half of it)
REPORT_Z_M = (50.0, 250.0, 550.0, 1050.0, 2050.0, 3050.0, 4050.0, 5050.0, 6050.0)


def cold_pool_record(case: str, span) -> dict:
    """CM1's boundary layer by quarter from its 3-D output (box_entrainment's thermodynamics): the share of columns
    in cold pools (COLD_POOL_K), the 50 m air's moist static energy in the warm sector and in the cold pools against the
    domain's (kJ/kg), the domain's, the warm sector's and the cold pools' moist static energy and density potential
    temperature by height against the domain's 50 m values, and the air rising without cloud (w above 0.5 m/s) by
    height."""
    from climate.crm import box_entrainment as be
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    info = be.case_info(case)
    cfg = info['record']['configuration']
    zh = info['zh']
    idx = [int(np.argmin(abs(zh - z))) for z in REPORT_Z_M]
    keys = ('dom_h', 'dom_th', 'warm_h', 'warm_th', 'warm_n', 'cold_h', 'cold_th', 'cold_n', 'rise_h', 'rise_n')
    acc = {q[0]: dict({k: np.zeros(zh.size) for k in keys}, snaps=0, cold=[], cold_half=[], h50=[]) for q in QUARTERS}
    for n in sorted(int(f.name[8:14]) for f in info['case'].glob('cm1out_t*_s.dat')):
        t = (n - 1) * info['tap']
        if not span[0] * DAY + 1.0 < t <= span[1] * DAY + 1.0:
            continue
        x = be.thermo(ra.read_snapshot(info['case'], n), zh, info['g'])
        a = acc[quarter_name(local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], t, solar_day_s()))]
        h, th = x['h'], x['thrho']
        anom = th[0] - th[0].mean()
        cold = anom < -COLD_POOL_K
        warm = ~cold
        a['snaps'] += 1
        a['dom_h'] += h.mean(axis=1)
        a['dom_th'] += th.mean(axis=1)
        for name, m in (('warm', warm), ('cold', cold)):
            a[f'{name}_h'] += h[:, m].sum(axis=1)
            a[f'{name}_th'] += th[:, m].sum(axis=1)
            a[f'{name}_n'] += m.sum()
        rise = (x['w'] > 0.5) & (x['qcl'] < ra.CLOUD_KG_KG)
        a['rise_h'] += (h * rise).sum(axis=1)
        a['rise_n'] += rise.sum(axis=1)
        a['cold'].append(cold.mean())
        a['cold_half'].append((anom < -0.5 * COLD_POOL_K).mean())
        a['h50'].append((h[0, warm].mean() - h[0].mean(), h[0, cold].mean() - h[0].mean() if cold.any() else np.nan))
    out = dict(heights_m=[float(zh[i]) for i in idx], cold_pool_k=COLD_POOL_K)
    for name, a in acc.items():
        if not a['snaps']:
            continue
        dom_h, dom_th = a['dom_h'] / a['snaps'], a['dom_th'] / a['snaps']
        per = lambda key: a[f'{key}_h'] / np.maximum(a[f'{key}_n'], 1)
        ths = lambda key: a[f'{key}_th'] / np.maximum(a[f'{key}_n'], 1)
        h50 = np.array(a['h50'], float)
        out[name] = dict(
            snapshots=a['snaps'], cold_pool_share=float(np.mean(a['cold'])),
            cold_pool_share_half=float(np.mean(a['cold_half'])),
            h50_warm_minus_domain_kj_kg=float(np.nanmean(h50[:, 0])) / 1e3,
            h50_cold_minus_domain_kj_kg=float(np.nanmean(h50[:, 1])) / 1e3 if np.isfinite(h50[:, 1]).any() else None,
            domain_h_minus_h50_kj_kg=[float(dom_h[i] - dom_h[0]) / 1e3 for i in idx],
            warm_h_minus_h50_kj_kg=[float(per('warm')[i] - dom_h[0]) / 1e3 for i in idx],
            cold_h_minus_h50_kj_kg=[float(per('cold')[i] - dom_h[0]) / 1e3 if a['cold_n'][i] else None for i in idx],
            rising_h_minus_h50_kj_kg=[float(per('rise')[i] - dom_h[0]) / 1e3 if a['rise_n'][i] > 20 else None
                                      for i in idx],
            domain_theta_rho_minus_50m_k=[float(dom_th[i] - dom_th[0]) for i in idx],
            warm_theta_rho_minus_domain_k=[float(ths('warm')[i] - dom_th[i]) for i in idx])
    return out


def cm1_targets(box: dict, samples: dict) -> dict:
    out = {}
    for name, s in samples.items():
        rows = [x for x in box['outputs']['rows'] if x['quarter'] == name]
        out[name] = dict(rain_ground_mm_day=box['micro'][name]['ground_mm_day'],
                         falling_formed_mm_day=box['micro'][name]['falling_formed_mm_day'],
                         condensed_mm_day=box['micro'][name]['condensed_mm_day'],
                         cloudy_updraught_mass_flux=np.interp(PROFILE_Z_M, s['zh'], s['updraft']).tolist(),
                         core_mass_flux=np.interp(PROFILE_Z_M, s['zh'], s['core']).tolist(),
                         cloud_top_median_m=s['top_median_m'], cloud_top_p90_m=s['top_p90_m'],
                         cloudy_share=s['cloudy_share'], top_above_16km_share=s['top_above_16km_share'],
                         cape_j_kg=float(np.mean([x['cape'] for x in rows])),
                         cin_j_kg=float(np.mean([x['cin'] for x in rows])))
    return out


def scm_record(boxes: dict, gcm: dict, grids: dict, prm=PARAMETERS, c=PLASIM) -> dict:
    """The single column on CM1's levels and the PlaSim grids, from each box's own first output for one lunar day
    and from the GCM's site mean for SCM_CYCLES lunar days, its temperature against CM1's mean at SCM_HEIGHTS_M."""
    from climate.gcm import levels as lv
    out = {}
    for tag, box in boxes.items():
        g = gcm[tag]
        w = g['w'] / g['w'].sum()
        t10, q10 = np.tensordot(w, g['t'], 1), np.tensordot(w, g['q'], 1)
        for gname in ('cm1', '10', '14', '15'):
            grid = None if gname == 'cm1' else grids[gname]
            if gname == '10':
                st = dict(t=t10, q=q10)
            elif gname in ('14', '15'):
                st = dict(t=lv.remap(t10, grids['10'], grids[gname]), q=lv.remap(q10, grids['10'], grids[gname]))
            else:
                x = box['outputs']['rows'][0]
                lnp10 = np.log(lv.full_levels(grids['10']) * x['p_half'][-1])
                st = dict(t=np.interp(np.log(x['p_full']), lnp10, t10), q=np.interp(np.log(x['p_full']), lnp10, q10))
            for label, start, ncyc in (('from CM1', None, 1), ('from the GCM', st, SCM_CYCLES)):
                s = scm(box, tag, grid, prm=prm, c=c, start=start, cycles=ncyc)
                idx = [int(np.argmin(abs(s['z_m'] - z))) for z in SCM_HEIGHTS_M]
                out[f'{tag}, {gname} levels, {label}'] = dict(
                    lunar_days=ncyc, rain_mm_day=s['rain_mm_day'],
                    height_m=[float(s['z_m'][k]) for k in idx],
                    start_minus_cm1_k=[float(s['t_start'][k] - s['t_cm1'][k]) for k in idx],
                    mean_minus_cm1_k_by_lunar_day=[[float(h['t_mean'][k] - s['t_cm1'][k]) for k in idx]
                                                   for h in s['cycles']],
                    vapour_mean_minus_cm1_g_kg=[float((s['q_mean'][k] - s['q_cm1'][k]) * 1e3) for k in idx])
    return out


def cost_record(boxes: dict, gcm: dict, grids: dict, c=PLASIM) -> dict:
    """The scheme's work per column against the simplified Betts-Miller scheme's, as exponentials evaluated (each
    saturation, mixing and decay takes one), on the GCM's site columns and on the boxes' columns laid on both grids;
    the Betts-Miller scheme's 16 steps a layer cost a quarter more model time a year than one step (climate/gcm/
    README.md), which scales the estimate."""
    import math as _math

    class Counter:
        n = 0

        def __getattr__(self, name):
            return getattr(_math, name)

        def exp(self, x):
            Counter.n += 1
            return _math.exp(x)
    global math
    saved = math
    counter = Counter()
    out = {}
    try:
        math = counter
        sbm.math = counter

        def count(fn, cols):
            Counter.n = 0
            for x in cols:
                fn(x)
            return Counter.n / max(len(cols), 1)
        for tag in ('sea', 'land'):
            g = gcm[tag]
            s, e = g['sigma'], g['edges']
            cols = [k for k in range(0, len(g['ps']), 11)]
            out[f'GCM columns, {tag}'] = dict(
                massflux=count(lambda k: column(g['t'][k], g['q'][k], s * g['ps'][k], e * g['ps'][k], g['constants']['ga'],
                                                STEP_S, g['shf'][k], g['lhf'][k], g['heat'][k], g['moist'][k],
                                                g['omega'][k], g['wind'][k], tag == 'land', c=g['constants']), cols),
                sbm_one_step=count(lambda k: sbm.column(g['t'][k], g['q'][k], s, np.diff(e), g['ps'][k], rh=0.7,
                                                        akap=0.286, c=g['constants']), cols),
                sbm_16_steps=count(lambda k: sbm.column(g['t'][k], g['q'][k], s, np.diff(e), g['ps'][k], rh=0.7,
                                                        eps=2.4, steps=16, akap=0.286, c=g['constants']), cols))
            for gname in ('10', '14', '15'):
                xs = [to_grid(x, grids[gname]) for x in boxes[tag]['outputs']['rows'][::4]]
                out[f'CM1 columns on {gname} levels, {tag}'] = dict(massflux=count(
                    lambda x: column(x['t'], x['q'], x['p_full'], x['p_half'], x['g'], STEP_S, x['shf'], x['lhf'],
                                     x['bl_dtdt'], x['bl_dqdt'], x['omega'], x['u_bl'], tag == 'land', c=c), xs))
    finally:
        math = saved
        sbm.math = saved
    return out


def analyse() -> dict:
    from climate.gcm import sbm_offline as so
    c = so.run_constants(GCM_RUN)
    grids = grid_sets()
    boxes, samples = {}, {}
    for tag, case, span in BOXES:
        boxes[tag] = dict(outputs=cm1_outputs(case, span), micro=cm1_microphysics(case, span))
        samples[tag] = cm1_samples(case, span)
    gcm = {'sea': gcm_site_columns(GCM_RUN, GCM_YEARS, 'box_near_side_sea'),
           'land': gcm_site_columns(GCM_RUN, GCM_YEARS, 'box_0e_design')}
    out = dict(schema=SCHEMA, quarters=[q[0] for q in QUARTERS], profile_heights_m=list(PROFILE_Z_M),
               bands=[f'{a / 1000:g}-{b / 1000:g} km' for a, b in BANDS_M], parameters=dict(PARAMETERS),
               constants={k: c[k] for k in ('gascon', 'ga', 'acpd', 'rdbrv', 'ALV', 'ALS')}, step_s=STEP_S,
               gamma=GAMMA, grids=dict(source=grids['source'], half_10=grids['10'].tolist(),
                                       half_14=grids['14'].tolist(), half_15=grids['15'].tolist()),
               boxes={tag: dict(case=case, span_days=list(span), inputs=boxes[tag]['outputs']['inputs'])
                      for tag, case, span in BOXES},
               gcm=dict(run=GCM_RUN, years=list(GCM_YEARS)),
               literature=[dict(source=a, state=b, where=w) for a, b, w in LITERATURE])
    out['cm1'] = {tag: cm1_targets(boxes[tag], samples[tag]) for tag in boxes}
    out['tau'] = {tag: tau_measure(boxes[tag], c=c) for tag in boxes}
    out['on_cm1_levels'] = {tag: run_box(tag, boxes[tag], None, c=c) for tag in boxes}
    out['variants_on_cm1_levels'] = {
        label: {tag: {qn: dict(type_share=q['type_share'], precipitation_mm_day=q['precipitation_mm_day'],
                               cloud=q['cloud'])
                      for qn, q in run_box(tag, boxes[tag], None, prm=dict(PARAMETERS, **change), c=c).items()}
                for tag in boxes}
        for label, change in VARIANTS.items()}
    out['on_plasim_levels'] = {gname: {tag: run_box(tag, boxes[tag], grids[gname], c=c) for tag in boxes}
                               for gname in ('10', '14', '15')}
    out['variants_on_15_levels'] = {
        label: {tag: {qn: dict(type_share=q['type_share'], precipitation_mm_day=q['precipitation_mm_day'],
                               cloud=q['cloud'])
                      for qn, q in run_box(tag, boxes[tag], grids['15'], prm=dict(PARAMETERS, **change), c=c).items()}
                for tag in boxes}
        for label, change in VARIANTS_15.items()}
    out['surface_layer_excess'] = {gname: {tag: surface_layer_excess(boxes[tag], grids[gname], c=c) for tag in boxes}
                                   for gname in ('10', '14', '15')}
    out['trigger_on_15_levels'] = {tag: trigger_record(boxes[tag], grids['15'], c=c) for tag in boxes}
    out['cold_pools'] = {tag: cold_pool_record(case, span) for tag, case, span in BOXES}
    out['on_gcm_columns'] = {tag: run_gcm(tag, gcm[tag]) for tag in ('sea', 'land')}
    out['single_column'] = scm_record(boxes, gcm, grids, c=c)
    out['cost'] = cost_record(boxes, gcm, grids, c=c)
    out['evidence'] = EVIDENCE
    out['reading_rule'] = READING_RULE
    return out


EVIDENCE = ('Offline single-column calculations, no new CM1 or GCM runs. CM1: one realization each of the sea box\'s '
            'third lunar day and the land box\'s second, their domain means every 3 hours (59 a quarter) and their 3-D '
            'output (cloud samples, cold pools); the scheme meets columns CM1\'s own convection shaped, so its response '
            'at once and CM1\'s realised convection are compared. GCM: 3-day means of five years of the rh 0.7 branch '
            'in its site cells, built by the simplified Betts-Miller scheme; the boundary layer\'s forcing there lacks '
            'the dynamics\' heating and takes the radiation and the dynamics\' moistening as single-step snapshots. The '
            'single column is driven by CM1\'s radiation and surface fluxes, not its own, with a mixed boundary layer '
            'in place of PlaSim\'s vertical diffusion. The adjustment time\'s factor and the entrainment are measured '
            'on these boxes, so their sea quarters are not independent tests of them. The fifteen levels are '
            'levels.py\'s l15. Costs are counts of exponentials; the Fortran routine\'s timing is in massflux.md.')
READING_RULE = ('Quarters of the lunar day at each site by hour angle from noon (night -180 to -90, morning -90 to 0, '
                'afternoon 0 to 90, evening 90 to 180). Precipitation in mm/day: formed, what the scheme forms in its '
                'updraught; handed_on, net of its downdraught\'s evaporation, in the phase PlaSim\'s mkrain gives '
                'it; ground, after mkrain\'s evaporation (its copy, massflux.evaporate, with the floor of nprc 3) at '
                'the site\'s CM1 constant (gamma_site, sea 0.027, land 0.046) and at PlaSim\'s 0.01. CM1\'s '
                'falling_formed is the water falling from its clouds (evaporated on the way down plus reaching the '
                'ground, cm1out_stats). Mass fluxes in kg/m2/s at profile_heights_m, the scheme\'s updraught and '
                'downdraught (negative) means over all outputs, CM1\'s sampled cloudy updraught and core over its '
                'domain. Heating (K/day) and moistening (g/kg/day) profiles and the band sums (W/m2, mm/day): the '
                'scheme\'s convection with mkrain\'s evaporation and melting, against CM1\'s microphysics, resolved '
                'transport and boundary-layer scheme. type_share: shares of outputs by type (None: no convection). '
                'Cloud tops: the scheme\'s where its updraught\'s kinetic energy runs out; CM1\'s the top of each '
                'column\'s lowest unbroken cloud with its base below 16 km. tau: CM1\'s removal of the scheme\'s PCAPE '
                '(Pa/s) by its own convective tendencies, the updraught held. surface_layer_excess (kJ/kg, median and '
                '10th and 90th percentiles by quarter, on each grid): CM1\'s 50 m air less the grid\'s lowest layer, '
                'what the surface parcel adds to that layer (the similarity, with the IFS\'s friction velocity and '
                'with CM1\'s, the Dyer-Hicks and neutral forms with CM1\'s, and the excess), and CM1\'s 50 m air less '
                'the parcel. trigger_on_15_levels: the surface '
                'test parcel\'s share reaching its cloud with energy to spare, the inhibition it meets (J/kg, '
                'median), the energy it starts with and the share of its buoyancy its kinetic energy takes, and '
                'CM1\'s own CIN. cold_pools: by quarter from the 3-D output, the share of columns whose density '
                'potential temperature at 50 m lies 1 K (and 0.5 K) below the mean, the warm sector\'s and the cold '
                'pools\' 50 m moist static energy against the domain\'s, and by height (heights_m) the domain\'s, '
                'the warm sector\'s, the cold pools\' and the rising clear air\'s (w above 0.5 m/s) moist static '
                'energy against the domain\'s 50 m air, and density potential temperature. single_column: '
                'temperature of each lunar day\'s mean minus CM1\'s mean over the span (K) at the listed heights; '
                'rain at the ground by quarter in the last lunar day.')


def compact(x, digits: int = 5):
    if isinstance(x, dict):
        return {str(k): compact(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [compact(v, digits) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(x) else float(f'{float(x):.{digits}g}')
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def main(argv=None) -> int:
    import argparse
    import hashlib
    import json
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args(argv)
    out = compact(analyse())
    out['producer'] = dict(domain='climate', files={f'gcm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'sbm.py', HERE / 'sbm_offline.py')
                                                    + ((HERE / 'levels.py',) if (HERE / 'levels.py').exists() else ())})
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
