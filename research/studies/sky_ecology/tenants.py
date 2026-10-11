"""Tenants of aerophytes: real estate and its light, what a body carries, domatia, guilds, nursery and contact.

Pure functions plus evaluate(). Masses per m² of projected aerophyte area unless labelled; light as
photosynthetic photons (400-700 nm). The light is estimated from the illumination domain's clear-sky
surface-light product: its downward direct and diffuse spectra over a Lambertian ground and its column
reflectance for light from below, from which the air below a body gets its reflectance (see
underside_factor). The bodies' masses and trim come from the aerophyte product and model.
"""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from shared.constants import AVOGADRO, JULIAN_DAY, MOON_SURFACE_GRAVITY, PLANCK, SPEED_OF_LIGHT
from research.studies.aerophytes import growth as aero_growth
from research.studies.aerophytes import run as aero
from research.studies.aerophytes import water as aero_water

PAR_NM = (400., 700.)
ZONES = dict(top_cap=(0., 60.), upper_flank=(60., 90.), lower_flank=(90., 120.), under_cap=(120., 180.))
ALBEDOS = (.05, .1, .2, .5)                  # the surface-light product's grid: dark to pale ground, and cloud
RAFT_GAP = 1-math.pi/(2*math.sqrt(3))        # open share of a one-layer hexagonal raft of touching spheres
DUMPABLE_WATER_KG_M2 = aero.DUMPABLE_WATER_KG_M2
FREE_WATER_KG_M2 = aero.Traits().free_water_kg_m2
BROMELIAD_TANK_L_PER_HA = 50_000.           # tank bromeliads hold up to this in a cloud forest (Sugden & Robins 1979)
EPIPHYTE_KG_DRY_PER_M2_BRANCH = (1.3, 6.0)  # epiphyte mass on central branches (Freiberg & Freiberg 2000)
FRINGE_M2_PER_M2 = (.5, 2.)                 # design guess: sorbent and fibre fringe per projected m² (the water
                                            # note's one square metre of sorbent sits inside it)
FOUNDER_KG_DRY_M2 = (.01, .1)               # design guess: canopy life a juvenile carries out of its nursery
CONTACT_SPEED_M_S = (.1, 1.)                # design guess: relative drift, height changes and sailing between bodies
CONTACT_MARGIN_M = 50.                      # design guess: a close pass, edge to edge
DOMATIUM_LIVING_KG_M2 = .1                  # design guess: the light living layer that keeps a domatium's skin
DOMATIUM_GPP = .2                           # design guess: the layer's own carbon fixation, kg C/m² a year
DOMATIUM_FOOD = (0., .1)                    # design guess: food bodies for guards, kg C/m² a year
DOMATIUM_SUPPORT_SHARE = (.1, .3)           # design guess: share of a photosynthetic module's surplus it sends
TENANT_GAS_SHARE = (.1, .3)                 # design guess: share of a century's surplus a body spends on gas for tenants
BUDGET_YEARS = 100.
ADDED_MASS_SPHERE = .5                      # potential flow round a sphere
RAFT_DIAMETER_M = 2000.                     # a 2-km sky reef
DRAG_SPHERE = (.2, .5)                      # design range: a sphere at Reynolds numbers of 1e6-1e7
DRAG_PLATE = 1.2                            # a flat plate broadside: a raft moving up or down
VALVE_AREA_SHARE = (1e-5, 1e-4)             # design guess: open valve area per projected m² venting ballonet air
DRAIN_AREA_SHARE = (1e-5, 1e-4)             # design guess: open drain area per projected m² dropping water
DRAIN_HEAD_M = (.1, 1.)                     # design guess: the head of water over the drains
DISCHARGE = .6                              # a sharp-edged orifice's discharge coefficient
WATER_KG_M3 = 1000.
GAMMA = 1.4                                 # heat-capacity ratio of hydrogen and of air near 290 K
CV_PER_R = 2.5                              # diatomic hydrogen and air near 290 K: c_v = (5/2) R, so rho c_v = (5/2) p / T
HEAT_TRANSFER_W_M2_K = (2., 10.)            # design guess: the envelope's exchange with the air, convection and radiation
SURFACE_PER_PROJECTED = dict(round=4., module=2*math.pi/math.sqrt(3))   # a sphere; a module per hexagonal cell


def _check(**values):
    for k, v in values.items():
        if not (isinstance(v, (int, float)) and math.isfinite(v) and v >= 0):
            raise ValueError(f'{k} must be finite and nonnegative')


# ---- light -------------------------------------------------------------------------------------

def layer_reflectance(column_reflectance, layer_share):
    """Reflectance of the air below a body, per wavelength, for light from above.

    The conservative two-stream reflectance of a layer of depth tau is (3/4)tau/(1 + (3/4)tau). Inverting it
    on the column's computed reflectance gives an effective depth that folds the column's absorption in;
    the layer below a body takes `layer_share` (its share of the air column's mass) of that depth.
    """
    r = np.asarray(column_reflectance, float)
    if np.any((r < 0) | (r >= 1)) or not 0 <= layer_share <= 1:
        raise ValueError('reflectance in [0, 1) and share in [0, 1] required')
    tau_eff = r/(.75*(1-r))
    x = .75*layer_share*tau_eff
    return x/(1+x)


def underside_factor(layer_r, albedo):
    """Upward over downward light at the body's height: the layer over a Lambertian ground (adding)."""
    if not 0 <= albedo < 1:
        raise ValueError('albedo in [0, 1)')
    t = 1-layer_r
    return layer_r+t*t*albedo/(1-albedo*layer_r)


def two_stream_column(tau550, centres_nm, exponent=4.05):
    """The conservative two-stream reflectance of a pure Rayleigh column, for checking against the product."""
    tau = tau550*(550./np.asarray(centres_nm, float))**exponent
    x = .75*tau
    return x/(1+x)


def height_factor(layer_r, albedo):
    """Downward light at the body's height over the downward light at the ground, by two-stream adding.

    The air below is a conservative layer of reflectance R and transmittance 1 - R over a Lambertian ground:
    the ground receives D (1 - R) / (1 - A R) of the D arriving at the body's height.
    """
    if not 0 <= albedo < 1:
        raise ValueError('albedo in [0, 1)')
    r = np.asarray(layer_r, float)
    return (1-albedo*r)/(1-r)


def par_bands(surface_light, layer_share, albedo, case='moon_1.2atm/design', column='moon_1.2atm'):
    """PAR photons (umol/m²/s) by Sun elevation at the body's height: direct and diffuse down on a level surface,
    and up.

    The product's ground spectra over the same albedo are carried up to the body's height by height_factor,
    the direct beam and the diffuse light alike; the upward light is that downward light times
    underside_factor.
    """
    lam = np.asarray(surface_light['centres_nm'], float)
    par = (lam >= PAR_NM[0]) & (lam <= PAR_NM[1])
    photons = lam*1e-9/(PLANCK*SPEED_OF_LIGHT)/AVOGADRO*1e6
    rcol = np.asarray(surface_light['columns'][column]['reflectance_from_below'], float)
    r_layer = layer_reflectance(rcol, layer_share)
    lift = height_factor(r_layer, albedo)
    up = underside_factor(r_layer, albedo)
    rows = []
    for key, s in surface_light['cases'][case]['spectra'].items():
        d = np.asarray(s['direct_w_m2_nm'], float)
        ground = (d+np.asarray(s['diffuse_black_w_m2_nm'], float))/(1-albedo*rcol)
        down, direct = ground*lift, d*lift
        rows.append(dict(sun_deg=float(key), direct=float((direct*photons)[par].sum()),
                         diffuse=float(((down-direct)*photons)[par].sum()),
                         up=float((down*up*photons)[par].sum()),
                         ground_down=float((ground*photons)[par].sum())))
    return sorted(rows, key=lambda r: r['sun_deg'])


def _fibonacci_directions(n):
    i = np.arange(n)+.5
    z = 1-2*i/n
    phi = math.pi*(1+5**.5)*i
    s = np.sqrt(1-z*z)
    return np.stack([s*np.cos(phi), s*np.sin(phi), z], axis=1)


def _neighbour_centres(raft):
    if not raft:
        return np.zeros((0, 3))
    ring = [(2*math.cos(math.radians(60*k)), 2*math.sin(math.radians(60*k)), 0.) for k in range(6)]
    ring += [(2*math.sqrt(3)*math.cos(math.radians(30+60*k)), 2*math.sqrt(3)*math.sin(math.radians(30+60*k)), 0.)
             for k in range(6)]
    ring += [(4*math.cos(math.radians(60*k)), 4*math.sin(math.radians(60*k)), 0.) for k in range(6)]
    return np.array(ring)


def _visible(points, normals, directions, centres):
    """[element, direction] True where a ray leaves the surface into the open, missing every neighbour."""
    front = normals@directions.T > 1e-9
    if len(centres) == 0:
        return front
    vis = front.copy()
    origin = points+1e-6*normals
    for c in centres:
        oc = origin-c
        b = oc@directions.T
        cc = (oc*oc).sum(axis=1)[:, None]-1.
        disc = b*b-cc
        hit = (disc > 0) & (-b-np.sqrt(np.maximum(disc, 0)) > 0)
        vis &= ~hit
    return vis


def zone_light(bands, raft=False, n_beta=18, n_phi=24, n_dir=2000, n_time=91):
    """Daytime-mean PAR on the zones of a sphere, alone or inside a one-layer hexagonal raft, at the equator.

    The Sun crosses the zenith from east to west over the lunar day (equinox, equator), its elevation
    uniform in hour angle and sampled at the midpoints of n_time equal steps. Diffuse light from the sky and
    from below is isotropic; in a raft the six nearest and twelve next modules block the direct beam and part
    of each surface's view (ray tests).
    """
    beta = (np.arange(n_beta)+.5)*math.pi/n_beta
    phi = (np.arange(n_phi)+.5)*2*math.pi/n_phi
    b, p = np.meshgrid(beta, phi, indexing='ij')
    normals = np.stack([np.sin(b)*np.cos(p), np.sin(b)*np.sin(p), np.cos(b)], axis=-1).reshape(-1, 3)
    weights = np.sin(b).reshape(-1)
    centres = _neighbour_centres(raft)
    dirs = _fibonacci_directions(n_dir)
    vis = _visible(normals, normals, dirs, centres)
    cosw = np.maximum(normals@dirs.T, 0.)*vis*(4*math.pi/n_dir)/math.pi
    view_sky = (cosw*(dirs[:, 2] > 0)).sum(axis=1)
    view_ground = (cosw*(dirs[:, 2] < 0)).sum(axis=1)
    elev = np.array([r['sun_deg'] for r in bands])
    def interp(key, a):
        return np.interp(a, elev, np.array([r[key] for r in bands]), left=0.)
    hours = (np.arange(n_time)+.5)*180./n_time-90.
    a = 90.-np.abs(hours)
    sun = np.stack([np.where(hours < 0, np.cos(np.radians(a)), -np.cos(np.radians(a))),
                    np.zeros_like(a), np.sin(np.radians(a))], axis=1)
    sin_a = np.maximum(np.sin(np.radians(a)), 1e-9)
    beam = interp('direct', a)/sin_a
    lit = _visible(normals, normals, sun, centres)
    direct = (np.maximum(normals@sun.T, 0.)*lit*beam[None, :]).mean(axis=1)
    sky = view_sky*interp('diffuse', a).mean()
    ground = view_ground*interp('up', a).mean()
    total = direct+sky+ground
    top = float((interp('direct', a)+interp('diffuse', a)).mean())
    deg = np.degrees(np.repeat(beta, n_phi))
    zones = {}
    for name, (lo, hi) in ZONES.items():
        m = (deg >= lo) & (deg < hi)
        w = weights[m]
        zones[name] = dict(par_umol_m2_s=float((total[m]*w).sum()/w.sum()),
                           share_of_level_top=float((total[m]*w).sum()/w.sum()/top),
                           direct=float((direct[m]*w).sum()/w.sum()), sky=float((sky[m]*w).sum()/w.sum()),
                           from_below=float((ground[m]*w).sum()/w.sum()))
    ground = float(interp('ground_down', a).mean()) if 'ground_down' in bands[0] else None
    return dict(raft=raft, level_top_par_umol_m2_s=top,
                level_top_mol_m2_day_24h=top*.5*JULIAN_DAY/1e6, ground_par_umol_m2_s=ground,
                height_gain=top/ground if ground else None, zones=zones)


def light_table(surface_light, layer_share):
    """Zone light for a lone sphere and a raft module over the product's albedos, with the check data."""
    out = []
    for albedo in ALBEDOS:
        bands = par_bands(surface_light, layer_share, albedo)
        for raft in (False, True):
            z = zone_light(bands, raft=raft)
            z.update(albedo=albedo, upward_share_overhead=bands[-1]['up']/(bands[-1]['direct']+bands[-1]['diffuse']),
                     height_gain_overhead=(bands[-1]['direct']+bands[-1]['diffuse'])/bands[-1]['ground_down'])
            for v in z['zones'].values():
                v['mol_m2_day_24h'] = v['par_umol_m2_s']*.5*JULIAN_DAY/1e6
            out.append(z)
    return out


# ---- real estate -------------------------------------------------------------------------------

def sphere_zone_areas():
    """Area of each zone per projected m² of a lone sphere: 4 in all, one per zone."""
    return {k: (math.cos(math.radians(lo))-math.cos(math.radians(hi)))*2 for k, (lo, hi) in ZONES.items()}


def raft_zone_areas():
    """Area of each module zone per m² of a one-layer hexagonal raft (cell 2 sqrt(3) r²)."""
    cell = 2*math.sqrt(3)
    return {k: v*math.pi/cell for k, v in sphere_zone_areas().items()}


def raft_pocket_volume_m3_per_m2(module_radius_m):
    """Volume between touching modules above their equator plane per m² of raft: (sqrt3 - pi/3) r³ per
    triangle, two triangles a cell."""
    _check(r=module_radius_m)
    return 2*(math.sqrt(3)-math.pi/3)*module_radius_m**3/(2*math.sqrt(3)*module_radius_m**2)


def real_estate(light):
    """Surfaces and water by class per projected m², with their light."""
    sphere = sphere_zone_areas()
    raft = raft_zone_areas()
    def light_of(raft_flag, albedo):
        return next(z for z in light if z['raft'] == raft_flag and z['albedo'] == albedo)['zones']
    rows = []
    for cls, areas, flag in (('lone_sphere', sphere, False), ('raft_module', raft, True)):
        rows.append(dict(body=cls, zone_m2_per_m2=areas, fringe_m2_per_m2=list(FRINGE_M2_PER_M2),
                         pool_water_l_per_m2=[FREE_WATER_KG_M2-DUMPABLE_WATER_KG_M2, FREE_WATER_KG_M2],
                         light_share_of_level_top={str(a): {k: v['share_of_level_top'] for k, v in light_of(flag, a).items()}
                                                   for a in ALBEDOS},
                         gap_share=RAFT_GAP if flag else None,
                         pocket_m3_per_m2={str(r): raft_pocket_volume_m3_per_m2(r) for r in (30., 50., 75.)} if flag else None))
    return dict(rows=rows, bromeliad_tank_l_per_m2=BROMELIAD_TANK_L_PER_HA/1e4)


# ---- what a body carries -------------------------------------------------------------------------

def capacity_rows(size_rows, area_factor, b_kg_per_h2, carbon_per_h2, gate=.7):
    """Steady tenant load within the 30% reserve, its gas carbon, and static trim mass equivalents.

    The sudden_* fields are static trim mass equivalents; response_table gives how fast a body answers them: how
    far it sinks, its period, overshoot and settling, and how long venting or dropping water takes.
    The sudden loads are upper bounds. The hull's spare superpressure is the day's water swing's allowance
    (aerophytes/water.md), so the whole of it is free only at the low point of the swing; the 4 kg/m² of
    water a body can drop is the reserve its count of gas cells answers damage with (aerophytes README).
    """
    out = []
    for r in size_rows:
        if not r['all_gates']:
            continue
        d = r['diameter_m']
        area = math.pi*(d/2)**2*area_factor
        gross = r['supported_kg_m2']/r['gas_fraction']
        steady = (gate-r['gas_fraction'])*gross
        h2 = steady/b_kg_per_h2
        surplus = r['surplus_c_kg_m2_year']
        # Tenants a century's budget lifts: the low end spends the small share on fermented gas, the high end
        # the large share on photolytic gas; never more than the lift allows.
        budget = [min(steady, TENANT_GAS_SHARE[0]*surplus*BUDGET_YEARS/max(carbon_per_h2)*b_kg_per_h2),
                  min(steady, TENANT_GAS_SHARE[1]*surplus*BUDGET_YEARS/min(carbon_per_h2)*b_kg_per_h2)]
        out.append(dict(diameter_m=d, area_m2=area, supported_kg_m2=r['supported_kg_m2'], gas_fraction=r['gas_fraction'],
                        gross_lift_kg_m2=gross, steady_tenants_kg_m2=steady, steady_tenants_t=steady*area/1e3,
                        hydrogen_kg_m2=h2, gas_carbon_kg_m2=[h2*c for c in carbon_per_h2],
                        years_of_surplus=[h2*c/surplus for c in carbon_per_h2] if surplus > 0 else None,
                        century_budget_tenants_kg_m2=budget, century_budget_tenants_t=[x*area/1e3 for x in budget],
                        sudden_trim_kg_m2=r['water_swing_kg_m2_at_spare'], sudden_trim_t=r['water_swing_kg_m2_at_spare']*area/1e3,
                        sudden_with_water_t=(r['water_swing_kg_m2_at_spare']+DUMPABLE_WATER_KG_M2)*area/1e3,
                        spare_superpressure_pa=r['spare_superpressure_pa'], surplus_c_kg_m2_year=surplus))
    return out


def landing(capacity, flyers):
    """Mass-equivalent flyer counts: the flyers each static allowance holds (response_table gives the timing)."""
    rows = []
    for c in capacity:
        for f in flyers:
            if f['edge'] != 'low':
                continue
            rows.append(dict(diameter_m=c['diameter_m'], flyer=f['flyer'], mass_kg=f['mass_kg'],
                             at_once_by_trim=c['sudden_trim_t']*1e3/f['mass_kg'],
                             at_once_dropping_water=c['sudden_with_water_t']*1e3/f['mass_kg'],
                             roosting_steadily=c['steady_tenants_t']*1e3/f['mass_kg']))
    return rows


def scale_height_of_density(air_profile, height_km=10.):
    """The air's density scale height at a height, from the sky-ship product's levels on either side (m)."""
    levels = sorted(air_profile, key=lambda a: a['height_km'])
    below = max((a for a in levels if a['height_km'] < height_km), key=lambda a: a['height_km'])
    above = min((a for a in levels if a['height_km'] > height_km), key=lambda a: a['height_km'])
    return (above['height_km']-below['height_km'])*1e3/math.log(below['density_kg_m3']/above['density_kg_m3'])


def _damping(rho, drag, amplitude_m, effective_mass_kg_m2):
    """Equivalent damping ratio of quadratic drag at an amplitude, per projected m²: (2/3 pi) rho C_d a / m_eff."""
    return 2/(3*math.pi)*rho*drag*amplitude_m/effective_mass_kg_m2


def lapse_rate(air_profile, height_km=10.):
    """The air's temperature lapse at a height, from the sky-ship product's levels on either side (K per m)."""
    levels = sorted(air_profile, key=lambda a: a['height_km'])
    below = max((a for a in levels if a['height_km'] < height_km), key=lambda a: a['height_km'])
    above = min((a for a in levels if a['height_km'] > height_km), key=lambda a: a['height_km'])
    return (below['temperature_c']-above['temperature_c'])/((above['height_km']-below['height_km'])*1e3)


def step_peak_sink(mass, added_ratio, rho, drag, h_rho_m, share, period_s, steps=600):
    """Deepest point below the start after a load arrives at once (m): fourth-order Runge-Kutta on the nonlinear
    float, buoyancy from the air's exponential density and quadratic drag, to the first turning point."""
    g = MOON_SURFACE_GRAVITY
    m_eff = mass*(1+share+added_ratio)

    def acc(z, v):
        e = math.exp(-z/h_rho_m)
        return (g*mass*(e-1-share)-.5*rho*e*drag*v*abs(v))/m_eff
    dt = period_s/steps
    z = v = 0.
    for _ in range(4*steps):
        k1z, k1v = v, acc(z, v)
        k2z, k2v = v+.5*dt*k1v, acc(z+.5*dt*k1z, v+.5*dt*k1v)
        k3z, k3v = v+.5*dt*k2v, acc(z+.5*dt*k2z, v+.5*dt*k2v)
        k4z, k4v = v+dt*k3v, acc(z+dt*k3z, v+dt*k3v)
        z_new, v_new = z+dt/6*(k1z+2*k2z+2*k3z+k4z), v+dt/6*(k1v+2*k2v+2*k3v+k4v)
        if v < 0 <= v_new:
            return -min(z, z_new)
        z, v = z_new, v_new
    return -z


def vertical_response(column_m, trim_kg_m2, spare_pa, air, h_rho_m, added_mass_ratio=ADDED_MASS_SPHERE,
                      drag=DRAG_SPHERE, lapse_k_m=None, surface_per_projected=4.):
    """How a fixed-volume body answers a load equal to its trim mass equivalent, per projected m².

    column_m is the hull volume per projected m² (the trim rule's V/A), so the body's mass at neutral buoyancy is
    the air it displaces, rho x column. Left alone it sinks to where the air is denser by dM/M,
    dz = H_rho ln(1 + dM/M), and the ambient pressure rises by rho g dz, which its spare superpressure must
    absorb while the gas keeps its temperature; as the gas warms to the lower air over its thermal time the cost
    falls back toward the trim rule's (both lapses: the sky-ship temperatures', and the one the density profile
    implies, which restores it exactly). The restoring force, M g / H_rho per metre, against the inertia
    M (1 + C_a) gives a period 2 pi sqrt((1 + C_a) H_rho / g). A load that arrives at once overshoots: the deepest
    point comes from integrating the nonlinear float with quadratic drag, and the largest step the spare
    superpressure survives from bisecting on it. The swing decays as a / (1 + zeta omega t) with zeta the
    equivalent damping at the static sink, so nine-tenths of it are gone after 9 / (zeta omega). Trimming takes the
    load at the same height: ballonet air dM / rho vented through valves at the spare superpressure, which within
    minutes expands the gas adiabatically and so takes at most 1/gamma of the trim equivalent until the gas warms
    again; or water dropped through drains under its own head, which costs no superpressure. The gas warms at fixed
    volume, so its thermal time is rho c_v x column over the envelope's exchange per projected m².
    """
    if column_m <= 0 or h_rho_m <= 0:
        raise ValueError('positive column and scale height required')
    _check(trim=trim_kg_m2, spare=spare_pa, added=added_mass_ratio)
    rho, g = air['density_kg_m3'], MOON_SURFACE_GRAVITY
    p, t_k = air['pressure_atm']*101325., air['temperature_c']+273.15
    mass = rho*column_m
    share = trim_kg_m2/mass
    sink = h_rho_m*math.log1p(share)
    used = rho*g*sink
    warmed = [used-p*lapse_k_m*sink/t_k if lapse_k_m is not None else None, p*math.log1p(share)]
    omega = math.sqrt(g/(h_rho_m*(1+added_mass_ratio)))
    period = 2*math.pi/omega
    zetas = [_damping(rho, c, sink, mass*(1+added_mass_ratio)) for c in drag]
    peaks = [step_peak_sink(mass, added_mass_ratio, rho, c, h_rho_m, share, period) for c in drag]

    def step_share(c):                                          # largest step load the spare survives, of the trim
        if rho*g*step_peak_sink(mass, added_mass_ratio, rho, c, h_rho_m, share, period) <= spare_pa:
            return 1.
        lo, hi = 0., 1.
        for _ in range(24):
            mid = (lo+hi)/2
            ok = rho*g*step_peak_sink(mass, added_mass_ratio, rho, c, h_rho_m, mid*share, period) <= spare_pa
            lo, hi = (mid, hi) if ok else (lo, mid)
        return lo
    vent_speed = math.sqrt(2*spare_pa/rho)
    vent = [trim_kg_m2/rho/(DISCHARGE*a*vent_speed) for a in VALVE_AREA_SHARE[::-1]]
    drain = [DISCHARGE*a*WATER_KG_M3*math.sqrt(2*g*h) for a, h in ((DRAIN_AREA_SHARE[1], DRAIN_HEAD_M[1]),
                                                                  (DRAIN_AREA_SHARE[0], DRAIN_HEAD_M[0]))]
    heat = CV_PER_R*p/t_k                                       # J/m³/K: the hull heats its contents at fixed volume
    thermal = [heat*column_m/(h*surface_per_projected) for h in HEAT_TRANSFER_W_M2_K[::-1]]
    return dict(mass_kg_m2=mass, trim_kg_m2=trim_kg_m2, load_share_of_mass=share, static_sink_m=sink,
                superpressure_used_pa=used, used_over_spare=used/spare_pa if spare_pa else None,
                superpressure_used_warmed_pa=warmed,
                passive_gradual_share_of_trim=min(1., spare_pa/used) if used else 1.,
                passive_gradual_share_warmed=[min(1., spare_pa/w) if w else 1. for w in warmed[::-1]],
                gas_thermal_time_s=thermal,
                period_s=period, added_mass_ratio=added_mass_ratio, damping_ratio=zetas,
                step_peak_sink_m=peaks, step_overshoot=[pk/sink-1 for pk in peaks],
                step_superpressure_used_pa=[rho*g*pk for pk in peaks],
                passive_step_share_of_trim=[step_share(c) for c in drag],
                settle_to_tenth_s=[9/(z*omega) for z in zetas[::-1]],
                vent_s=vent, vent_at_once_share_of_trim=min(1., 1/GAMMA),
                vent_at_once_s=[t*min(1., 1/GAMMA) for t in vent],
                drop_trim_water_s=[trim_kg_m2/d for d in drain],
                drop_all_water_s=[DUMPABLE_WATER_KG_M2/d for d in drain])


def response_table(capacity, air, air_profile):
    """Vertical response to the trim mass equivalent: lone round bodies, lone reef modules, and reef modules in a
    2-km raft that shares the load through its ties. The raft's added mass is a solid disc's, (8/3) rho R^3, an
    upper bound: its 9.3% gaps and finite thickness lower it, and with it the raft's period."""
    h_rho = scale_height_of_density(air_profile)
    lapse = lapse_rate(air_profile)
    module_factor = 2*math.sqrt(3)/math.pi
    rows = []
    for kind, keep in (('round', (50., 100., 200., 500.)), ('round_strong', (1000.,)),
                       ('colony_module', (60., 100., 150.))):
        for c in capacity[kind]:
            if c['diameter_m'] not in keep:
                continue
            d = c['diameter_m']
            column = 2/3*d if kind != 'colony_module' else 4/3*(d/2)/module_factor
            r = vertical_response(column, c['sudden_trim_kg_m2'], c['spare_superpressure_pa'], air, h_rho,
                                  lapse_k_m=lapse, surface_per_projected=SURFACE_PER_PROJECTED[
                                      'module' if kind == 'colony_module' else 'round'])
            r.update(body=kind, diameter_m=d, load_t=c['sudden_trim_t'], area_m2=c['area_m2'], alone=True)
            rows.append(r)
            if kind == 'colony_module':
                raft_area = math.pi*(RAFT_DIAMETER_M/2)**2
                shared = c['sudden_trim_t']*1e3/raft_area             # the module's load spread over the raft
                added = 8/3*air['density_kg_m3']*(RAFT_DIAMETER_M/2)**3/(air['density_kg_m3']*column*raft_area)
                rr = vertical_response(column, shared, c['spare_superpressure_pa'], air, h_rho, added_mass_ratio=added,
                                       drag=(DRAG_PLATE, DRAG_PLATE))
                # At swings of centimetres drag hardly damps; the air's own motions set the swing, so only the
                # undamped bound of a step, twice the static sink, is kept, and the module's own valves and drains
                # are the lone module's.
                for k in ('damping_ratio', 'step_overshoot', 'settle_to_tenth_s', 'vent_s', 'drop_trim_water_s',
                          'drop_all_water_s', 'vent_at_once_s', 'superpressure_used_warmed_pa',
                          'passive_gradual_share_warmed', 'gas_thermal_time_s'):
                    rr[k] = None
                rr.update(body='colony_module_in_2km_raft', diameter_m=d, load_t=c['sudden_trim_t'],
                          area_m2=raft_area, alone=False, step_peak_sink_m=[2*rr['static_sink_m']]*2,
                          step_superpressure_used_pa=[2*rr['superpressure_used_pa']]*2,
                          passive_step_share_of_trim=[min(1., c['spare_superpressure_pa']/(2*rr['superpressure_used_pa']))]*2)
                rows.append(rr)
    return dict(density_scale_height_m=h_rho, lapse_k_per_m=lapse, rows=rows,
                design_guesses=dict(added_mass_sphere=ADDED_MASS_SPHERE, drag_sphere=list(DRAG_SPHERE),
                                    drag_plate=DRAG_PLATE, valve_area_share=list(VALVE_AREA_SHARE),
                                    drain_area_share=list(DRAIN_AREA_SHARE), drain_head_m=list(DRAIN_HEAD_M),
                                    discharge=DISCHARGE, raft_diameter_m=RAFT_DIAMETER_M, gamma=GAMMA,
                                    cv_per_r=CV_PER_R,
                                    heat_transfer_w_m2_k=list(HEAT_TRANSFER_W_M2_K)))


def domatium(diameter_m, air, structural, food_c=0., spare_pa=aero.SPARE_TRIM_PA):
    """A colony module grown for tenants: no green top, the round bodies' spare pressure, a full water store.

    Its light living layer keeps the skin and makes the module's gas by photolysis; the colony feeds the rest
    of its upkeep through the ties, as a siphonophore's zooids share one stem.
    """
    fam = aero.design_families(10, structural)
    t = replace(aero.colony_module(diameter_m/2, fam, spare_pa), name='domatium', living_dry_kg_m2=DOMATIUM_LIVING_KG_M2,
                community_dry_host_fraction=0., gpp_c_kg_m2_year=DOMATIUM_GPP, consumer_food_c_kg_m2_year=food_c)
    c, t, _ = aero.settle_case(diameter_m, air, t, cells=1)
    area = math.pi*(diameter_m/2)**2*t.cell_area_factor
    gross = c['mass']['full_gross_supported_mass_kg']/area
    f = c['mass']['required_lifting_mix_volume_fraction']
    column = 4/3*(diameter_m/2)/t.cell_area_factor
    swing = aero_water.allowed_swing(spare_pa, column, air['density_kg_m3'], air['pressure_atm']*101325.)
    return dict(diameter_m=diameter_m, food_to_tenants_c=food_c, supported_kg_m2=c['mass']['supported_kg_m2'],
                gas_fraction=f, steady_tenants_kg_m2=(.7-f)*gross, steady_tenants_t=(.7-f)*gross*area/1e3,
                sudden_trim_kg_m2=swing, sudden_trim_t=swing*area/1e3,
                sudden_with_water_t=(swing+DUMPABLE_WATER_KG_M2)*area/1e3,
                carbon_balance_c_kg_m2_year=c['carbon']['remaining_assimilate_c_kg_m2_year'],
                area_m2=area, mass_gate=c['gates']['mass_and_reserve'], solar_gas=c['gates']['solar_h2_allocation'])


def domatia_table(air, structural, module_rows):
    """Domatia at three sizes beside the photosynthetic modules that pay for them."""
    out = []
    green = {r['diameter_m']: r for r in module_rows if r['all_gates']}
    for d in (60., 100., 150.):
        for food in DOMATIUM_FOOD:
            dm = domatium(d, air, structural, food)
            g = green.get(d)
            if g:
                deficit = -min(0., dm['carbon_balance_c_kg_m2_year'])
                dm['green_modules_per_domatium'] = [deficit/(s*g['surplus_c_kg_m2_year']) for s in DOMATIUM_SUPPORT_SHARE[::-1]]
                dm['green_module_steady_tenants_kg_m2'] = (.7-g['gas_fraction'])*g['supported_kg_m2']/g['gas_fraction']
            out.append(dm)
    return out


# ---- guilds ---------------------------------------------------------------------------------------

def guilds(flyers, food_web_chain, light):
    """Each tenant guild with its Earth precedent and its costs and benefits to the host."""
    reef_under = next(z for z in light if z['raft'] and z['albedo'] == .1)['zones']['under_cap']['share_of_level_top']
    sphere_under = next(z for z in light if not z['raft'] and z['albedo'] == .1)['zones']['under_cap']['share_of_level_top']
    albatross = next(f for f in flyers if f['edge'] == 'low')
    tol = food_web_chain['tolerance']
    return [
        dict(guild='guards', where='domatia and flanks', precedent='ant-plants: Macaranga triloba puts about 5% of its '
             'daily above-ground production (9% of construction costs) into food bodies for its ants (Heil et al. 1997); '
             'ant debris and respiration give Dischidia major 29% of its nitrogen and 39% of its leaf carbon (Treseder et al. 1995)',
             light='none: they live in domatia', weight_kg_m2=[.001, .01], damage='the guards keep grazers and egg-layers '
             'off the green top (proposed; efficacy unmeasured)',
             mass_basis='unspecified in the initial screen; collective mass per footprint area, not individual size or a capacity limit',
             nutrients='debris and respiration in the domatia return nitrogen, carbon and some '
             'phosphorus', cost_c_kg_m2_year=[0., .1], cost_note='the food bodies are the reference organism\'s 0.1 kg C '
             'per m² a year shared among all tenants, not a separate allocation per guild'),
        dict(guild='roosting and nesting flyers', where='the flanks of round giants and landing domatia',
             precedent='seabird islands: colonies receive 99 Gg of phosphorus a year in droppings (Otero et al. 2018); '
             'islands without rats carry 760 times the seabirds and 251 times the nitrogen input, and their reefs 48% more '
             'fish (Graham et al. 2018); guano falls at 0.32-0.37 g/m² a day within 50-320 m of Spitsbergen colonies '
             '(Zwolicki et al. 2013)', light='droppings foul the green top unless the roosts are on the flanks',
             weight_kg_per_flyer=albatross['mass_kg'], damage='claws and nests wear the skin; the loads arrive at once',
             nutrients=f"a {albatross['mass_kg']/1e3:.1f}-t sea-feeding flyer excretes about "
                       f"{albatross['p_excreted_kg_year']:.0f} kg of phosphorus a year, a share of it at the roost"),
        dict(guild='pool food webs', where='the water store, held in the crevices between a reef\'s modules and in '
             'the catchments of round bodies', precedent='phytotelmata: tank bromeliads hold up to 50,000 litres a hectare '
             'in a cloud forest (Sugden & Robins 1979); a jumping spider living on a bromeliad supplies 18% of its host\'s '
             'nitrogen and lengthens its leaves 15% (Romero et al. 2006)', light='none', weight_kg_m2=[1., 5.],
             weight_note='the water is the host\'s own store, carried already', damage='mosquito-like larvae and '
             'pathogens breed in still water', nutrients='decomposers mineralise litter and frass in the pools and the '
             'host absorbs the phosphate, the retention the sky needs'),
        dict(guild='canopy soils', where='the flanks of giants', precedent='canopy soils: rain-forest trees grow roots '
             'into the soils on their own branches (Nadkarni 1981); epiphytes weigh 1.3-6.0 kg dry per m² of central '
             'branch (Freiberg & Freiberg 2000)', light='the flanks\' own light', weight_kg_m2=[2.6, 12.],
             weight_note='dry, on the two m² of flank per projected m², and more when wet: within the spare lift of round '
             'bodies from 150 m (47 kg/m²) up, against 25 kg/m² at 100 m', damage='weight and water on the skin',
             nutrients='hold litter and water; the host roots into them'),
        dict(guild='drifting-weed community', where='fringes and undersides', precedent='holopelagic Sargassum carries '
             'epiphytes, fungi, more than 100 invertebrate and over 100 fish species and four sea turtles (Coston-Clements '
             'et al. 1991); its nitrogen-fixing epiphytes may be the largest source of reactive nitrogen to the weed and its '
             'community (Johnson et al. 2023, after Phlips et al. 1986); 1,205 species are known to raft (Thiel & Gutow 2005)',
             light=f'the undersides\' {sphere_under:.0%} (lone body) and {reef_under:.0%} (raft) of the top\'s light',
             weight_kg_m2=[.05, .5], damage='epiphytes shade whatever they cover', nutrients='nitrogen fixers feed the host'),
        dict(guild='grazers', where='the green top', precedent='tropical forest canopies lose 12-19% of their leaf '
             'production to herbivores (Metcalfe et al. 2014); forests and shrublands lose 7.9% of their above-ground '
             'production (Cebrian 2004, in Metcalfe et al. 2014)', light='they remove green tissue',
             cost_share_of_surplus=[t['share_of_surplus'] for t in tol], damage='bites and the loss of unresorbed nutrients',
             nutrients='frass carries leaf phosphorus before resorption; caught in pools it returns, falling it is lost'),
        dict(guild='collectors', where='fringes and tethers', precedent='orb-weavers eat the pollen their webs catch '
             '(Eggs & Sanders 2013)', light='none', weight_kg_m2=[.001, .01], damage='none',
             nutrients='they turn airborne pollen, spores and insects into frass on the host'),
        dict(guild='bud dispersers', where='visitors', precedent='the canopy nursery\'s partnership: aerophytes and their '
             'megaforest hosts (aerophyte growth note)', light='none', damage='none',
             nutrients='they carry the host\'s buds to nurseries and are paid in food'),
    ]


# ---- nursery inoculation -------------------------------------------------------------------------

def inoculation(round_rows):
    """Lift a juvenile has for the canopy life it carries out of the nursery, against the founder load."""
    out = []
    for r in round_rows:
        if r['diameter_m'] not in (50., 60., 80.):
            continue
        gross = r['supported_kg_m2']/r['gas_fraction']
        spare = (.7-r['gas_fraction'])*gross
        out.append(dict(diameter_m=r['diameter_m'], spare_lift_kg_m2=spare, founder_kg_dry_m2=list(FOUNDER_KG_DRY_M2),
                        founder_wet_kg_m2=[x/.25 for x in FOUNDER_KG_DRY_M2], passes=r['all_gates']))
    return out


# ---- contact and disease -------------------------------------------------------------------------

def encounter_rate_per_day(diameter_m, cover, relative_speed_m_s, margin_m=CONTACT_MARGIN_M):
    """Close passes a day for one body among equal bodies at a cover: density x speed x swept width."""
    _check(d=diameter_m, u=relative_speed_m_s, margin=margin_m)
    if not 0 < cover <= 1:
        raise ValueError('cover in (0, 1]')
    density = cover/(math.pi*(diameter_m/2)**2)
    return density*relative_speed_m_s*2*(diameter_m+margin_m)*JULIAN_DAY


def contact_table(gathering, settling, start_km=10.):
    """Encounter rates, how fast passive routes gather, and how long spores stay aloft.

    The first screen's residence times are for a mixed layer of surface air; a particle shed at the body's
    height falls the whole start_km at its Stokes speed, which follows the air's viscosity and barely its
    density.
    """
    rows = []
    for d in (10., 100., 1000., 2000.):
        for cover in (.001, .01, .02):
            rows.append(dict(diameter_m=d, cover=cover,
                             passes_per_day=[encounter_rate_per_day(d, cover, u) for u in CONTACT_SPEED_M_S]))
    from_height = {name: start_km*1e3/p['fall_m_s']['moon']/JULIAN_DAY for name, p in settling['particles'].items()}
    return dict(encounters=rows, gathering=gathering, spores=settling,
                fall_days_from_body_height=dict(start_km=start_km, days=from_height),
                diadema=dict(area_km2=3.5e6, months=13, carrier='surface currents (Lessios 1988)'))


def evaluate(surface_light, air10, structural, aero_product, flyers, food_web_chain, layer_share, gathering,
             settling, air_profile):
    """The tenants' numbers for the product."""
    light = light_table(surface_light, layer_share)
    check = dict(product=[float(x) for x in np.asarray(surface_light['columns']['moon_1.2atm']['reflectance_from_below'])],
                 two_stream=[float(x) for x in two_stream_column(surface_light['columns']['moon_1.2atm']['rayleigh_optical_depth_550nm'],
                                                                 surface_light['centres_nm'])])
    lam = np.asarray(surface_light['centres_nm'])
    diff = np.array(check['product'])-np.array(check['two_stream'])
    blue_green = (lam >= 400) & (lam <= 560)
    par = (lam >= PAR_NM[0]) & (lam <= PAR_NM[1])
    worst = int(np.argmax(np.where(blue_green, np.abs(diff), -1.)))
    check = dict(max_abs_difference_400_560_nm=float(np.abs(diff)[blue_green].max()), at_worst_nm=float(lam[worst]),
                 share_of_400_560_nm_bins_within_0_01=float((np.abs(diff)[blue_green] <= .01).mean()),
                 median_abs_difference_400_560_nm=float(np.median(np.abs(diff)[blue_green])),
                 largest_excess_over_two_stream_par=float(diff[par].max()),
                 at_nm={str(n): dict(product=check['product'][int(np.argmin(abs(lam-n)))],
                                     two_stream=check['two_stream'][int(np.argmin(abs(lam-n)))])
                        for n in (400, 450, 500, 550, 600, 650, 700)})
    rows = aero_product['size_limits']['10_km']['rows']
    b = aero_product['trim_cycle']['net_supported_kg_per_hydrogen_kg']
    def costs(gpp):
        return [aero_growth.photo_carbon_per_hydrogen(gpp, aero.Traits().incident_cycle_mean_w_m2,
                                                      aero.Traits().assumed_full_solar_to_h2_lhv_efficiency),
                aero_growth.fermentation_carbon_per_hydrogen(2.1)]
    module_factor = 2*math.sqrt(3)/math.pi
    capacity = dict(round=capacity_rows(rows['round'], 1., b, costs(2.)),
                    round_strong=capacity_rows(rows['round_strong'], 1., b, costs(2.)),
                    round_strong_canopy_light=capacity_rows(rows['round_strong_gpp5'], 1., b, costs(5.)),
                    colony_module=capacity_rows(rows['colony_module'], module_factor, b, costs(2.)))
    lands = dict(round=landing(capacity['round'], flyers), round_strong=landing(capacity['round_strong'], flyers),
                 round_strong_canopy_light=landing(capacity['round_strong_canopy_light'], flyers),
                 colony_module=landing(capacity['colony_module'], flyers))
    return dict(
        evidence='Light estimated from the clear-sky surface-light product (ground spectra carried to the body\'s '
                 'height by two-stream adding through the air below, whose reflectance is scaled from the column\'s; '
                 'isotropic diffuse light; ray tests among raft neighbours); masses, trim and gas from the aerophyte '
                 'product and model; guild costs and benefits from Earth precedents and stated design guesses. No '
                 'tenant community is simulated.',
        reading_rule='Per m² of projected aerophyte area (a raft\'s hexagonal cell) unless labelled. Light is '
                     'photosynthetic photons at 10 km above the equator, a daytime mean over the lunar day (multiply by '
                     'half the 86,400 s of a day for the 24-hour mean); share_of_level_top is over a level surface on top. '
                     'Steady tenants fill the gas to 70% of the body. The sudden/at_once fields are static trim mass '
                     'equivalents; vertical_response gives their timing: the sink, period, overshoot and settling of a '
                     'fixed-volume body, and the time to vent ballonet air or drop water. Trim uses the full pressure '
                     'allowance at the assumed water-swing state; trim plus water also uses the 4 kg/m² dumpable store, '
                     'spending the damage reserve. Starting trim and control direction matter.',
        light=dict(layer_share_below_10_km=layer_share, albedos=list(ALBEDOS), cases=light, check=check),
        real_estate=real_estate(light), net_lift_kg_per_kg_h2=b,
        gas_carbon_per_kg_h2=dict(reference_light=costs(2.), canopy_light=costs(5.)),
        capacity=capacity, landing=lands, vertical_response=response_table(capacity, air10, air_profile),
        domatia=domatia_table(air10, structural, rows['colony_module']),
        guilds=guilds(flyers, food_web_chain, light),
        inoculation=inoculation(rows['round']),
        contact=contact_table(gathering, settling),
        design_guesses=dict(fringe_m2_per_m2=list(FRINGE_M2_PER_M2), founder_kg_dry_m2=list(FOUNDER_KG_DRY_M2),
                            tenant_gas_share_of_century_surplus=list(TENANT_GAS_SHARE),
                            contact_speed_m_s=list(CONTACT_SPEED_M_S), contact_margin_m=CONTACT_MARGIN_M,
                            domatium_living_kg_m2=DOMATIUM_LIVING_KG_M2, domatium_gpp=DOMATIUM_GPP,
                            domatium_food=list(DOMATIUM_FOOD), domatium_support_share=list(DOMATIUM_SUPPORT_SHARE)),
        unresolved=['Cloud below a body is the albedo 0.5 case only; the cloud field aloft is not mapped.',
                    'A raft\'s shading takes rigid touching spheres; real rafts flex and part.',
                    'Wet canopy-soil and guano masses, pathogen loads and the guards\' effect on grazing are not measured.'])
