"""The lakes' methane as a tracer of its own, carried offline through the lake model's saved columns.

The lake model (biosphere/ecology/lake.py) carries oxygen and a debt of reduced substances, whose reoxidation takes
oxygen where it returns; it carries no methane. lake_anoxia read methane as a share of the debt. Here methane is a
tracer with its own fate, run on the columns the model saved (columns.py), step by step in the model's own order:

- mixing: every step the mixed layer the model recorded mixes the tracer over its layers and the share of the next
  that lake._mix entrained, as it mixes the oxygen;
- diffusion: below it the model's own diffusivity, the mixing efficiency times the dissipation over the squared
  buoyancy frequency from the saved temperature and solute (lake.py's density table, gravity and cap), plus its
  molecular diffusivity for dissolved matter, through lake._diffuse; the temperature and solute are the day's saved
  column, mixed each step over the mixed layer as the model mixed its own, with the mixed layer at the surface
  temperature the model saved for that step;
- production in the water: in anoxic water, under 0.5 g/m3 of oxygen, the share of the water's own demand that
  respiration cannot meet, K/(K + O2), times the methane share, over 64 g of oxygen a mole (lake_anoxia as committed
  counted water to 2 g/m3, where the model reoxidises its debt in place; that rule is kept for reproducing it);
- the floors: each floor's sediment makes the methane share of its demand (the trophic state's sediment demand and
  the part of the deep water's that lake.demand_profile lays on floors), whatever the water above, since the mud is
  anoxic below its surface. Anaerobic oxidation in the mud takes what its acceptors allow; of the rest the share
  lake_anoxia's steady porewater balance gives at the floor's depth, temperature and the methane in the water above
  leaves as bubbles, each keeping to the surface the share lake_anoxia's rising bubble keeps and leaving the rest in
  the layers it crosses; what diffuses out enters anoxic water at the share its oxygen cannot meet, the rest
  oxidised at the mud's surface, as the model's Monod form meets the floor's demand;
- aerobic oxidation: methanotrophs with Michaelis-Menten kinetics in methane and oxygen, their stock growing on what
  they oxidise and lost at a rate of its own, mixed and diffused with the water, never below a seed stock;
- the oxygen: the model's oxygen already carries the reoxidation of its debt, and its debt holds the methane its
  anoxic water makes. Each methane lane carries that methane as the model holds it (made by the same rule, moved
  with the water, reoxidised with the debt at the model's rate in proportion to its part of the debt, never leaving
  the water) and the difference between the oxygen the model spends on it and the oxygen the lane's methanotrophs
  use, moved with the water and given back to the air at the surface's oxygen transfer; the methanotrophs see the
  model's oxygen plus that difference. Methane made in anoxic water that escapes, bubbles away or is oxidised
  anaerobically so returns the oxygen the model spent on it, and there the oxygen and methane budgets agree to first
  order (the model's other demands are left as the model had them). A floor under oxic water carries no such lane
  methane: what its mud bubbles off or oxidises with iron keeps the oxygen the model spent on that floor's demand,
  which leaves the water above slightly short of oxygen and its methane oxidation, limited by methane, all but
  unchanged. Where oxygen is short, the methane takes its part of what the water holds in proportion to its demand
  beside the debt's other substances;
- anaerobic oxidation in anoxic water: first order in methane, capped by its acceptor's supply (biology.py);
- exchange: from the mixed layer to the air at the model's surface-renewal transfer velocity for O2 carried to
  methane by the Schmidt numbers' ratio (Sc^-1/2), toward equilibrium with the held air's methane, as the model
  relaxes its oxygen.

Two clones check the transport against the model's own fields: the water's age (a day a day everywhere, zero in the
mixed layer), which the model carries through the same mixing and diffusion, and the debt itself (made where
respiration falls short, reoxidised at the model's rate).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from atmosphere.radiative_convective.thermodynamics import earthlike_air
from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_anoxia import run as A
from research.studies.lake_oxygen import run as R

DAY = 86400.0
YEAR = 365.25 * DAY
ANOXIC = R.ANOXIC                         # g/m3: lake_oxygen's anoxic mark; methane is made below it
PRE_FIX_MADE_BELOW = 2.0                  # g/m3: lake_anoxia as committed (schema /1) counted methane below 2 g/m3


def lake_anoxia_made_below(product=None):
    """The water lake_anoxia counts methane in: under 2 g/m3 as committed (schema /1), under its own ANOXIC since."""
    if product is not None and str(product.get('schema', '')).endswith('/1'):
        return PRE_FIX_MADE_BELOW
    return float(getattr(A, 'ANOXIC', PRE_FIX_MADE_BELOW))
HALF_SATURATION = A.HALF_SATURATION       # g/m3: the model's respiration and debt slow toward it
CH4_TO_O2 = 64.0                          # g of O2 a mole of methane takes: CH4 + 2 O2 -> CO2 + 2 H2O
SCHMIDT_RATIO = G.SCHMIDT_20C['CH4'] / L.SCHMIDT_O2_20C   # the two Schmidt numbers scale alike with temperature
AIR_CH4 = 1.87e-6                         # the held air's methane (domain synthesis), mole fraction
OVERTURN_AFTER_DAYS = 30.0                # an overturn's window runs this long past its last full mixing
N2_FRACTION = float(earthlike_air(R.DESIGN_PRESSURE_PA, R.DESIGN_CO2_PPM).fractions['N2'])
POROSITY, TORTUOSITY = 0.9, 1.0 - math.log(0.9 ** 2)      # lake_anoxia's ebullition_share defaults


@dataclass
class Lane:
    """One traced column: the oxygen column it reads and what happens to its methane.
    kind: 'methane', or the checks 'age' and 'debt'."""
    column: int
    kind: str = 'methane'
    share: float = 0.5                    # methane's share of the anoxic mineralization (lake_anoxia's central)
    # methanotrophs (biology.py gives the literature's values)
    mu20_per_h: float = 0.5 / 24.0        # maximum specific growth at 20 C
    yield_c: float = 0.3                  # mol of biomass carbon a mole of methane oxidised
    k_ch4_um: float = 5.0                 # half-saturation in methane, uM
    k_o2_um: float = 5.0                  # half-saturation in oxygen, uM
    loss_per_day: float = 0.1             # the stock's loss (grazing, lysis, settling)
    seed_c_um: float = 0.35               # the stock that stays in all water, uM of carbon
    activation_ev: float = 0.49           # temperature dependence of the uptake from 20 C
    aerobic: bool = True
    # bubbles
    ebullition: bool = True
    # anaerobic oxidation (biology.py)
    aom_floor_mol_m2_s: float = 0.0       # the mud's: the methane its acceptors' supply can oxidise, per m2 of floor
    aom_water_per_s: float = 0.0          # in anoxic water: first-order rate
    aom_water_cap_mol_m2_s: float = 0.0   # and its cap per m2 of lake surface
    label: str = ''


class Tracer:
    def __init__(self, col, lanes, bubbles=None, air_ch4=AIR_CH4, oxygen_rule='proportional',
                 mixed_layer_oxygen=True, made_below=ANOXIC, oxygen_budget=True):
        self.col, self.lanes = col, list(lanes)
        self.oxygen_rule = oxygen_rule
        self.ml_oxygen = mixed_layer_oxygen and col.ml_o2 is not None
        self.made_below, self.oxygen_budget = made_below, oxygen_budget
        b, ph = col.basin, col.physics
        self.n, self.J = b.n, len(self.lanes)
        self.v = b.volume
        self.cum_v = np.cumsum(b.volume)
        self.a0 = float(b.area[0])
        self.dt = float(col.dt)
        self.exchange_ad = b.area[1:-1] / np.diff(b.centre) * self.dt
        self.n2_factor = ph.gravity / (L.RHO0 * np.diff(b.centre))
        self.table = L._density_table(ph.salinity)
        self.osborn = ph.mixing_efficiency * ph.dissipation
        self.beta = ph.solute_density
        # the fixed model adds the debt's anaerobic products to the density (lake.Lake._debt_density)
        self.debt_dense = bool(getattr(ph, 'debt_density', 0.0) or getattr(ph, 'iron_density_per_debt', 0.0)) \
            and col.coupled_debt is not None
        self._ph_ns = SimpleNamespace(physics=ph)
        self._diffusivity = diffusivity_of(getattr(col, 'form', 'osborn'))
        v_ = b.volume
        self._ln_ns = SimpleNamespace(basin=b, physics=ph,
                                      _centre_of_volume=float((b.centre * v_).sum() / v_.sum()))
        if getattr(ph, 'deep_mixing', 'dissipation') == 'wind' and col.dissipation is None:
            raise ValueError(f'{col.name}: wind mixing needs the dissipation of every step')
        lanes = self.lanes
        kind = np.array([l.kind for l in lanes])
        self.is_m, self.is_age, self.is_debt = kind == 'methane', kind == 'age', kind == 'debt'
        self.cols = np.array([l.column for l in lanes])
        scen = [col.scenarios[c] for c in self.cols]
        sed = b.sediment
        self.demand = np.stack([L.demand_profile(b, s) for s in scen], 1) / DAY             # g O2 m-3 s-1
        floor = np.stack([floor_demand(b, s) for s in scen], 1) / DAY                        # g O2 m-2 s-1 of floor
        self.floor_g = floor * sed[:, None]                                                  # g O2 s-1 a layer
        self.water_demand = np.maximum(self.demand - self.floor_g / self.v[:, None], 0.0)   # g O2 m-3 s-1
        self.share = np.array([l.share for l in lanes])
        self.debt_rate = np.array([s.debt_rate for s in scen]) / DAY
        # the euphotic water's respiration from a labile pool where the fixed model keeps one and the run saved it
        self.labile_k = np.array([0.0 if getattr(s, 'labile_days', None) is None else 1.0 / (s.labile_days * DAY)
                                  for s in scen])
        self.pooled = bool(np.any(self.labile_k > 0.0)) and col.pool is not None
        self.y = np.array([l.yield_c for l in lanes])
        self.q20 = np.array([l.mu20_per_h / 3600.0 for l in lanes]) / self.y                # mol CH4 (mol C)-1 s-1
        self.kc = np.array([l.k_ch4_um for l in lanes]) * 1e-3                              # mol/m3
        self.ko = np.array([l.k_o2_um for l in lanes]) * 32e-3                               # g/m3 of O2
        self.loss = np.array([l.loss_per_day for l in lanes]) / DAY
        self.seed = np.array([l.seed_c_um for l in lanes]) * 1e-3                            # mol C/m3
        self.ea = np.array([l.activation_ev for l in lanes]) / L.BOLTZMANN_EV
        self.aerobic = np.array([l.aerobic for l in lanes], dtype=float)
        self.ebul = np.array([l.ebullition for l in lanes]) & self.is_m
        self.aom_floor = np.array([l.aom_floor_mol_m2_s for l in lanes])
        self.aom_k = np.array([l.aom_water_per_s for l in lanes])
        self.aom_cap = np.array([l.aom_water_cap_mol_m2_s for l in lanes]) * self.a0          # mol/s, whole lake
        self.floor_mol = self.floor_g * self.share / CH4_TO_O2                               # mol s-1 a layer
        self.k_ratio = SCHMIDT_RATIO ** -0.5
        self.bubbles = bubbles
        self.k_full, self.frac = mixed_layers(b, col.mixed_depth)
        t_s, p_atm = col.surface_t, col.pressure / G.ATM
        self.c_eq = G.solubility('CH4', t_s) * air_ch4 * (p_atm - L.vapour_pressure_atm(t_s))
        self.events = overturn_events(b, col.mixed_depth, self.dt)
        self.window = overturn_windows(b, col.mixed_depth, self.dt)
        self.event_of = np.full(col.steps, -1, dtype=np.int64)
        after = int(round(OVERTURN_AFTER_DAYS * DAY / self.dt))
        for e, (a_, b_) in enumerate(self.events):
            self.event_of[a_:min(b_ + after + 1, col.steps)] = e
        self.event_start = {a_: e for e, (a_, _) in enumerate(self.events)}
        self.snap_set = set(int(x) for x in col.snap_step)

    # --- the fields a step reads -----------------------------------------------------------------------------
    def fields(self, s):
        """A snapshot's fields for the steps it stands for: its temperature and solute as working copies that each
        step mixes over its mixed layer, its oxygen and debt as saved, which already carry the day's reoxidation
        (lake_anoxia reads them so), with the mixed layer's oxygen each step where the model saved it."""
        col = self.col
        ox = np.maximum(col.oxygen[s][:, self.cols].astype(float), 0.0)
        return dict(t=col.t[s].astype(float), sol=col.solute[s].astype(float), ox_saved=ox, ox=ox.copy(),
                    debt=np.maximum(col.debt[s][:, self.cols].astype(float), 0.0),
                    cdebt=None if col.coupled_debt is None else col.coupled_debt[s].astype(float),
                    pool=None if not self.pooled else np.maximum(col.pool[s][:, self.cols].astype(float), 0.0))

    def rates(self, w, step=None):
        """What a step needs from the working fields: the model's diffusivity (its deep mixing: the step's own
        dissipation under the wind, the Lake number from the step's friction velocity, or the fixed dissipation), the
        methane made, the share of a floor's diffusing methane that enters the water, the model's debt reoxidation,
        the temperature factor."""
        t, sol, ox, debt = w['t'], w['sol'], w['ox'], w['debt']
        ph, col = self.col.physics, self.col
        rho = np.interp(t, L._TABLE_T, self.table) + (self.beta * sol if self.beta else 0.0)
        if self.debt_dense:
            rho = rho + L.Lake._debt_density(self._ph_ns, w['cdebt'])
        n2 = np.diff(rho) * self.n2_factor
        mixing = getattr(ph, 'deep_mixing', 'dissipation')
        if mixing == 'wind':
            k_turb = self._diffusivity(n2, ph, float(col.dissipation[step]))
        elif mixing == 'lake_number':
            u2 = float(col.u_star[step]) ** 2
            k_turb = L.lake_number_diffusivity(n2, ph, L.Lake.lake_number(self._ln_ns, rho, u2, ph.gravity))
        else:
            k_turb = self._diffusivity(n2, ph)
        coeff = (k_turb + L.DIFFUSIVITY_O2) * self.exchange_ad
        unmet = HALF_SATURATION / (HALF_SATURATION + ox)
        below = ox < self.made_below
        water_demand, demand = self.water_demand, self.demand
        if self.pooled:
            water_demand = water_demand + self.labile_k[None, :] * w['pool']
            demand = demand + self.labile_k[None, :] * w['pool']
        made_water = np.where(below, water_demand * unmet, 0.0) * self.share / CH4_TO_O2            # mol m-3 s-1
        made_debt = demand * unmet * self.share / CH4_TO_O2                                            # the clone's
        enters = np.where(below, unmet, 0.0)                   # share of a floor's diffusing methane that enters
        f = 1.0 - 1.0 / (1.0 + self.dt * self.debt_rate[None, :] * ox / (HALF_SATURATION + ox))
        reox = np.minimum(debt * f, ox)                                     # g O2 m-3 the model's debt takes a step
        tfac = np.exp(self.ea[None, :] * (1.0 / 293.15 - 1.0 / (t[:, None] + 273.15)))
        return dict(t=t, coeff=coeff, ox=ox, debt=debt, made_water=made_water, made_debt=made_debt, enters=enters,
                    f=f, reox=reox, tfac=tfac, anoxic=ox < ANOXIC, k_turb=k_turb)

    # --- the run ---------------------------------------------------------------------------------------------
    def initial(self, s=0):
        """The tracer at the first snapshot: methane as lake_anoxia reads it, the debt's methane share over 64 g;
        the debt clone the debt's; the age the saved age; the methanotrophs at their seed stock; no oxygen
        difference."""
        col = self.col
        debt = np.maximum(col.debt[s][:, self.cols].astype(float), 0.0)
        c = np.where((self.is_debt | self.is_m)[None, :], self.share * debt / CH4_TO_O2, 0.0)
        c[:, self.is_age] = col.age[s].astype(float)[:, None]
        bio = np.broadcast_to(self.seed, (self.n, self.J)).copy()
        return dict(c=c, b=bio, d=np.zeros((self.n, self.J)), model=c.copy())

    def run(self, state, start=None, stop=None, record=True, profiles=False):
        """Steps start..stop-1 (default the whole record) from state (initial() or a previous run's end), which it
        carries on. Returns the budgets in moles (oxygen in g) over the whole lake, by phase: 0 between overturns, 1
        in an overturn's window, and each overturn's own."""
        col = self.col
        start = 0 if start is None else start
        stop = col.steps if stop is None else stop
        n, J = self.n, self.J
        c, bio = state['c'], state['b']
        d = state.setdefault('d', np.zeros((n, J)))
        mm = state.setdefault('model', c.copy())                 # the methane the model's debt holds
        v, cum_v, dt, a0 = self.v, self.cum_v, self.dt, self.a0
        keys = ('made_water', 'made_floor', 'floor_entered', 'floor_total', 'aom_floor', 'floor_surface_oxidised',
                'bubbles',
                'bubbles_air', 'bubbles_dissolved', 'evaded', 'invaded', 'oxidised', 'aom_water', 'o2_used_g',
                'o2_model_g', 'o2_returned_to_air_g')
        bud = {k: np.zeros((2, J)) for k in keys}
        n_ev = len(self.events)
        ev_keys = ('evaded', 'oxidised', 'bubbles_air', 'aom_water', 'made_water', 'made_floor', 'o2_used_g')
        ev_bud = {k: np.zeros((n_ev, J)) for k in ev_keys}
        ev_before = np.zeros((n_ev, J))
        snaps, snap_steps, snap_inv = [], [], []
        mean_c, mean_b = np.zeros((n, J)), np.zeros((n, J))
        evaded = np.zeros((stop - start, J)) if record else None
        s_cur, w, fe = -1, None, None
        stack = np.empty((n, 4 * J))
        m = np.flatnonzero(self.is_m)
        d_cols = np.flatnonzero(self.is_debt)
        a_cols = np.flatnonzero(self.is_age)
        bub = self.bubbles
        use_bub = bool(m.size) and bub is not None and bool(self.ebul[m].any())
        aom_w = bool(m.size) and bool(np.any(self.aom_k[m] > 0))
        sed = col.basin.sediment
        floor_cap = (self.aom_floor[None, m] * sed[:, None]) if m.size else None                  # mol/s a layer
        budget_on = 1.0 if self.oxygen_budget else 0.0
        for step in range(start, stop):
            s = col.snapshot_index(step)
            if s != s_cur:
                w, s_cur = self.fields(s), s
                if use_bub:
                    rate = np.maximum(self.floor_mol[:, m] - floor_cap, 0.0) / np.maximum(sed[:, None], 1e-300)
                    fe = bub.shares(rate, c[:, m]) * self.ebul[m]
            phase = 1 if self.window[step] else 0
            ev = int(self.event_of[step])
            if step in self.event_start:
                ev_before[self.event_start[step]] = c.T @ v
            k_full, frac = int(self.k_full[step]), float(self.frac[step])
            kk = max(1, k_full)
            # --- the working fields, the tracer, mixing and diffusion ----------------------------------
            mix(w['t'], v, cum_v, k_full, frac)
            mix(w['sol'], v, cum_v, k_full, frac)
            if w['cdebt'] is not None:
                mix(w['cdebt'], v, cum_v, k_full, frac)
            w['t'][:kk] = col.surface_t[step]                    # the mixed layer's own temperature this step
            if self.ml_oxygen:
                # the mixed layer's oxygen at this step, which the model saved step by step; the saved column below
                w['ox'][:] = w['ox_saved']
                w['ox'][:kk] = col.ml_o2[step][self.cols][None, :]
            r = self.rates(w, step)
            stack[:, :J], stack[:, J:2 * J], stack[:, 2 * J:3 * J], stack[:, 3 * J:] = c, bio, d, mm
            mix(stack, v, cum_v, k_full, frac)
            L._diffuse(stack, r['coeff'], v)
            c[:], bio[:], d[:], mm[:] = stack[:, :J], stack[:, J:2 * J], stack[:, 2 * J:3 * J], stack[:, 3 * J:]
            # --- the clones ----------------------------------------------------------------------------
            if a_cols.size:
                c[:, a_cols] += dt / DAY
            if d_cols.size:
                c[:, d_cols] += r['made_debt'][:, d_cols] * dt
                gone = np.minimum(c[:, d_cols] * r['f'][:, d_cols],
                                  r['ox'][:, d_cols] * self.share[d_cols] / CH4_TO_O2)
                c[:, d_cols] -= gone
            if m.size:
                # --- the water's methane and the floors' -----------------------------------------------
                water = r['made_water'][:, m] * v[:, None] * dt                            # mol a layer this step
                floor = self.floor_mol[:, m] * dt
                aom_f = np.minimum(floor, floor_cap * dt)
                rest = floor - aom_f
                if use_bub:
                    gas = rest * fe
                    to_air = gas * bub.survival[:, None]
                    dissolved = bub.weights.T @ gas
                else:
                    gas = to_air = dissolved = np.zeros_like(rest)
                diffusing = rest - gas
                enters = diffusing * r['enters'][:, m]
                c[:, m] += (water + enters + dissolved) / v[:, None]
                made_floor = floor * r['enters'][:, m]
                mm[:, m] += (water + made_floor) / v[:, None]            # the model's debt takes all of it
                bud['made_water'][phase, m] += water.sum(0)
                bud['made_floor'][phase, m] += made_floor.sum(0)
                bud['floor_total'][phase, m] += floor.sum(0)
                bud['floor_entered'][phase, m] += enters.sum(0)
                bud['aom_floor'][phase, m] += aom_f.sum(0)
                bud['floor_surface_oxidised'][phase, m] += (diffusing - enters).sum(0)
                bud['bubbles'][phase, m] += gas.sum(0)
                bud['bubbles_air'][phase, m] += to_air.sum(0)
                bud['bubbles_dissolved'][phase, m] += dissolved.sum(0)
                if ev >= 0:
                    ev_bud['made_water'][ev, m] += water.sum(0)
                    ev_bud['made_floor'][ev, m] += made_floor.sum(0)
                    ev_bud['bubbles_air'][ev, m] += to_air.sum(0)
                # --- aerobic oxidation by methanotrophs, on the model's oxygen and the difference --------------
                cm, bm, dm, mmm = c[:, m], bio[:, m], d[:, m], mm[:, m]
                debt_m = r['debt'][:, m]
                # the model's reoxidation of its methane this step: its part of the debt's (mol a m3)
                per_debt = np.where(debt_m > 1e-12, r['reox'][:, m] / np.maximum(debt_m, 1e-12), 0.0)
                model_ox = np.minimum(mmm * per_debt, mmm)
                mmm = mmm - model_ox
                ox_m = np.maximum(r['ox'][:, m] + budget_on * dm, 0.0)
                q = self.q20[m] * r['tfac'][:, m] * (ox_m / (self.ko[m] + ox_m)) * self.aerobic[m]
                c_new = cm / (1.0 + dt * q * bm / (self.kc[m] + cm))
                if self.oxygen_rule == 'all':
                    allowed = ox_m / CH4_TO_O2
                else:
                    need = CH4_TO_O2 * cm
                    rest_of_debt = np.maximum(debt_m - CH4_TO_O2 * mmm, 0.0)
                    allowed = ox_m / CH4_TO_O2 * need / np.maximum(need + rest_of_debt, 1e-300)
                ox_mol = np.minimum(cm - c_new, allowed)
                cm = cm - ox_mol
                bm = bm + self.y[m] * ox_mol
                bm = self.seed[m] + (bm - self.seed[m]) / (1.0 + self.loss[m] * dt)
                model_m = CH4_TO_O2 * model_ox                                     # g O2 m-3 on the model's methane
                dm = dm + budget_on * (model_m - CH4_TO_O2 * ox_mol)
                done = ox_mol.T @ v
                bud['oxidised'][phase, m] += done
                bud['o2_used_g'][phase, m] += CH4_TO_O2 * done
                bud['o2_model_g'][phase, m] += model_m.T @ v
                if ev >= 0:
                    ev_bud['oxidised'][ev, m] += done
                    ev_bud['o2_used_g'][ev, m] += CH4_TO_O2 * done
                # --- anaerobic oxidation in anoxic water ---------------------------------------------------
                if aom_w:
                    rate = np.where(r['anoxic'][:, m], self.aom_k[m] * cm, 0.0)                  # mol m-3 s-1
                    whole = rate.T @ v
                    scale = np.where(whole > self.aom_cap[m], self.aom_cap[m] / np.maximum(whole, 1e-300), 1.0)
                    took = np.minimum(rate * scale * dt, cm)
                    cm = cm - took
                    bud['aom_water'][phase, m] += took.T @ v
                    if ev >= 0:
                        ev_bud['aom_water'][ev, m] += took.T @ v
                c[:, m], bio[:, m], d[:, m], mm[:, m] = cm, bm, dm, mmm
            # --- exchange with the air -------------------------------------------------------------------
            if col.open_water[step]:
                if m.size:
                    share_v = a0 * dt / cum_v[kk - 1]
                    relax = 1.0 - math.exp(-col.k_o2[step] * self.k_ratio * share_v)
                    before = c[:kk, m].copy()
                    c[:kk, m] += (self.c_eq[step] - c[:kk, m]) * relax
                    out = (before - c[:kk, m]).T @ v[:kk]
                    bud['evaded'][phase, m] += np.maximum(out, 0.0)
                    bud['invaded'][phase, m] += np.maximum(-out, 0.0)
                    if ev >= 0:
                        ev_bud['evaded'][ev, m] += out
                    if record:
                        evaded[step - start, m] = out
                    if budget_on:
                        keep = math.exp(-col.k_o2[step] * share_v)          # the oxygen's own transfer
                        bud['o2_returned_to_air_g'][phase, m] += (d[:kk, m] * (1.0 - keep)).T @ v[:kk]
                        d[:kk, m] *= keep
                if a_cols.size:
                    c[:kk, a_cols] = 0.0
            if record:
                mean_c += c
                mean_b += bio
                if step in self.snap_set:
                    snap_inv.append(c.T @ v)
                    snap_steps.append(step)
                    if profiles:
                        snaps.append(c.astype(np.float32))
        state['c'], state['b'], state['d'], state['model'] = c, bio, d, mm
        out = dict(budget=bud, steps=(start, stop), events=dict(budget=ev_bud, inventory_before=ev_before),
                   o2_difference_g=(d.T @ v))
        if record:
            out['snap_inventory'] = np.array(snap_inv)
            out['snap_steps'] = np.array(snap_steps, dtype=np.int64)
            out['evaded_series'] = evaded
            out['mean_c'] = mean_c / max(stop - start, 1)
            out['mean_b'] = mean_b / max(stop - start, 1)
            if profiles:
                out['snap_c'] = np.array(snaps)
        return out

    def inventory(self, c):
        return c.T @ self.v

    def floor_production_max(self):
        """The largest floor production any methane lane gives, mol m-2 s-1 (for Bubbles)."""
        sed = self.col.basin.sediment
        rate = np.where(sed[:, None] > 0, self.floor_mol / np.maximum(sed[:, None], 1e-300), 0.0)
        return float(rate[:, self.is_m].max()) if self.is_m.any() else 0.0


def diffusivity_of(form):
    """The diffusivity a column's run took below its mixed layer: lake.py's Osborn rule, or Shih and Bouffard's three
    regimes of the buoyancy Reynolds number as lake_oxygen's mixing_forms runs them."""
    if form in (None, 'osborn'):
        return L.diffusivity
    if form == 'shih_bouffard':
        from research.studies.lake_oxygen import mixing_forms as MF
        return MF.shih_bouffard
    raise ValueError(f'unknown diffusivity form {form!r}')


# --- bubbles -------------------------------------------------------------------------------------------------------

class Bubbles:
    """Bubbles from each floor layer: the steady porewater balance of lake_anoxia's gases.ebullition_share, in arrays,
    and the share of each bubble's methane reaching the air from lake_anoxia's rising bubble at the layer's mean
    temperature, with the rest left in the layers it crosses. A bubble rising from depth z0 keeps at depth z the share
    s(z0)/s(z) of its methane, s being the share a bubble starting at that depth keeps to the surface (the remaining
    bubble taken as one starting there), so the layers between the floor and the surface take what it loses."""

    def __init__(self, basin, t_mean, gravity, production_max, surface_atm=A.SURFACE_ATM,
                 size_mm=A.BUBBLE_MM['central'], interface='rigid', decay_per_m=A.SEDIMENT_DECAY_PER_M,
                 n2_fraction=N2_FRACTION):
        """production_max: the largest floor production (mol m-2 s-1) any lane gives; floors deeper than where it
        can bubble with no methane in the water above are left out. The mud's N2 is the air's at saturation, its share
        n2_fraction of the dry air (the design air's by default)."""
        self.basin, self.gravity, self.surface_atm = basin, gravity, surface_atm
        self.size_mm, self.interface, self.decay = size_mm, interface, decay_per_m
        self.centre, self.z = basin.centre, basin.z
        n = basin.n
        self.t_mean = tm = np.asarray(t_mean, float)
        # the porewater balance's constants at each layer's mean temperature (gases.ebullition_share)
        p = G.local_pressure_atm(self.centre, surface_atm, gravity) - L.vapour_pressure_atm(tm)
        n2_atm = n2_fraction * (surface_atm - L.vapour_pressure_atm(tm))
        self.c_crit = G.solubility('CH4', tm) * np.maximum(p - n2_atm, 0.0)
        self.f_scale = POROSITY * G.diffusivity('CH4', tm) / TORTUOSITY * decay_per_m
        self.active = (basin.sediment > 0) & (self.f_scale * self.c_crit < production_max)
        u = np.concatenate([np.linspace(0.0, 1.0, 4001)[:-1] ** 2, np.linspace(1.0, 200.0, 20000)])
        self._u, self._r = u, 1.0 - (1.0 + u) * np.exp(-u)
        t_at = lambda depth: float(np.interp(depth, self.centre, tm))
        n2_mol = lambda tc: float(G.solubility('N2', tc)) * n2_fraction * (surface_atm
                                                                        - float(L.vapour_pressure_atm(tc)))

        def surv(depth):
            if depth <= 0.0:
                return 1.0
            tc = t_at(depth)
            return G.bubble_survival(depth, size_mm / 1000.0, tc, surface_atm, gravity, dict(N2=n2_mol(tc)),
                                     A.BUBBLE_CH4_FRACTION, interface)[0]

        last = int(np.flatnonzero(self.active).max()) if self.active.any() else -1
        s_z = np.ones(n + 1)
        s_c = np.ones(n)
        for i in range(last + 1):
            s_z[i] = surv(self.z[i])
            s_c[i] = surv(self.centre[i])
        # survival never rises with depth: interleave the interfaces and centres down the column
        if last >= 0:
            seq = np.empty(2 * (last + 1))
            seq[0::2] = s_z[:last + 1]
            seq[1::2] = s_c[:last + 1]
            seq = np.maximum(np.minimum.accumulate(seq), 1e-12)
            s_z[:last + 1], s_c[:last + 1] = seq[0::2], seq[1::2]
        self.survival = np.where(self.active, s_c, 0.0)
        wts = np.zeros((n, n))
        for i in range(last + 1):
            if not self.active[i]:
                continue
            wts[i, i] = 1.0 - s_c[i] / s_z[i]
            for k in range(i):
                wts[i, k] = s_c[i] / s_z[k + 1] - s_c[i] / s_z[k]
        self.weights = wts
        self.deepest_m = float(self.z[last + 1]) if last >= 0 else 0.0

    def shares(self, production, bottom_ch4):
        """The ebullition share of each floor layer's production (mol m-2 s-1, (n, J)) under water holding
        bottom_ch4 (mol/m3, (n, J)): exp(-u) with (1 + u) exp(-u) = 1 - F*/P (gases.ebullition_share)."""
        f_star = self.f_scale[:, None] * np.maximum(self.c_crit[:, None] - bottom_ch4, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(production > 0, f_star / production, np.inf)
        u = np.interp(np.clip(ratio, 0.0, 1.0), self._r, self._u)
        out = np.where(ratio < 1.0, np.exp(-u), 0.0)
        return np.where(self.active[:, None], out, 0.0)


# --- pieces -----------------------------------------------------------------------------------------------------

def floor_demand(basin, s):
    """Each layer's floor demand per m2 of its floor (g O2 m-2 d-1): the trophic state's sediment demand, and on the
    floors within the deep demand's reach below the euphotic zone the part of the deep water's demand that finds no
    water there, as lake.demand_profile lays it (tested to give back that profile)."""
    z_eu = min(4.6 / s.attenuation, basin.depth)
    top, bottom, centre = basin.z[:-1], basin.z[1:], basin.centre
    lower = min(z_eu + s.deep_scale_m, basin.depth)
    deep = np.clip((np.minimum(bottom, lower) - np.maximum(top, z_eu)) / basin.dz, 0.0, 1.0)
    v_deep = float((deep * basin.volume).sum())
    a_eu = float(np.interp(z_eu, basin.z, basin.area))
    volumetric = s.deep_water / s.deep_scale_m
    missing = max(0.0, s.deep_water * a_eu - volumetric * v_deep)
    floor_in_range = basin.sediment * ((centre >= z_eu) & (centre < lower))
    rate = np.full(basin.n, float(s.sediment))
    if missing > 0.0:
        weight = floor_in_range if floor_in_range.sum() > 0 else basin.sediment * (np.arange(basin.n) == basin.n - 1)
        safe = np.maximum(basin.sediment, 1e-300)
        rate = rate + np.where(basin.sediment > 0, missing * weight / weight.sum() / safe, 0.0)
    return np.where(basin.sediment > 0, rate, 0.0)


def mixed_layers(basin, depth):
    """The fully mixed layers and the entrained share of the next at each step, from the mixed layer's recorded
    depth, which lake._mix sets to the bottom of the last full layer or to the top of the next plus its entrained
    share; the saved depths are float32, so a depth within its rounding of an interface sits on it."""
    z = basin.z
    d = np.asarray(depth, float)
    i = np.clip(np.searchsorted(z, d, side='left'), 1, z.size - 1)
    tol = np.maximum(1e-4, 3e-7 * z[i])
    on = np.abs(z[i] - d) <= tol
    near_top = ~on & (np.abs(d - z[i - 1]) <= tol)
    k_full = np.where(on, i, i - 1)
    frac = np.where(on | near_top, 0.0, (d - z[i - 1]) / (z[i] - z[i - 1]))
    k_full = np.maximum(k_full, 1)
    frac = np.where(k_full >= basin.n, 0.0, np.clip(frac, 0.0, 1.0))
    return k_full.astype(np.int64), frac


def mix(x, v, cum_v, k_full, frac):
    """lake._mix's mixing of a dissolved quantity, for columns x (n, J) in place."""
    n = v.size
    v_ml = float(cum_v[k_full - 1])
    w = v[:k_full]
    if frac > 0.0 and k_full < n:
        vf = frac * v[k_full]
        mean = (w @ x[:k_full] + vf * x[k_full]) / (v_ml + vf)
        x[k_full] = frac * mean + (1.0 - frac) * x[k_full]
        x[:k_full] = mean
    else:
        x[:k_full] = (w @ x[:k_full]) / v_ml


def overturn_events(basin, mixed_depth, dt):
    """lake_anoxia's overturns: the steps at which the mixed layer reaches the bottom layer, grouped into events
    when less than 15 days apart; (first, last) step of each."""
    md = np.asarray(mixed_depth, float)
    hit = np.flatnonzero(md >= basin.z[-2] - 1e-4)
    if not hit.size:
        return []
    events, first = [], hit[0]
    for a, b in zip(hit[:-1], hit[1:]):
        if (b - a) * dt > 15 * DAY:
            events.append((int(first), int(a)))
            first = b
    events.append((int(first), int(hit[-1])))
    return events


def overturn_windows(basin, mixed_depth, dt, after_days=OVERTURN_AFTER_DAYS):
    """Steps within an overturn's window: from each event's first full mixing to after_days past its last."""
    out = np.zeros(np.asarray(mixed_depth).size, bool)
    after = int(round(after_days * DAY / dt))
    for a, b in overturn_events(basin, mixed_depth, dt):
        out[a:min(b + after + 1, out.size)] = True
    return out
