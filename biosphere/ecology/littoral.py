"""A shallow littoral community through the lunar month: producers, decomposers, grazers, crayfish and fish, with
their oxygen, carbon, phosphorus, nitrogen and calcite.

The littoral is one box per member of an ensemble: a mixed water column of depth h over a lit floor, trading a share
of its water each day with the open water beside it. Every pool is areal, per m2 of floor: carbon in g C, phosphorus
in mg P, nitrogen in mg N, oxygen in g O2, fines in g; a water-borne pool's concentration is its amount over h. The
model steps every member at once (numpy arrays over the ensemble) on its forcing's step, about an hour.

Producers. Phytoplankton grow on the column's mean light and the dissolved phosphorus, respire, die, sink to the floor
at a speed that scales with gravity (Stokes's law) and are traded with the open water. Periphyton (the floor's
microalgae and cyanobacterial mats) grow on the floor's light under the meadow's shade, limited by space. A stonewort
meadow keeps a structure and a carbohydrate reserve: its gross photosynthesis follows the floor's light through a
self-shading canopy; its upkeep runs at all hours, optionally idled to a share of itself while the floor is dark (the
land plants' carbon trigger, research/studies/lunar_cycle_ecology B1); photosynthate fills the reserve to a ceiling
and then builds structure, which needs phosphorus and nitrogen; photosynthate that nutrients cannot build into
structure past a full reserve is not fixed. A deficit draws the reserve down and then the structure. The meadow lays
calcite with its net photosynthesis at a stated molar ratio and carries phosphorus down with it. Phosphorus limits
growth outright; nitrogen fixers fill whatever the dissolved nitrogen cannot, at a stated carbon cost.

Decomposers. Detritus (carbon with its phosphorus and nitrogen) decays at a rate rising with temperature, aerobically
while the water holds oxygen and more slowly without. Methanogenesis takes a share of the anaerobic decay, turning half
its carbon to methane; a share of the methane leaves as bubbles past the methanotrophs, which oxidise a share of the
rest while oxygen lasts. A share of the nitrogen set free is lost to denitrification, with a yield of N2O. A share of
the detritus is buried with its nutrients each day.

Consumers. Each consumer is one pool or, in De Roos et al.'s (2008) stage-structured biomass form, juveniles and
adults. A stage eats every food in its diet on a multi-food Holling type II response, at a maximum rate that follows
its thermal performance and falls as oxygen drops; it assimilates each food at its own efficiency; its upkeep rises
with temperature; its net production is the least its assimilated carbon, phosphorus and nitrogen can build at its
own fixed ratios, the excess carbon respired and the excess nutrients excreted. Juveniles grow and mature at De Roos's
rate, gamma = (nu - mu) / (1 - z^(1 - mu/nu)), for net production nu, mortality mu and the ratio z of mass at birth to
mass at maturation; adults turn net production into offspring at a stated efficiency, times a cue factor for how far
the lunar month's signals let them breed. A starving stage burns its own biomass at its deficit. Juvenile prey shelter
in the meadow (their availability to predators halves at a stated meadow biomass). Heat beyond a stage's lethal limit
and oxygen below its lethal level add mortality. Crayfish clip more stonewort than they eat and stir fines from the
floor; grazers' insects emerge and leave; fishing takes a share of the harvested adults.

Oxygen. Photosynthesis adds it; respiration, decay and the oxidation of methane and reduced products take it; the air
exchanges it at the transfer velocity against saturation; the open water's exchange carries it in and out. Where the
water's oxygen tension exceeds the air's pressure plus the hydrostatic pressure of the column above the floor, the
excess leaves as bubbles: a 2 m column bubbles at 115% saturation at lunar gravity, against about 190% on Earth. An
oxygen deficit below zero stands for reduced substances that take oxygen back as it returns.

Each pool's outflows over a step are scaled down together where they would take more than it holds, with predation's
claim on a consumer's stage met first; carbon, phosphorus, nitrogen and oxygen stay conserved through every exchange
the model names, which its tests check.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from biosphere.ecology import lake as L
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

BOLTZMANN_EV = 8.617333262e-5          # eV/K
T_REF_K = 293.15                       # rates are given at 20 C
O2_PER_C = 32.0 / 12.011               # g O2 per g C respired or fixed
CACO3_PER_C = 100.087 / 12.011         # g CaCO3 per g C at a molar ratio of one
WATER_DENSITY = 998.0
WATER_POOLS = ('phyto', 'zoo')         # consumers' foods held in the water column


# ---- traits ------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Thermal:
    """A rate's response to temperature: Arrhenius from 20 C, falling above an optimum as a Sharpe-Schoolfield
    inactivation (half the rate lost at t_half_c), and a mortality above a lethal limit."""
    e_a: float = 0.65                  # eV, activation energy of the rising limb
    t_half_c: float = 34.0             # C
    e_h: float = 3.0                   # eV, steepness of the inactivation
    lethal_c: float = 37.0             # C
    lethal_rate: float = 1.0           # per day of the stage's biomass, per kelvin above the lethal limit


def arrhenius(t_c, e_a):
    return np.exp(e_a / BOLTZMANN_EV * (1.0 / T_REF_K - 1.0 / (np.asarray(t_c, float) + 273.15)))


def performance(t_c, th: Thermal):
    tk = np.asarray(t_c, float) + 273.15
    return arrhenius(t_c, th.e_a) / (1.0 + np.exp(th.e_h / BOLTZMANN_EV * (1.0 / (th.t_half_c + 273.15) - 1.0 / tk)))


def heat_mortality(t_c, th: Thermal):
    return th.lethal_rate * np.maximum(np.asarray(t_c, float) - th.lethal_c, 0.0)


@dataclass(frozen=True)
class Phyto:
    mu_max: float = 1.46               # gross growth, per day at 20 C (Eppley 1972)
    thermal: Thermal = Thermal(e_a=0.32, t_half_c=38.0, lethal_c=42.0)
    light_half: float = 60.0           # umol m-2 s-1
    p_half: float = 1.5                # mg P m-3
    respiration: float = 0.05          # per day at 20 C (e_a 0.65)
    mortality: float = 0.05            # per day
    sinking_earth: float = 0.5         # m/day under Earth's gravity
    chl_per_c: float = 0.02            # g chlorophyll per g C
    c_to_p: float = 60.0               # by mass
    c_to_n: float = 6.6


@dataclass(frozen=True)
class Periphyton:
    mu_max: float = 1.0
    thermal: Thermal = Thermal(e_a=0.32, t_half_c=38.0, lethal_c=42.0)
    light_half: float = 100.0
    p_half: float = 0.75               # mats on the floor reach its pore water
    respiration: float = 0.05
    mortality: float = 0.05
    capacity: float = 10.0             # g C m-2: space and the mat's self-shading
    c_to_p: float = 60.0
    c_to_n: float = 6.6


@dataclass(frozen=True)
class Meadow:
    p_max: float = 0.15                # gross photosynthesis, g C per g C of structure per day at 20 C, saturating light
    thermal: Thermal = Thermal(e_a=0.32, t_half_c=34.0, lethal_c=38.0, lethal_rate=0.2)
    light_half: float = 50.0
    shade: float = 0.05                # m2 per g C: the canopy's self-shading, and the shade it casts on the mat
    upkeep: float = 0.015              # per day at 20 C (e_a 0.65): a tenth of p_max
    idle_share: float = 1.0            # upkeep in the dark as a share of itself (1: no idling)
    idle_light: float = 5.0            # umol m-2 s-1 at the floor, below which a meadow idles
    reserve_ceiling: float = 0.3       # g C of reserve per g C of structure
    growth_efficiency: float = 0.75    # structure built per g C of photosynthate spent on it
    senescence: float = 0.005          # per day
    p_half: float = 0.5                # mg P m-3: rhizoids and the floor's pore water
    c_to_p: float = 200.0
    c_to_n: float = 20.0
    calcite_ratio: float = 1.0         # mol CaCO3 per mol of net photosynthesis (McConnaughey 1991)
    calcite_p: float = 0.2             # mg P carried down per g CaCO3
    stabilises: float = 20.0           # g C m-2 of meadow that halves the floor's resuspension


@dataclass(frozen=True)
class Floor:
    decay: float = 0.02                # per day at 20 C (e_a 0.65), aerobic
    anaerobic_factor: float = 0.3      # anaerobic decay's rate over aerobic
    oxygen_half: float = 1.0           # g O2 m-3
    anaerobic_floor: float = 0.25      # share of decay anaerobic even under oxic water (the sediment below its skin)
    methanogenic_share: float = 0.6    # share of anaerobic decay through methanogenesis (half its carbon to methane)
    bubble_share: float = 0.5          # share of the methane leaving as bubbles
    methanotroph_max: float = 0.9      # share of the dissolved methane oxidised while oxygen lasts
    methanotroph_o2_half: float = 0.3  # g O2 m-3
    burial: float = 0.0005             # per day of the detritus, with its nutrients
    denitrified: float = 0.3           # share of the nitrogen set free lost as N2 and N2O
    n2o_yield: float = 0.005           # N2O-N per N denitrified
    fixation_cost: float = 7.0         # g C respired per g N fixed
    wind_resuspension: float = 2.0     # g fines m-2 per day at 3 m/s on a bare floor, rising as the wind's square
    fines_sinking_earth: float = 5.0   # m/day under Earth's gravity
    kd_water: float = 0.08             # 1/m, clear water's PAR attenuation
    kd_chl: float = 0.016              # m2 per mg chlorophyll
    kd_fines: float = 0.05             # m2 per g of fines


@dataclass(frozen=True)
class Consumer:
    name: str
    structured: bool
    imax: tuple                        # (juveniles, adults): g C per g C per day at 20 C
    upkeep: tuple                      # per day at 20 C
    mortality: tuple                   # per day
    half_saturation: float             # food density for half the maximum intake (g C m-2, or g C m-3 if on_water)
    diet: tuple                        # ({food: preference} of juveniles, {food: preference} of adults)
    assimilation: dict                 # food -> efficiency
    thermal: Thermal = Thermal()
    upkeep_e_a: float = 0.65
    o2_low: float = 3.0                # g/m3: intake falls linearly below this
    o2_lethal: float = 1.0             # g/m3: below this, hypoxic mortality
    hypoxic_rate: float = 0.5          # per day below the lethal level
    c_to_p: float = 40.0
    c_to_n: float = 4.5
    z: float = 0.001                   # mass at birth over mass at maturation
    reproduction: float = 0.8          # offspring biomass per g C of adult net production
    cue: float = 1.0                   # share of reproduction the lunar month's signals allow
    on_water: bool = False             # feeds on concentrations (filter feeders), not areal densities
    vulnerability: tuple = (1.0, 1.0)  # multiplies predators' preference for its juveniles and adults
    refuge: float = 20.0               # g C m-2 of meadow that halves its juveniles' availability to predators
    clip: float = 0.0                  # g C of stonewort destroyed per g C eaten
    stir: float = 0.0                  # g fines resuspended per g C of the population per day
    emergence: float = 0.0             # per day: emerging insects leaving the water
    shell: float = 0.0                 # g CaCO3 laid per g C of net production
    harvested: bool = False            # its adults are fished
    planktonic: bool = False           # lives in the water and is traded with the open water
    open_share: float = 0.0            # its open-water carbon as a share of the open water's phytoplankton carbon


@dataclass
class Community:
    phyto: Phyto = field(default_factory=Phyto)
    periphyton: Periphyton = field(default_factory=Periphyton)
    meadow: Meadow = field(default_factory=Meadow)
    floor: Floor = field(default_factory=Floor)
    consumers: tuple = ()


@dataclass
class Setting:
    """Each member's physical setting (scalars or arrays over the ensemble)."""
    depth: object                      # m
    exchange: object                   # per day
    tp_open: object                    # mg P m-3 in the open water
    dip_share_open: float = 0.3        # dissolved inorganic share of the open water's phosphorus
    phyto_share_open: float = 0.5      # share in phytoplankton
    din_open: object = 20.0            # mg N m-3
    fines_open: object = 0.3           # g m-3
    gravity: float = MOON_SURFACE_GRAVITY
    pressure_pa: float = 121590.0
    o2_fraction: float = 0.1746
    salinity: float = 0.142
    harvest: object = 0.0              # per year, on harvested adults, from harvest_from_year
    harvest_from_year: float = 1e9


def _arr(x, n):
    a = np.asarray(x, dtype=float)
    return np.full(n, float(a)) if a.ndim == 0 else a.astype(float).copy()


# The accounts a run collects (per m2, summed over its steps); a key appears once its flux has occurred.
ACCOUNTS = ('gpp_c', 'resp_c', 'calcite_meadow', 'calcite_shell', 'p_calcite', 'burial_c', 'burial_p', 'burial_n',
            'emerged_c', 'emerged_p', 'emerged_n', 'harvest_c', 'harvest_p', 'harvest_n', 'ch4_c', 'ch4_bubbles_c',
            'ch4_oxidised_c', 'n2o_n', 'denitrified_n', 'n_fixed', 'o2_air', 'o2_bubbles', 'o2_exchange', 'o2_bio',
            'p_exchange', 'n_exchange', 'c_exchange', 'births_c', 'eaten_chara_c', 'clipped_chara_c')


class Littoral:
    """The community model for an ensemble of n members sharing one Community (whose traits may be scalars or arrays
    over the members) and one Setting."""

    def __init__(self, community: Community, setting: Setting, n: int):
        self.c, self.s, self.n = community, setting, n
        self.cons = list(community.consumers)
        self.stages = []                 # (pool name, consumer index, stage: 0 juveniles, 1 adults or unstructured)
        for k, con in enumerate(self.cons):
            if con.structured:
                self.stages += [(con.name + '_j', k, 0), (con.name + '_a', k, 1)]
            else:
                self.stages += [(con.name, k, 1)]
        self.stage_of = {name: (k, st) for name, k, st in self.stages}
        self.depth = _arr(setting.depth, n)
        self.exchange = _arr(setting.exchange, n)
        self.tp_open = _arr(setting.tp_open, n)
        self.din_open = _arr(setting.din_open, n)
        self.fines_open = _arr(setting.fines_open, n)
        self.harvest = _arr(setting.harvest, n)
        g = setting.gravity / STANDARD_GRAVITY
        self.sink_phyto = community.phyto.sinking_earth * g
        self.sink_fines = community.floor.fines_sinking_earth * g
        # the bubble point: oxygen's tension may exceed saturation by the hydrostatic pressure above the floor
        self.bubble_ratio = 1.0 + WATER_DENSITY * setting.gravity * self.depth / (setting.o2_fraction * setting.pressure_pa)

    # ---- pools ----------------------------------------------------------------------------------------------------------
    def initial(self, t_c=29.0):
        n, s, h = self.n, self.s, self.depth
        state = dict(dip=self.tp_open * s.dip_share_open * h, din=self.din_open * h,
                     oxygen=L.oxygen_saturation(t_c, s.pressure_pa / 101325.0, s.salinity, s.o2_fraction) * h,
                     fines=self.fines_open * h, phyto=np.zeros(n), peri=np.zeros(n), chara=np.zeros(n),
                     reserve=np.zeros(n), det_c=np.zeros(n), det_p=np.zeros(n), det_n=np.zeros(n))
        for name, _, _ in self.stages:
            state[name] = np.zeros(n)
        return state

    def seed(self, state, pool, amount):
        """Adds a seed biomass (g C m-2) to a living pool, or bare carbon to the detritus (det_c); a meadow's seed
        carries a tenth of its mass as reserve. Returns the phosphorus and nitrogen it brings (mg m-2)."""
        a = _arr(amount, self.n)
        state[pool] += a
        if pool == 'det_c':                                   # bare organic carbon, with no nutrients
            return np.zeros(self.n), np.zeros(self.n)
        cp, cn = self.ratios(pool)
        if pool == 'chara':
            state['reserve'] += 0.1 * a
        return a * 1000.0 / cp, a * 1000.0 / cn

    def ratios(self, pool):
        """(C:P, C:N) by mass of a living pool."""
        c = self.c
        if pool == 'phyto':
            return c.phyto.c_to_p, c.phyto.c_to_n
        if pool == 'peri':
            return c.periphyton.c_to_p, c.periphyton.c_to_n
        if pool == 'chara':
            return c.meadow.c_to_p, c.meadow.c_to_n
        con = self.cons[self.stage_of[pool][0]]
        return con.c_to_p, con.c_to_n

    def living(self):
        return ['phyto', 'peri', 'chara'] + [name for name, _, _ in self.stages]

    def totals(self, x):
        """Total organic carbon (g), phosphorus (mg) and nitrogen (mg) per m2 in the box, and its oxygen (g)."""
        c = x['det_c'] + x['reserve']
        p = x['dip'] + x['det_p']
        nn = x['din'] + x['det_n']
        for pool in self.living():
            cp, cn = self.ratios(pool)
            c = c + x[pool]
            p = p + x[pool] * 1000.0 / cp
            nn = nn + x[pool] * 1000.0 / cn
        return c, p, nn, x['oxygen']

    # ---- one step -------------------------------------------------------------------------------------------------------
    def step(self, x, f, dt, year, acc):
        """Advance every member by dt days. f holds the step's forcing per member: t (C), t_open (C), par0
        (umol m-2 s-1 just below the surface), k (O2 transfer velocity, m/s) and wind (m/s at 10 m). acc collects the
        step's fluxes per m2. Returns the step's diagnostics."""
        c, s, n, h = self.c, self.s, self.n, self.depth
        ph, pe, me, fl = c.phyto, c.periphyton, c.meadow, c.floor
        t = np.asarray(f['t'], float)
        zero = np.zeros(n)
        add = lambda k, v: acc.__setitem__(k, acc.get(k, 0.0) + v)

        # the setting
        chl = x['phyto'] / h * 1000.0 * ph.chl_per_c                        # mg chl m-3
        kd = fl.kd_water + fl.kd_chl * chl + fl.kd_fines * x['fines'] / h
        floor_light = f['par0'] * np.exp(-kd * h)
        column_light = f['par0'] * (1.0 - np.exp(-kd * h)) / np.maximum(kd * h, 1e-9)
        sat = L.oxygen_saturation(t, s.pressure_pa / 101325.0, s.salinity, s.o2_fraction)
        o2c = np.maximum(x['oxygen'] / h, 0.0)
        dipc = np.maximum(x['dip'] / h, 0.0)
        warm = arrhenius(t, 0.65)

        out = {k: zero.copy() for k in x}                                     # every pool's outflow over the step
        inn = {k: zero.copy() for k in x}                                     # and inflow
        o2_made, o2_used = zero.copy(), zero.copy()
        resp_c = zero.copy()

        # ---- the plankton and the mat: potential growth (gross less respiration) ---------------------------------
        gross_a = ph.mu_max * performance(t, ph.thermal) * column_light / (column_light + ph.light_half) * x['phyto'] * dt
        gross_a *= dipc / (dipc + ph.p_half)
        resp_a = ph.respiration * warm * x['phyto'] * dt
        mort_a = (ph.mortality + heat_mortality(t, ph.thermal)) * x['phyto'] * dt
        sink_a = self.sink_phyto / h * x['phyto'] * dt
        shade = np.exp(-me.shade * x['chara'])
        peri_light = floor_light * shade
        gross_b = (pe.mu_max * performance(t, pe.thermal) * peri_light / (peri_light + pe.light_half) * x['peri']
                   * np.maximum(1.0 - x['peri'] / pe.capacity, 0.0) * dipc / (dipc + pe.p_half) * dt)
        resp_b = pe.respiration * warm * x['peri'] * dt
        mort_b = (pe.mortality + heat_mortality(t, pe.thermal)) * x['peri'] * dt

        # ---- the meadow's own balance --------------------------------------------------------------------------
        canopy = (1.0 - np.exp(-me.shade * x['chara'])) / me.shade
        gm = me.p_max * performance(t, me.thermal) * canopy * floor_light / (floor_light + me.light_half) * dt
        dark = floor_light < me.idle_light
        um = me.upkeep * warm * x['chara'] * np.where(dark, me.idle_share, 1.0) * dt
        reserve = x['reserve'] + gm - um
        deficit = np.maximum(-reserve, 0.0)
        reserve = np.maximum(reserve, 0.0)
        sen = (me.senescence + heat_mortality(t, me.thermal)) * x['chara'] * dt
        catabolised = np.minimum(deficit, np.maximum(x['chara'] - sen, 0.0))
        um = um - (deficit - catabolised)                                     # upkeep not met is not respired
        ceiling = me.reserve_ceiling * x['chara']
        spare = np.maximum(reserve - ceiling, 0.0)
        potential_build = spare * me.growth_efficiency

        # ---- phosphorus shared by the producers' growth ------------------------------------------------------------
        grow_a = gross_a - resp_a
        grow_b = gross_b - resp_b
        need_p = (np.maximum(grow_a, 0.0) * 1000.0 / ph.c_to_p + np.maximum(grow_b, 0.0) * 1000.0 / pe.c_to_p
                  + potential_build * dipc / (dipc + me.p_half) * 1000.0 / me.c_to_p)
        have_p = np.maximum(x['dip'], 0.0)
        p_scale = np.where(need_p > have_p, have_p / np.maximum(need_p, 1e-30), 1.0)
        grow_a = np.where(grow_a > 0, grow_a * p_scale, grow_a)
        grow_b = np.where(grow_b > 0, grow_b * p_scale, grow_b)
        build = potential_build * dipc / (dipc + me.p_half) * p_scale
        spent = build / me.growth_efficiency
        reserve = reserve - spent
        excess = np.maximum(reserve - ceiling, 0.0)                           # photosynthate nothing could take
        reserve = reserve - excess
        gm = gm - excess
        # the plankton's and mat's losses to themselves: a negative growth burns their own biomass
        loss_a_self = np.maximum(-grow_a, 0.0) + mort_a + sink_a
        loss_b_self = np.maximum(-grow_b, 0.0) + mort_b

        # ---- consumers' intake ---------------------------------------------------------------------------------
        intake = {}
        prey_out = {k: zero.copy() for k in x}
        clip_out = zero.copy()
        for name, k, stage in self.stages:
            con = self.cons[k]
            diet = con.diet[min(stage, len(con.diet) - 1)]
            dens = {}
            for food, pref in diet.items():
                pool = 'det_c' if food == 'detritus' else food
                amount = x[pool]
                if pool in self.stage_of:
                    kk, st = self.stage_of[pool]
                    prey = self.cons[kk]
                    amount = amount * prey.vulnerability[st]
                    if st == 0:
                        amount = amount / (1.0 + x['chara'] / prey.refuge)
                if con.on_water:
                    amount = amount / h if pool in WATER_POOLS else zero
                dens[food] = pref * amount
            total = sum(dens.values()) if dens else zero
            o2f = np.clip((o2c - con.o2_lethal) / max(con.o2_low - con.o2_lethal, 1e-9), 0.0, 1.0)
            rate = con.imax[min(stage, len(con.imax) - 1)] * performance(t, con.thermal) * o2f * x[name] * dt
            intake[name] = {food: rate * d / (con.half_saturation + total) for food, d in dens.items()}
            for food, amt in intake[name].items():
                pool = 'det_c' if food == 'detritus' else food
                prey_out[pool] = prey_out[pool] + amt
                if pool == 'chara':
                    clip_out = clip_out + amt * con.clip

        # ---- the floor's decay and burial ------------------------------------------------------------------------
        aerobic = (1.0 - fl.anaerobic_floor) * o2c / (o2c + fl.oxygen_half)
        decay = x['det_c'] * fl.decay * warm * (aerobic + (1.0 - aerobic) * fl.anaerobic_factor) * dt
        anaerobic_c = x['det_c'] * fl.decay * warm * (1.0 - aerobic) * fl.anaerobic_factor * dt
        burial = x['det_c'] * fl.burial * dt

        # ---- scale every pool's outflows to what it holds --------------------------------------------------------
        def fit(have, *flows):
            total = sum(flows)
            k = np.where(total > have, np.maximum(have, 0.0) / np.maximum(total, 1e-30), 1.0)
            return k
        k_phyto = fit(x['phyto'], loss_a_self, prey_out['phyto'])
        k_peri = fit(x['peri'], loss_b_self, prey_out['peri'])
        k_chara = fit(np.maximum(x['chara'] - sen - catabolised, 0.0), prey_out['chara'], clip_out)
        k_det = fit(x['det_c'], decay, burial, prey_out['det_c'])
        scale = dict(phyto=k_phyto, peri=k_peri, chara=k_chara, det_c=k_det)
        for name, _, _ in self.stages:
            scale[name] = fit(x[name], prey_out[name])                        # predation has first claim
        for name in intake:
            for food in intake[name]:
                pool = 'det_c' if food == 'detritus' else food
                intake[name][food] = intake[name][food] * scale[pool]
        loss_a_self, gr_a_neg = loss_a_self * k_phyto, np.maximum(-grow_a, 0.0) * k_phyto
        mort_a, sink_a = mort_a * k_phyto, sink_a * k_phyto
        loss_b_self, gr_b_neg = loss_b_self * k_peri, np.maximum(-grow_b, 0.0) * k_peri
        mort_b = mort_b * k_peri
        clip_out = clip_out * k_chara
        decay, anaerobic_c, burial = decay * k_det, anaerobic_c * k_det, burial * k_det

        # ---- producers: carbon, nutrients, oxygen --------------------------------------------------------------
        det_c_in, det_p_in, det_n_in = zero.copy(), zero.copy(), zero.copy()

        def to_detritus(carbon, cp, cn):
            nonlocal det_c_in, det_p_in, det_n_in
            det_c_in = det_c_in + carbon
            det_p_in = det_p_in + carbon * 1000.0 / cp
            det_n_in = det_n_in + carbon * 1000.0 / cn

        dip_in, dip_out, din_in, din_out = zero.copy(), zero.copy(), zero.copy(), zero.copy()
        pos_a, pos_b = np.maximum(grow_a, 0.0), np.maximum(grow_b, 0.0)
        inn['phyto'] += pos_a
        out['phyto'] += gr_a_neg + mort_a + sink_a
        inn['peri'] += pos_b
        out['peri'] += gr_b_neg + mort_b
        to_detritus(mort_a + sink_a, ph.c_to_p, ph.c_to_n)
        to_detritus(mort_b, pe.c_to_p, pe.c_to_n)
        dip_in += gr_a_neg * 1000.0 / ph.c_to_p + gr_b_neg * 1000.0 / pe.c_to_p
        din_in += gr_a_neg * 1000.0 / ph.c_to_n + gr_b_neg * 1000.0 / pe.c_to_n
        gpp = (pos_a + np.minimum(resp_a, gross_a) + pos_b + np.minimum(resp_b, gross_b)) + gm
        resp_prod = np.minimum(resp_a, gross_a) + gr_a_neg + np.minimum(resp_b, gross_b) + gr_b_neg + um + (spent - build)
        dip_out += pos_a * 1000.0 / ph.c_to_p + pos_b * 1000.0 / pe.c_to_p + build * 1000.0 / me.c_to_p
        n_need = pos_a * 1000.0 / ph.c_to_n + pos_b * 1000.0 / pe.c_to_n + build * 1000.0 / me.c_to_n
        n_take = np.minimum(n_need, np.maximum(x['din'], 0.0))
        n_fixed = n_need - n_take
        din_out += n_take
        fix_cost = n_fixed / 1000.0 * fl.fixation_cost
        # the fixation's carbon comes out of the producers' new growth, pro rata
        grown = pos_a + pos_b + build
        fix_share = np.where(grown > 0, np.minimum(fix_cost / np.maximum(grown, 1e-30), 1.0), 0.0)
        fix_cost = fix_share * grown
        inn['phyto'] -= pos_a * fix_share
        inn['peri'] -= pos_b * fix_share
        build_kept = build * (1.0 - fix_share)
        # nutrients of the growth the fixation's cost took back are released again
        dip_in += (pos_a * fix_share * 1000.0 / ph.c_to_p + pos_b * fix_share * 1000.0 / pe.c_to_p
                   + build * fix_share * 1000.0 / me.c_to_p)
        din_in += (pos_a * fix_share * 1000.0 / ph.c_to_n + pos_b * fix_share * 1000.0 / pe.c_to_n
                   + build * fix_share * 1000.0 / me.c_to_n)
        resp_prod = resp_prod + fix_cost
        # the meadow: structure gains what was built, loses senescence, catabolism, grazing and clipping
        inn['chara'] += build_kept
        out['chara'] += sen + catabolised
        to_detritus(sen, me.c_to_p, me.c_to_n)
        dip_in += catabolised * 1000.0 / me.c_to_p
        din_in += catabolised * 1000.0 / me.c_to_n
        new_reserve = reserve
        o2_made += gpp * O2_PER_C
        o2_used += resp_prod * O2_PER_C
        resp_c += resp_prod
        add('gpp_c', gpp)
        # calcite with the meadow's net photosynthesis, carrying phosphorus down
        calcite = np.maximum(gm - um, 0.0) * me.calcite_ratio * CACO3_PER_C
        p_calc = np.minimum(calcite * me.calcite_p, np.maximum(x['dip'] - dip_out + dip_in, 0.0))
        dip_out += p_calc
        add('calcite_meadow', calcite)
        add('p_calcite', p_calc)

        # ---- consumers ------------------------------------------------------------------------------------------
        built, burned = {}, {}
        shell = zero.copy()
        for name, k, stage in self.stages:
            con = self.cons[k]
            a_c, a_p, a_n = zero.copy(), zero.copy(), zero.copy()
            for food, amt in intake[name].items():
                pool = 'det_c' if food == 'detritus' else food
                out[pool] = out[pool] + amt
                eff = con.assimilation.get(food, 0.0)
                if pool == 'det_c':
                    cp_f = x['det_c'] * 1000.0 / np.maximum(x['det_p'], 1e-30)
                    cn_f = x['det_c'] * 1000.0 / np.maximum(x['det_n'], 1e-30)
                    cp_f = np.where(x['det_p'] > 0, cp_f, 1e30)
                    cn_f = np.where(x['det_n'] > 0, cn_f, 1e30)
                else:
                    cp_f, cn_f = self.ratios(pool)
                a_c += eff * amt
                a_p += eff * amt * 1000.0 / cp_f
                a_n += eff * amt * 1000.0 / cn_f
                det_c_in += (1.0 - eff) * amt
                det_p_in += (1.0 - eff) * amt * 1000.0 / cp_f
                det_n_in += (1.0 - eff) * amt * 1000.0 / cn_f
                if pool == 'chara':
                    add('eaten_chara_c', amt)
            up = con.upkeep[min(stage, len(con.upkeep) - 1)] * arrhenius(t, con.upkeep_e_a) * x[name] * dt
            surplus = a_c - up
            b_c = np.minimum(np.maximum(surplus, 0.0), np.minimum(a_p * con.c_to_p, a_n * con.c_to_n) / 1000.0)
            # respiration: upkeep met from intake, and the surplus that nutrients could not build
            r = np.minimum(a_c, up) + np.maximum(surplus, 0.0) - b_c
            o2_used += r * O2_PER_C
            resp_c += r
            dip_in += a_p - b_c * 1000.0 / con.c_to_p
            din_in += a_n - b_c * 1000.0 / con.c_to_n
            built[name] = b_c
            burned[name] = np.maximum(-surplus, 0.0)                         # the deficit, burned from biomass
            shell += b_c * con.shell
        # stage losses within what predation left; maturation, births, emergence, harvest
        emerged = zero.copy()
        harvested = zero.copy()
        births = zero.copy()
        for name, k, stage in self.stages:
            con = self.cons[k]
            left = np.maximum(x[name] - prey_out[name] * scale[name], 0.0)
            mu = (con.mortality[min(stage, len(con.mortality) - 1)] + heat_mortality(t, con.thermal)
                  + np.where(o2c < con.o2_lethal, con.hypoxic_rate, 0.0))
            dead = x[name] * (1.0 - np.exp(-mu * dt))
            starve = burned[name]
            em = con.emergence * x[name] * dt
            hv = (self.harvest / 365.25 * x[name] * dt if (con.harvested and stage == 1 and year >= s.harvest_from_year)
                  else zero)
            k_left = fit(left, dead, starve, em, hv)
            dead, starve, em, hv = dead * k_left, starve * k_left, em * k_left, hv * k_left
            out[name] += dead + starve + em + hv
            to_detritus(dead, con.c_to_p, con.c_to_n)
            o2_used += starve * O2_PER_C
            resp_c += starve
            dip_in += starve * 1000.0 / con.c_to_p
            din_in += starve * 1000.0 / con.c_to_n
            emerged += em
            add('emerged_p', em * 1000.0 / con.c_to_p)
            add('emerged_n', em * 1000.0 / con.c_to_n)
            harvested += hv
            add('harvest_p', hv * 1000.0 / con.c_to_p)
            add('harvest_n', hv * 1000.0 / con.c_to_n)
            if not con.structured:
                inn[name] += built[name]
                continue
            if stage == 0:
                rest = np.maximum(left - dead - starve - em - hv, 0.0)
                nu = built[name] / np.maximum(x[name], 1e-30) / dt
                mu_all = mu + (prey_out[name] * scale[name] + starve) / np.maximum(x[name], 1e-30) / dt
                mature = np.minimum(maturation(nu, mu_all, con.z) * x[name] * dt, rest)
                inn[name] += built[name]
                out[name] += mature
                inn[name[:-2] + '_a'] += mature
            else:
                b = built[name] * con.reproduction * con.cue
                inn[name[:-2] + '_j'] += b
                births += b
                lost = built[name] - b                                         # conversion losses, respired
                o2_used += lost * O2_PER_C
                resp_c += lost
                dip_in += lost * 1000.0 / con.c_to_p
                din_in += lost * 1000.0 / con.c_to_n
        add('emerged_c', emerged)
        add('harvest_c', harvested)
        add('births_c', births)
        add('calcite_shell', shell)
        out['chara'] += clip_out
        to_detritus(clip_out, me.c_to_p, me.c_to_n)
        add('clipped_chara_c', clip_out)
        # structure eaten, clipped or senescent takes its share of the reserve to the detritus (carbohydrate only)
        lost = np.clip((out['chara'] - catabolised) / np.maximum(x['chara'], 1e-30), 0.0, 1.0)
        reserve_lost = new_reserve * lost
        new_reserve = new_reserve - reserve_lost
        det_c_in += reserve_lost

        # ---- the floor: decay's products, methane, nitrogen, burial ------------------------------------------------
        frac = decay / np.maximum(x['det_c'], 1e-30)
        bfrac = burial / np.maximum(x['det_c'], 1e-30)
        efrac = out['det_c'] / np.maximum(x['det_c'], 1e-30)                  # eaten so far
        p_rel, n_rel = x['det_p'] * frac, x['det_n'] * frac
        bur_p, bur_n = x['det_p'] * bfrac, x['det_n'] * bfrac
        det_p_eaten, det_n_eaten = x['det_p'] * efrac, x['det_n'] * efrac
        out['det_c'] += decay + burial
        inn['det_c'] += det_c_in
        out['det_p'] += p_rel + bur_p + det_p_eaten
        inn['det_p'] += det_p_in
        out['det_n'] += n_rel + bur_n + det_n_eaten
        inn['det_n'] += det_n_in
        dip_in += p_rel
        denit = fl.denitrified * n_rel
        din_in += n_rel - denit
        methanogenic = anaerobic_c * fl.methanogenic_share
        methane = 0.5 * methanogenic
        bubbled = methane * fl.bubble_share
        oxidised = (methane - bubbled) * fl.methanotroph_max * o2c / (o2c + fl.methanotroph_o2_half)
        emitted = methane - oxidised
        # oxygen of decay: every carbon oxidised to CO2 takes O2, methanogenesis's own CO2 takes none, and methane
        # oxidised takes two O2 for each carbon
        o2_used += (decay - methanogenic) * O2_PER_C + oxidised * 2.0 * O2_PER_C
        resp_c += decay - emitted
        add('burial_c', burial)
        add('burial_p', bur_p)
        add('burial_n', bur_n)
        add('ch4_c', emitted)
        add('ch4_bubbles_c', bubbled)
        add('ch4_oxidised_c', oxidised)
        add('denitrified_n', denit)
        add('n2o_n', denit * fl.n2o_yield)
        add('n_fixed', n_fixed)

        # ---- apply the biology ----------------------------------------------------------------------------------
        new = {}
        for k_ in x:
            new[k_] = x[k_] + inn[k_] - out[k_]
        new['reserve'] = new_reserve
        new['dip'] = x['dip'] + dip_in - dip_out
        new['din'] = x['din'] + din_in - din_out
        oxy = x['oxygen'] + o2_made - o2_used
        add('o2_bio', o2_made - o2_used)
        add('resp_c', resp_c)

        # ---- the air, bubbles, fines and the open water ------------------------------------------------------------
        sat_amt = sat * h
        air = (sat_amt - oxy) * (1.0 - np.exp(-f['k'] * 86400.0 / h * dt))
        oxy = oxy + air
        bubbles = np.maximum(oxy - sat_amt * self.bubble_ratio, 0.0)
        oxy = oxy - bubbles
        hold = 1.0 / (1.0 + x['chara'] / me.stabilises)
        stir = sum((self.cons[k].stir * x[name] for name, k, _ in self.stages), zero)
        raised = (fl.wind_resuspension * (np.asarray(f['wind'], float) / 3.0) ** 2 + stir) * hold * dt
        settled = x['fines'] * (1.0 - np.exp(-self.sink_fines / h * dt))
        new['fines'] = x['fines'] + raised - settled
        ex = 1.0 - np.exp(-self.exchange * dt)
        o2_open = L.oxygen_saturation(f['t_open'], s.pressure_pa / 101325.0, s.salinity, s.o2_fraction) * h
        ex_o2 = (o2_open - oxy) * ex
        new['oxygen'] = oxy + ex_o2
        phyto_open = self.tp_open * s.phyto_share_open * ph.c_to_p / 1000.0 * h
        ex_dip = (self.tp_open * s.dip_share_open * h - new['dip']) * ex
        ex_din = (self.din_open * h - new['din']) * ex
        ex_phyto = (phyto_open - new['phyto']) * ex
        new['dip'] += ex_dip
        new['din'] += ex_din
        new['phyto'] += ex_phyto
        new['fines'] += (self.fines_open * h - new['fines']) * ex
        p_ex = ex_dip + ex_phyto * 1000.0 / ph.c_to_p
        n_ex = ex_din + ex_phyto * 1000.0 / ph.c_to_n
        c_ex = ex_phyto
        for name, k, _ in self.stages:
            con = self.cons[k]
            if con.planktonic:
                z_open = con.open_share * phyto_open
                ex_z = (z_open - new[name]) * ex
                new[name] += ex_z
                p_ex = p_ex + ex_z * 1000.0 / con.c_to_p
                n_ex = n_ex + ex_z * 1000.0 / con.c_to_n
                c_ex = c_ex + ex_z
        add('o2_air', air)
        add('o2_bubbles', bubbles)
        add('o2_exchange', ex_o2)
        add('p_exchange', p_ex)
        add('n_exchange', n_ex)
        add('c_exchange', c_ex)
        for k_ in new:
            if k_ != 'oxygen':
                new[k_] = np.maximum(new[k_], 0.0)
        x.update(new)
        return dict(floor_light=floor_light, o2=np.maximum(x['oxygen'], 0.0) / h, kd=kd, saturation=sat)


def maturation(nu, mu, z):
    """De Roos et al.'s (2008) maturation rate (per day) of a juvenile stage with net production nu and mortality mu
    (per day), z the ratio of mass at birth to mass at maturation; zero where nu is not positive."""
    nu = np.asarray(nu, float)
    mu = np.asarray(mu, float) * np.ones_like(nu)
    pos = nu > 1e-12
    r = np.where(pos, 1.0 - mu / np.where(pos, nu, 1.0), 0.0)
    near = pos & (np.abs(r) < 1e-6)
    far = pos & ~near
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        g = np.where(far, (nu - mu) / (1.0 - np.power(z, np.where(far, r, 1.0))), 0.0)
    g = np.where(near, nu / -math.log(z), g)
    return np.maximum(np.nan_to_num(g, nan=0.0, posinf=0.0), 0.0)


# ---- driving the model ---------------------------------------------------------------------------------------------------
def drive(model: Littoral, state, forcing, dt, steps, schedule=(), record_every=24, start_year=0.0):
    """Step the model through its forcing and return what it did.

    forcing maps t, t_open, par0, k and wind to arrays over the steps (one value per step, or one row per step with a
    column per member); series shorter than the run repeat from their start. schedule lists (year, pool, g C m-2)
    seedings, each made at the first step at or after its year. Returns the recorded times (days), each pool's series
    and the oxygen's (g/m3: the value at each record and the lowest within the interval before it), the accounts
    summed over the run, and the carbon, phosphorus and nitrogen the seeds brought.
    """
    names = list(state)
    n_rec = steps // record_every
    series = {k: np.zeros((n_rec, model.n)) for k in names}
    o2, o2_low = np.zeros((n_rec, model.n)), np.zeros((n_rec, model.n))
    floor_light = np.zeros((n_rec, model.n))
    time = np.zeros(n_rec)
    acc = {}
    seeds = dict(c=np.zeros(model.n), p=np.zeros(model.n), n=np.zeros(model.n))
    pending = sorted(schedule, key=lambda e: e[0])
    arrays = {k: np.asarray(v, dtype=float) for k, v in forcing.items()}
    low = np.full(model.n, np.inf)
    light_sum = np.zeros(model.n)
    r = 0
    for i in range(steps):
        year = start_year + i * dt / 365.25
        while pending and pending[0][0] <= year:
            _, pool, amount = pending.pop(0)
            p_in, n_in = model.seed(state, pool, amount)
            a = _arr(amount, model.n)
            seeds['c'] += a * (1.1 if pool == 'chara' else 1.0)
            seeds['p'] += p_in
            seeds['n'] += n_in
        f = {k: v[i % v.shape[0]] for k, v in arrays.items()}
        d = model.step(state, f, dt, year, acc)
        low = np.minimum(low, d['o2'])
        light_sum += d['floor_light']
        if (i + 1) % record_every == 0 and r < n_rec:
            time[r] = (i + 1) * dt
            for k in names:
                series[k][r] = state[k]
            o2[r], o2_low[r] = d['o2'], low
            floor_light[r] = light_sum / record_every
            low = np.full(model.n, np.inf)
            light_sum = np.zeros(model.n)
            r += 1
    return dict(time_day=time, series=series, o2=o2, o2_low=o2_low, floor_light=floor_light, accounts=acc, seeds=seeds)


def example_community(idle_share=1.0):
    """An illustrative community for the tests and for trying the model: plankton, mat and meadow with a zooplankter,
    a grazing snail, a crayfish and a predatory fish. Every trait is a round placeholder of a plausible order of
    magnitude, chosen so that the pieces interact; none is a species' measured value or a design."""
    warm = Thermal(e_a=0.65, t_half_c=36.0, lethal_c=38.0)
    zoo = Consumer('zoo', False, imax=(1.2,), upkeep=(0.12,), mortality=(0.05,), half_saturation=0.25,
                   diet=({'phyto': 1.0},), assimilation={'phyto': 0.6}, thermal=warm, on_water=True, planktonic=True,
                   open_share=0.1)
    snail = Consumer('snail', False, imax=(0.4,), upkeep=(0.015,), mortality=(0.01,), half_saturation=1.0,
                     diet=({'peri': 1.0, 'detritus': 0.3},), assimilation={'peri': 0.6, 'detritus': 0.2},
                     thermal=warm, shell=3.0, c_to_p=60.0, c_to_n=5.0)
    cray = Consumer('cray', True, imax=(0.25, 0.08), upkeep=(0.012, 0.005), mortality=(0.01, 0.003),
                    half_saturation=1.5,
                    diet=({'peri': 1.0, 'detritus': 0.5, 'snail': 0.5},
                          {'chara': 0.5, 'detritus': 1.0, 'snail': 1.0, 'peri': 0.5, 'cray_j': 0.2}),
                    assimilation={'peri': 0.6, 'detritus': 0.3, 'snail': 0.7, 'chara': 0.4, 'cray_j': 0.7},
                    thermal=Thermal(t_half_c=33.0, lethal_c=36.0), z=0.001, clip=1.5, stir=0.5, harvested=True,
                    vulnerability=(1.0, 0.3))
    fish = Consumer('fish', True, imax=(0.2, 0.05), upkeep=(0.012, 0.005), mortality=(0.01, 0.002),
                    half_saturation=0.5,
                    diet=({'zoo': 1.0, 'snail': 0.3}, {'cray_j': 1.0, 'cray_a': 0.3, 'snail': 0.3, 'fish_j': 0.2}),
                    assimilation={'zoo': 0.7, 'snail': 0.7, 'cray_j': 0.7, 'cray_a': 0.7, 'fish_j': 0.7},
                    thermal=warm, z=0.0001, vulnerability=(1.0, 0.0))
    return Community(meadow=Meadow(idle_share=idle_share), consumers=(zoo, snail, cray, fish))


def lunar_forcing(steps, dt, lat_deg=10.0, noon_par=1100.0, t_mean=29.5, t_swing=3.0, k=1.5e-5, wind=2.0,
                  period_day=29.53059):
    """A plain lunar day for trying the model: light as the Sun's cosine above the horizon, the water warmest late in
    the day and coolest late in the night, a steady wind and transfer velocity."""
    t = (np.arange(steps) + 0.5) * dt
    phase = 2.0 * math.pi * t / period_day
    par0 = noon_par * math.cos(math.radians(lat_deg)) * np.maximum(np.cos(phase), 0.0)
    temp = t_mean + t_swing * np.cos(phase - 0.25 * math.pi)
    full = lambda v: np.full(steps, float(v))
    return dict(t=temp, t_open=full(t_mean), par0=par0, k=full(k), wind=full(wind))


def main(years=8.0, tp_open=20.0, dt=1.0 / 12.0):
    """The example community seeded in sequence over three years, in four boxes (1, 2 and 4 m trading three tenths of
    their water a day, and 2 m trading a twentieth), without and with a meadow that idles in the dark: each pool near
    the end of every year. A demonstration that the pieces interact; its numbers belong to the placeholder traits."""
    steps = int(years * 365.25 / dt)
    lunar_day = int(round(29.53059 / dt))
    setting = Setting(depth=np.array([1.0, 2.0, 4.0, 2.0]), exchange=np.array([0.3, 0.3, 0.3, 0.05]), tp_open=tp_open)
    schedule = [(0.0, 'phyto', 0.05), (0.0, 'peri', 0.2), (0.0, 'chara', 1.0), (0.0, 'zoo', 0.01), (0.5, 'snail', 0.05),
                (1.5, 'cray_a', 0.05), (3.0, 'fish_a', 0.02)]
    pools = ('chara', 'peri', 'phyto', 'det_c', 'zoo', 'snail', 'cray_j', 'cray_a', 'fish_j', 'fish_a')
    for idle in (1.0, 0.25):
        model = Littoral(example_community(idle_share=idle), setting, 4)
        state = model.initial()
        out = drive(model, state, lunar_forcing(lunar_day, dt), dt, steps, schedule, record_every=lunar_day)
        year = out['time_day'] / 365.25
        print(f'meadow upkeep in the dark: {idle:g} of its upkeep in the light; g C per m2 in each box')
        for y in np.arange(0.9, years, 1.0):
            i = int(np.argmin(np.abs(year - y)))
            print(f'  year {year[i]:.1f}  ' + '  '.join(
                f"{k} {np.array2string(out['series'][k][i], precision=2, suppress_small=True)}" for k in pools))
        last = year > years - 1.0
        print('  lowest oxygen in the last year, g/m3:', out['o2_low'][last].min(0).round(2))


if __name__ == '__main__':
    main()
