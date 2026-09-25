# Middle atmosphere: methods

This model finds the temperature and ozone of a column above a fixed surface
temperature when the temperature, the photochemistry and the sunlight that
drives both are mutually consistent. It is one-dimensional and represents the
global and diurnal mean. The Earth control runs through the same code as the
Moon.

## Column

Levels are equally spaced in log pressure from the surface to 0.1 Pa (100
layers). Heights follow hydrostatic balance with gravity GM/r², as in
`atmosphere/radiative_convective/thermodynamics.py`. On the Moon 0.1 Pa lies
about 600 km up, where spherical shells matter. Dry air is Earth-like, with
Earth's O2 partial pressure: 17.5% O2 at 1.2 atm, 21% at 1.0 atm. CO2 is 400 ppm.

The troposphere follows the saturated pseudoadiabat (Ding and Pierrehumbert
2016) from the surface, with the Manabe–Wetherald relative-humidity profile.
Above the tropopause, water vapour keeps its tropopause mixing ratio. The
exception is the Earth control, which prescribes the observed stratospheric
5 ppm; one lunar sensitivity case does the same. Every level above the
tropopause is in radiative equilibrium. The tropopause is the lowest level at
which that equilibrium is stable against the moist adiabat and the top
convective layer still cools radiatively. The solver moves it one level at a
time until both conditions hold.

## Thermal radiation

Thermal radiation uses the correlated-k scheme in
`atmosphere/radiative_convective/ck.py`: 16 bands from 1 to 3000 cm⁻¹, with 16
g-points each (8 Gauss points on [0, 0.9] and 8 on [0.9, 1]). The
k-distributions come from the repository's line-by-line spectroscopy (HITRAN
lines with Voigt profiles, the MT_CKD 4.3 water continuum, HITRAN
collision-induced absorption). They are computed on a node grid of 27
pressures (0.1 Pa to 3.2×10⁵ Pa, four per decade), 12 temperatures (100–375 K)
and five water fractions (0–0.5). Each node spectrum is sampled at a quarter of
the local Voigt half-width, between 0.0002 and 0.01 cm⁻¹, so the Doppler cores
that dominate cooling at low pressure are represented. Resolution tests change
k by less than 2% at any g-point. Gases combine by random overlap with
resorting and rebinning (Amundsen et al. 2017). Fluxes use the linear-in-τ
source function of `longwave.py` with four Gauss angles. Accuracy against the
line-by-line model is recorded in `atmosphere/radiative_convective/validation.json`.

## Sunlight

- **202–500 nm** (`photolysis.py`): 1-nm bins of the TSIS-1 HSRS solar
  spectrum, multiplied by the shield's transmission. Absorbers are O2
  (Herzberg), O3 (temperature-dependent between 218 and 295 K), NO2, NO3, N2O5,
  N2O, HNO3, H2O2 and the O2–O2/O2–N2 collision-induced bands. Rayleigh
  scattering follows Bodhaine et al. (1999). The radiative transfer is the
  delta-two-stream adding solver of `shortwave.py`, with its pseudo-spherical
  direct beam and six Gauss nodes in the solar zenith cosine over the sunlit
  hemisphere.
- **500 nm to 5 µm**: the line-by-line solar calculation of
  `radiative_convective` (H2O, CO2 and O2 lines, continuum, collision-induced
  absorption, Rayleigh scattering), with the ozone Chappuis and Wulf absorption
  added. It runs on every second level. The heating per unit mass is then
  interpolated in log pressure onto the fine layers, keeping the column total.
  Sharing each coarse layer's absorption evenly would give a stepwise heating
  profile, which the thin upper layers turn into alternating level
  temperatures.
- All absorbed energy heats locally. Sunlight shortward of 202 nm is not
  included; every shield considered blocks it.

## Photochemistry

Fourteen species are solved: O, O(¹D), O3, H, OH, HO2, H2O2, H2, NO, NO2, NO3,
N2O5, HNO3 and N2O. There are 42 reactions and 10 photolysis channels, with
NASA/JPL 19-5 rate coefficients (Burkholder et al. 2019). O + O + M comes from
Campbell and Gray (1973). Cross-sections and quantum yields are the JPL
recommendations from the MPI-Mainz UV/VIS Spectral Atlas (`inputs.json`).
Photolysis rates use the actinic flux: the direct beam averaged over each
layer's exponential attenuation, plus twice the diffuse hemispheric fluxes.

Transport is eddy diffusion of mixing ratio through spherical shells. The
Earth profile is 10⁵ cm² s⁻¹ in the troposphere and 3×10³ cm² s⁻¹ at the
tropopause, growing as p^-1/2 above it (Massie and Hunten 1981) up to a cap of
10⁷. For the Moon the whole profile is multiplied by 36, the square of the
gravity ratio. That keeps the Earth mixing time per scale height, since the
lunar scale height is six times larger. One sensitivity case keeps the Earth
values. How fast the Moon's middle atmosphere actually mixes depends on its
circulation, which only a 3-D model can supply.

Lower boundary conditions:

- fixed surface mixing ratios for N2O (330 ppb, or 0) and H2 (0.53 ppm);
- deposition velocities: O3 0.1, HNO3 1, H2O2 0.5, N2O5 1 and NO2 0.1 cm s⁻¹;
- washout of HNO3, H2O2 and N2O5 in the troposphere, at precipitation divided
  by column water, with precipitation 1 m yr⁻¹.

No species crosses the top. The steady state is reached by backward-Euler steps
of growing length, each solved by Newton iteration with an analytic Jacobian
(tested against finite differences). Photolysis rates are refreshed from the
evolving composition.

## Coupling and outputs

Each outer iteration solves the chemistry for the current temperatures. It then
recomputes the sunlight absorbed in each layer and solves radiative equilibrium
above the tropopause by Newton's method, holding the ozone mixing ratio and the
absorbed sunlight fixed. The Jacobian comes from finite differences of the
correlated-k heating. Iteration stops when temperatures change by less than
0.3 K and the ozone profile by less than 1% of its peak (or 10 ppb).

Energy balance is posed on layers while temperatures live on levels. An
alternating level pattern therefore leaves every layer mean, flux and rate
almost unchanged. It is the Jacobian's weakest mode, localised just above the
tropopause, where radiative relaxation takes 100–400 days. Each Newton step
therefore adds a weak penalty (10⁻³ K/day per K) on the second difference of
the level temperatures, excluding the physical kink at the tropopause. The
smooth profile this selects leaves a heating residual of up to 0.01 K/day,
worth about 1–2 K in the lowest stratosphere. Reported temperatures are layer
means, and the tropopause test compares layer means with the adiabat.

The top-of-atmosphere imbalance at the fixed surface temperature is reported.
It is the forcing the clear column lacks to hold that temperature.

A linear, radiation-only estimate of the day–night swing above the tropopause
comes from the equatorial local heating over the solar day: 29.5 days on the
Moon, 1 day for the Earth control. Its first harmonic drives the
radiative-equilibrium Jacobian, with exchange between levels included; no
dynamics are involved. For the whole column, the mean outgoing longwave over
one night divided by the column's heat capacity gives the cooling that
radiation alone would cause if no heat were carried in.

Surface UV indices use the CIE erythemal action spectrum. They are given for
the sub-solar point, a 45° solar elevation and the global mean.

Balance temperatures (`balance.py`) come from the same columns solved at
surface temperatures of 268, 278, 288, 298 and 308 K (278–298 K for the Earth
control). Each state has its own ozone, tropopause and stratosphere. The
top-of-atmosphere surplus is interpolated linearly between solved states to
where it cancels an assumed cloud effect.

## Upper air out of local thermodynamic equilibrium

Above about 5 Pa, collisions no longer keep CO2's bending mode populated at
the local temperature. The `_nonlte` cases therefore use a two-level source
function in the 500–820 cm⁻¹ bands above 50 Pa (`ck.CKLongwave.fluxes_nonlte`):
S = (J + εB)/(1 + ε). Here J is the absorption-weighted mean intensity of the
layer, and ε is collisional deactivation of CO2(01101) by N2, O2 and O divided
by the Einstein A of 1.5 s⁻¹. The rates are those compiled by López-Puertas
and Taylor (2001); the one for atomic oxygen is uncertain by about a factor of
two.

The band cores stay optically thick up to the model top, so S comes from one
linear solve per band rather than iteration. The change non-LTE makes in
isothermal-layer transfer is added to the LTE fluxes, so where ε is large the
result is exactly LTE. Near-infrared sunlight absorbed above 50 Pa is taken
two ways:
- as heat in full;
- in the collisionally deactivated fraction ε/(1 + ε) only (`_nir_collisional`).

The two bound how much of it heats the air. On the Earth control the
collisional bound stays within 6 K of the US Standard Atmosphere (1976) from 10
to 0.3 Pa, while heat in full is 40 K too warm at 0.3 Pa.

## Exobase

`escape.py` gives the thermal column its base temperature at 0.3 Pa, where the
correlated-k and line-by-line temperatures agree within about 1 K. It also
repeats the calculation from the top layer (0.1 Pa) and from 10 Pa as
sensitivities.

Heat deposited above the base is swept. It is also estimated for a filter that
passes 0.1% of sunlight below its edge, from the WHI 2008 reference spectrum
(Woods et al. 2009, near solar minimum; about 2.5 times more near solar
maximum):
- everything below 121 nm and the Schumann–Runge continuum (122.5–175 nm) is
  absorbed above the base;
- Lyman-α is absorbed there in the fraction the O2 column allows;
- a heating efficiency of 0.4 converts absorbed energy to heat.

Atomic oxygen, which the thermal column does not hold, is carried on each
solution as a trace gas starting from the chemistry's fraction at the base.
Above the base its fraction f follows the zero-flux balance of eddy diffusion K
and its molecular diffusion D through the air:

d ln f/dr = D/(D + K) · (m_air − m_O) g/(kT),

with D = b/n and b = 9.69×10¹⁶ T^0.774 cm⁻¹ s⁻¹ for O through N2 (Banks and
Kockarts 1973). Three cases are reported:
- well mixed (K ≫ D to the exobase), the lower bound;
- separated from the base up (K = 0), the upper bound;
- K from the chemistry's eddy profile (Massie and Hunten, times the lunar
  scaling) continued above the model top, so that separation sets in at the
  homopause, where D = K. This is the central estimate, and it rests on that
  continuation.

Jeans escape of the oxygen at the column's exobase is reported for each. An
oxygen-rich upper thermosphere would also cool by 63-µm emission, which the
column lacks, so the estimates overstate the exobase temperature and the loss
where oxygen is abundant. Where the separated fraction at the exobase is large
(tens of percent), oxygen would change the column itself, and the trace-gas
treatment no longer holds.

## Radiation benchmarks

`benchmarks.py` recomputes, on four final profiles:
- thermal fluxes, line by line at 0.01 cm⁻¹ (H2O, CO2 and O3 lines, MT_CKD 4.3,
  collision-induced absorption), beside the correlated-k fluxes;
- solar direct, diffuse-down and upward fluxes at every level, for zenith
  cosines 1.0 and 0.5 and the global mean.

A 3-D model's radiation code can be checked against them on identical columns.

## Limits

- One-dimensional global and diurnal mean. Clouds, circulation, waves and
  tides are absent. The month-long lunar day drives strong day–night
  circulations that this model cannot represent. The linear swing estimate is
  what radiation alone would do.
- The standard cases treat CO2 emission in local thermodynamic equilibrium,
  which fails above roughly 1–10 Pa; the `_nonlte` cases replace it with the
  two-level approximation above. The correlated-k CO2 cooling adds its own
  error there: 20–40% above about 30 Pa (Earth) or 5 Pa (Moon), and several
  times in the top layer. `lbl_check.py` measures the resulting temperature
  change on finished cases.
- No chlorine, bromine, methane, CO or hydrocarbon chemistry, HO2NO2, aerosol
  surface chemistry or ions. Tropospheric ozone production is therefore absent.
- No sunlight below 202 nm, so NO, H2O and Schumann–Runge O2 photolysis are
  missing. On the Moon behind any of these shields this is the physical
  situation. For the Earth control it is not: its NOy has no upper sink, which
  is the main reason the control underestimates Earth's ozone column.
- The Earth-based eddy diffusion, its lunar scaling, the washout rate and the
  stratospheric humidity rule are assumptions that the sensitivity cases vary.

## References

- Amundsen D. S. et al. (2017) Astron. Astrophys. 598, A97 (random overlap with resorting and rebinning).
- Banks P. M. & Kockarts G. (1973) Aeronomy, Part B, Academic Press (molecular diffusion coefficients).
- Bodhaine B. A. et al. (1999) J. Atmos. Oceanic Technol. 16, 1854 (Rayleigh scattering).
- Burkholder J. B. et al. (2019) Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies, Evaluation No. 19, JPL Publication 19-5.
- Campbell I. M. & Gray C. N. (1973) Chem. Phys. Lett. 18, 607 (O + O + M).
- Coddington O. M. et al. (2021) Geophys. Res. Lett. 48, e2020GL091709 (TSIS-1 HSRS).
- Ding F. & Pierrehumbert R. T. (2016) Astrophys. J. 822, 24 (non-dilute moist adiabat).
- Johnston H. S. et al. (1996) J. Phys. Chem. 100, 4713 (NO3 photolysis yields).
- Keller-Rudek H. et al. (2013) Earth Syst. Sci. Data 5, 365 (MPI-Mainz UV/VIS Spectral Atlas).
- López-Puertas M. & Taylor F. W. (2001) Non-LTE Radiative Transfer in the Atmosphere, World Scientific.
- Manabe S. & Wetherald R. T. (1967) J. Atmos. Sci. 24, 241.
- Massie S. T. & Hunten D. M. (1981) J. Geophys. Res. 86, 9859 (eddy diffusion).
- Matsumi Y. et al. (2002) J. Geophys. Res. 107, 4024 (O(1D) quantum yields).
- McKinlay A. F. & Diffey B. L. (1987) CIE Journal 6, 17 (erythemal action spectrum).
- Sander S. P. et al. (2011) JPL Publication 10-6 (recommended cross-sections).
- Troe J. (2000) Z. Phys. Chem. 214, 573 (NO2 photolysis yields).
- U.S. Standard Atmosphere, 1976. NOAA, NASA and USAF, Washington, D.C.
- Woods T. N. et al. (2009) Geophys. Res. Lett. 36, L01101 (WHI 2008 reference solar spectrum).

These references record where equations and data come from. Full-source reading
for manuscript use remains outstanding under the ensemble charter.
