# Radiative–convective column: methods

A one-dimensional, line-by-line calculation of the clear-sky energy budget of a
prescribed atmospheric column, for the Open Moon and for an Earth control run
through the same code. It answers one question: for a given surface temperature,
how much thermal radiation leaves the top of the atmosphere (OLR) and how much
sunlight the planet absorbs (ASR)? Comparing the two curves shows where a column
could balance, and the largest OLR a water-rich column can emit is its
runaway-greenhouse limit. This is the "inverse" method of the habitable-zone
literature (Kasting 1988; Kopparapu et al. 2013; Goldblatt et al. 2013).

## What is prescribed and what is computed

**Prescribed (inputs):** gravity and radius of the planet; the dry surface
pressure and composition (N2, O2, Ar, CO2); the surface temperature; the
relative-humidity profile; the stratospheric temperature; the surface albedo.
Clouds, ozone, methane, nitrous oxide and aerosols are absent.

**Computed:** the column's temperature and water profile, its heights and
layer masses, spectral optical depths, and the thermal and solar fluxes.

The temperature structure is not solved from radiative equilibrium. Below the
tropopause it follows a saturated pseudoadiabat; above, it is isothermal at the
prescribed stratospheric temperature. That is the standard structure for
habitable-zone limits and a known simplification. Where a result depends on the
stratospheric temperature or the humidity profile, the sweeps vary them.

## Column

- **Saturation pressure:** Wagner & Pruß (2002) over liquid water above
  273.16 K and Wagner et al. (2011) over ice below. The latent heat is the
  Clausius–Clapeyron slope of that curve for an ideal vapour. It is 1.6% above
  the tabulated value at 373 K; the adiabat and the vapour curve stay consistent
  with each other.
- **Moist adiabat:** Ding & Pierrehumbert (2016), eq. 12. It is non-dilute, so
  it holds from vapour-poor air to a steam-dominated column, and is integrated in
  ln p by fixed-step RK4. Tests check the dilute textbook limit and the
  pure-vapour Clausius–Clapeyron limit.
- **Humidity:** Manabe & Wetherald (1967) relative humidity,
  RH = 0.77 (p/p_s − 0.02)/0.98, or saturation (RH = 1) for runaway limits.
  The stratosphere keeps the tropopause water mole fraction.
- **Hydrostatics:** geopotential from the virtual temperature, converted to
  height with g = GM/r². Each layer's mass is Δp/g(z). The Moon's tall column
  therefore holds about 5% more air per unit area than p_s/g_s.
- **Composition:** oxygen is held at Earth's partial pressure of 21.2 kPa. At
  1.2 atm that is 17.5% O2, the repository's design case; at 1 atm it is 20.95%.
  Argon keeps Earth's Ar/N2 ratio, and CO2 is a stated dry mole fraction.

## Spectroscopy

Inputs are fetched and hash-checked by `fetch_inputs.py` (see `inputs.json`) and
are not committed:

- **Lines:** HITRAN as served by hitran.org on 2026-09-24, H2O (7
  isotopologues), CO2 (12) and O2 (3), with HITRAN partition sums.
  - Intensity floor: 10⁻²⁸ cm/molecule at 296 K for H2O and CO2; all O2 lines.
  - Line strength: scaled in temperature with the partition-sum ratio, lower-state
    energy and stimulated emission.
  - Widths: air and self broadening with n_air, plus the air pressure shift.
  - Line shape: Voigt, cut 25 cm⁻¹ from centre.
- **Water continuum:** MT_CKD 4.3 (Mlawer et al. 2012, 2023; © AER). It is the
  AER formulation ported to Python: density scaling, self temperature
  exponents, the radiation term and AER's four-point interpolation. It
  reproduces AER's own example output to 5×10⁻⁷. Consistent with MT_CKD, water
  lines have their value at 25 cm⁻¹ (the pedestal) subtracted.
- **Collision-induced absorption:** HITRAN CIA for N2–N2, O2–N2, O2–O2 and
  N2–H2O, interpolated linearly in temperature and clamped outside the tabulated
  range. Where two datasets from different sources cover the same band at the
  same temperatures (O2–O2 near 20,000–29,800 cm⁻¹), only the first is used, so
  nothing is counted twice.
- **Rayleigh scattering:** Bodhaine et al. (1999) for air. Test value: τ =
  0.0973 at 550 nm for 1013.25 hPa. Water vapour is taken to scatter 0.75× as
  much per molecule. The composition change from a few per cent less O2 is
  ignored; it is under 1%.

### Line-by-line evaluation

The longwave grid is 1–3000 cm⁻¹ at 0.01 cm⁻¹. The solar grid is 2000–20,000
cm⁻¹ at 0.05 cm⁻¹, plus 20,000–49,500 cm⁻¹ at 2 cm⁻¹, where only Rayleigh
scattering and O2 collision absorption act. Each line is summed in three parts
that add up to its cut-off profile:

1. **Core:** the exact profile on the fine mesh within about 1 cm⁻¹ of centre.
   The Voigt function is used within 50 Gaussian widths and the Lorentzian
   beyond, where they agree to 0.1%.
2. **Near wing:** the Lorentzian on a mesh 10× coarser (4× for the solar grid)
   out to 6 cm⁻¹.
3. **Far wing:** from 4 to 25 cm⁻¹, blended in smoothly. It is computed for all
   lines at once by FFT convolution, using the expansion
   (γ/π)(1/d² − γ²/d⁴), whose neglected term is below 3×10⁻⁴.

Against a brute-force Voigt sum on the fine mesh, the method is within 0.6% at
every grid point and within 10⁻⁴ in integrated strength, for H2O and CO2 from
3 Pa to 1 bar.

CO2 and O2 broaden with air only, so their cross-sections depend only on
pressure and temperature. They are computed once at nodes p = 10^(i/10) Pa and
T = 10j K, cached, and interpolated log-linearly. Against direct calculation
the interpolation error is under 1% at any wavenumber and 0.3% in band
integrals. For an Earth column, OLR moves by 0.01 W/m².

## Radiative transfer

- **Thermal:** layers are non-scattering and plane-parallel, with a source
  function linear in optical depth; four Gauss angles per hemisphere. Tests:
  - The Planck integral recovers σT⁴ to within 10⁻⁴.
  - An isothermal column is in exact equilibrium.
  - A grey slab reproduces the exponential-integral solution (2E₃) to 10⁻⁵
    with 12 angles.
  - The linear source is exactly consistent under subdivision.
- **Solar layers:** the two-stream "practical improved flux method" (Zdunkowski
  et al. 1980), delta-scaled. It gives non-negative reflectance for strongly
  absorbing layers, where the plain Eddington form goes negative. Layers are
  combined by adding over a Lambertian surface. Tests: energy conservation to
  10⁻⁶ for conservative scattering; zero reflection for pure absorbers; one
  layer equals two half-layers to 10⁻¹⁰.
- **Solar beam geometry:** the direct beam follows straight rays through
  spherical shells; each layer's effective cosine reproduces that attenuation,
  so the adding stays energy-conserving.
- **Rayleigh reflectance:** against a Monte Carlo reference for a Rayleigh slab
  over a black surface (τ = 0.1–1.5, μ₀ = 0.25–1), the two-stream reflectance is
  low by 0.1–3%. At the Moon's Rayleigh optical depths that is roughly 0.003 in
  planetary albedo, about 1 W/m².
- **Planetary means:** six Gauss nodes in μ₀ over the sunlit disk.
- **Spectrum and grid limits:** the solar spectrum is TSIS-1 HSRS (202–2730 nm)
  joined to a 5772 K blackbody beyond, scaled to 1361 W/m². Sunlight longward
  of 5 μm, about 1.9 W/m² in the global mean, is counted as absorbed.

## Validation record

[validation.json](validation.json) holds the numbers.

- **Independent code:** PyRADS (Koll & Cronin 2018), run on its own saturated
  profiles at Earth and lunar gravity from 250 to 360 K. Given the same data
  (MT_CKD 3.2 and the main water isotopologue), this code reproduces its OLR
  within 0.4%.
- **Data spread:** with the default data (MT_CKD 4.3, all isotopologues), OLR is
  up to 7 W/m² higher in hot saturated columns, almost all from the continuum
  version. That spread is a data uncertainty in the runaway limit, not a code
  error.
- **Convergence:** spectral step, angles, line floor, solar zenith nodes and
  solar levels each move the budget by less than 0.03 W/m². The thermal vertical
  grid converges at first order and the standard grid reads OLR about 0.6 W/m²
  low.
- **Rayleigh reflectance:** two-stream against Monte Carlo, as above.

## Limits of interpretation

- Clear sky only. Clouds dominate Earth's albedo and would do so on the Moon;
  the results give the clear-sky budget, and any cloud effect is a separate,
  stated assumption.
- The column is one global mean. The Moon's 29.5-day day, its slow rotation, the
  night side and the poles are outside this calculation.
- The temperature structure is prescribed (adiabat plus isothermal
  stratosphere), and the humidity follows an Earth-derived profile. Neither is
  predicted for the Moon.
- There is no ozone. With the spectral shield blocking O2-photolysing
  ultraviolet, an ozone layer is itself uncertain.
- The geometric effect of the Moon's tall atmosphere on how much sunlight it
  intercepts and emits is represented only in the solar beam. Emission from
  large heights (up to (r/R)² ≈ 1.1 more area) is not.
- Numerical convergence and agreement with an independent code are not
  empirical validation.

## References

- Bodhaine B. A. et al. (1999) J. Atmos. Oceanic Technol. 16, 1854.
- Coddington O. M. et al. (2021) Geophys. Res. Lett. 48, e2020GL091709 (TSIS-1 HSRS).
- Ding F. & Pierrehumbert R. T. (2016) Astrophys. J. 822, 24.
- Gordon I. E. et al. (2022) J. Quant. Spectrosc. Radiat. Transfer 277, 107949 (HITRAN2020).
- Goldblatt C. et al. (2013) Nature Geosci. 6, 661.
- Joseph J. H., Wiscombe W. J. & Weinman J. A. (1976) J. Atmos. Sci. 33, 2452.
- Karman T. et al. (2019) Icarus 328, 160 (HITRAN CIA).
- Kasting J. F. (1988) Icarus 74, 472.
- Koll D. D. B. & Cronin T. W. (2018) PNAS 115, 10293.
- Kopparapu R. K. et al. (2013) Astrophys. J. 765, 131.
- Manabe S. & Wetherald R. T. (1967) J. Atmos. Sci. 24, 241.
- Mlawer E. J. et al. (2012) Phil. Trans. R. Soc. A 370, 2520 (MT_CKD).
- Wagner W. & Pruß A. (2002) J. Phys. Chem. Ref. Data 31, 387.
- Wagner W. et al. (2011) J. Phys. Chem. Ref. Data 40, 043103.
- Zdunkowski W. G., Welch R. M. & Korb G. (1980) Beitr. Phys. Atmos. 53, 147.

These references record where equations and data come from. Full-source reading
for manuscript use remains outstanding under the ensemble charter.
