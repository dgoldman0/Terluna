# Böhm hydrodynamics: implementation comparison

The new adapter reproduces a restricted part of McSnow's Böhm implementation and exposes a problem with using it as a simple reduction of the current snow collision kernel. Its drag calculation agrees with the authors' archived Earth graupel examples under the stated surface-air assumption. Its collision calculation agrees with independently compiled author code, but often returns cross-section ratios greater than one for the prescribed aggregate geometries. Full-column adoption is withheld while that behavior and the geometry assumptions are investigated. Existing snow producers and products remain reproducible.

## Calculation and interface

`boehm.py` accepts prescribed mass, outer diameter and effective projected area for a rough particle with a spherical envelope. It evaluates the terminal drag, the two-body boundary-layer approximation and the potential-flow correction described by [Welss, Siewert and Seifert (2024)](https://doi.org/10.1029/2023MS003805). Its numerical reference is [McSnow revision 1e428d7e901cbc5b91f2040e575f736ec5963a74](https://gitlab.dkrz.de/mcsnow/mcsnow/-/tree/1e428d7e901cbc5b91f2040e575f736ec5963a74), whose original equations and code discrepancies are retained in the source record.

```python
properties = particle_properties(D, m, g, state, q, alpha=1, sigma=0)
Ec = collision_efficiency(left, right, state, diagnostics=False)
```

Lengths are metres, mass kilograms and gravity m/s². `state['rho']` is air density in kg/m³ and `state['eta']` dynamic viscosity in Pa s, matching the winter-column interface. Inputs broadcast with NumPy; every returned particle property is an array, including zero-dimensional arrays for scalar inputs. To calculate a pair matrix, give each left field shape `(n,1)` and each right field shape `(1,n)`. Air density, viscosity and gravity must agree between the two particles. Diagnostics include the raw efficiency, its boundary-layer and potential-flow factors, reduced Reynolds number, F/H, b cy and a flag for values outside [0,1]. Finite values are returned unchanged; unsupported geometry and nonfinite calculations raise errors.

Here alpha is the **fall-axis aspect ratio**, fixed at one. It is distinct from a horizontal image aspect ratio. The supplied q is the effective flow area divided by the outer circle area: `Aeff=q*pi*D²/4`, with `0<q<=1`. The source uses `r=D/2`, `r_e=r*sqrt(q)` and `alpha_e=1/sqrt(q)` inside its collision calculation; its reduced two-body aspect ratio generally differs from one. These definitions describe a prescribed porous envelope. They do not reconstruct a real aggregate's orientation, structure or contact history.

The outer gravitational kernel remains explicit in the caller:

```
K = pi*(D1+D2)**2/4 * abs(v1-v2) * Ec
```

The geometric collision radius is D/2 even though effective projected area enters the drag and collision-state transformations. Sticking and breakup are separate operations. Multiplying Ec into the earlier effective-area kernel would change its normalization without documenting the change. The adapter supplies no ambient turbulence or mechanical failure model.

SIG=0 specifies deterministic internal speeds; SIG=.25 reproduces the source's internal rough-ice speed-width setting. In the source, that width changes the expected relative speed inside Ec while the usual outer kernel still uses `abs(v1-v2)`. It is a model parameter, not a measured lunar turbulence intensity. Both controls are retained explicitly.

## Verification completed

The independent reference compiles exact retained source slices with GNU Fortran, using a small driver that supplies particle and air states and sets SIG explicitly. The Python adapter is never imported when generating the Fortran expectations. The original license, copied source slices, generator and their hashes are retained under `references/`; the 122 pair cases are in `results/boehm_reference.json`.

The cases include Earth and lunar gravity, q=.15/.35/1, mixed q and SIG, comparable and disparate sizes, and 24 cases using the local lunar column and its corrected mass/area closures. Both F/H branches and all three potential-flow regimes are exercised. The Python result differs from the scalar Fortran efficiency by at most 3.36×10⁻¹⁰ relatively; the largest absolute difference is 4.10×10⁻⁶ where Ec is approximately 37,734. Drag speed and Reynolds number agree to approximately 4×10⁻¹⁴ relatively in the initial synthetic cases. Larger relative cancellation errors in the inertial drag remainder are covered explicitly by the test tolerance. These comparisons establish implementation agreement.

Five archived Earth graupel ground cases from the authors' [Zenodo dataset](https://doi.org/10.5281/zenodo.10077724) provide an additional check. They use alpha=1, density 800 kg/m³ and diameters 20,50,100,200,500 µm. With the source's default saturated surface—273.15 K, 101325 Pa and vapour pressure 610.78 Pa—the adapter agrees with archived speeds within 4.1×10⁻⁷ relatively and Reynolds numbers within 6.6×10⁻⁷. Source gas constants and the original viscosity formula, including its default-precision literals, are used for this check. The original run namelist was not present in the retrieved material; the result is conditional on that explicit surface state. It does not establish exact identity of the full-height humidity profile or archived source revision.

Independent collision-trajectory validation remains unperformed. The retrieved Wang–Ji curves describe nonspherical primary crystals collecting water droplets; the archived McSnow drop-pair cases also require slightly oblate liquid particles and smooth-surface parameters. Those cases lie outside this adapter's alpha=1 rough-ice scope. None of these checks validates lunar aggregate strength, decimetre geometry or a complete snowfall distribution.

The focused test suite contains 16 passing cases, including the parameterized input checks. It verifies source and fixture identities, scalar and broadcast results, pair-matrix symmetry, equality behavior, retention of efficiencies above one and the archived Earth drag comparison. To regenerate the independent fixture without writing production products:

```
python atmosphere/snow/references/build_boehm_reference.py --output research/runs/snow_scale/boehm_reference.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python -m pytest atmosphere/tests/test_snow_boehm.py -q -p no:cacheprovider
```

## What the local lunar comparison reveals

At the current lunar column's −10°C state, the representative dendrite-area closure uses `m=.0121*D^1.9` and q=.26019796807371026. The following values come from the independent Fortran fixtures; all these particles satisfy the simple solid-ice envelope mass bound, although the large-aggregate geometry remains extrapolated.

| Outer diameters | Ec, SIG=0 | Ec, SIG=.25 |
|---|---:|---:|
| 0.2 mm and 15 cm | 0.6478 | 0.6478 |
| 0.2 mm and 60 cm | 0.2937 | 0.2937 |
| 1 mm and 2 mm | 1.3857 | 1.2032 |
| 1 cm and 2 cm | 5.0833 | 1.4674 |
| 15 cm and 60 cm | 16.5108 | 3.0286 |
| 15 cm and 15.001 cm | 37,734.4 | 2.4141 |

This combination suppresses some small-particle capture and enhances some other encounter kernels. Ec is a ratio of collision cross sections; exceeding one alone does not prove a physically impossible encounter rate. It does show that this implementation cannot be treated as the desired universal suppression factor. Its agreement with the author code does not settle whether the porous spherical-envelope mapping and source approximations apply to these aggregate pairs.

For SIG=0, the source's internal relative speed approaches zero as distinct terminal speeds converge. In the large-F/H branch, F and G then scale as `1/abs(dv)`, so Ec grows approximately as `1/abs(dv)`. Multiplying by the outer `abs(dv)` yields a finite nonzero limiting kernel. At exact equality, Fortran `SIGN(1,LOG(1))` returns +1 for both terms, giving an internal speed 2v, while the outer speed difference is zero. Thus the deterministic control has a discontinuity at equality; its large Ec does not mean an infinite kernel. SIG=.25 makes the internal speed difference finite and removes that particular discontinuity, while still producing efficiencies above one.

A diagnostic sweep over all five existing geometry choices, 181 diameters from 20 µm to 10 m, found the same behavior. Its small end is not physically interpretable under the aggregate mass law: the corrected coefficient exceeds the mass of a solid-ice sphere below about 66 µm, and the legacy coefficient below about 97 µm (using 917 kg/m³ as a diagnostic ice density). The 20 µm ratios are approximately 3.72 and 5.68. The millimetre and centimetre examples in the table show that the efficiency behavior persists above that invalid lower range. There is no justified new giant-flake coefficient or upper-size cutoff in this calculation.

## Source choices and remaining limits

The adapter follows the source's F/H threshold of **10**, its exponent on both boundary-layer branches, its complete exponential factor in the potential-flow term, `b=3*v1/(r_e1+delta1)`, and kstar in G. These differ from omissions or ambiguous parentheses in the displayed 2024 equations. Both source and paper retain `(F+G)/(2F)` in the large-F/H approximation; a direct expansion of the other branch includes an additional one in its numerator. That approximation is preserved. The default-precision Fortran literals 9.18 and 3.53 are represented with their float32 rounding, which matters for numerical parity. Stable algebra avoids unnecessary overflow; the b cy=.5 and 1 endpoints use their analytic limits.

The retained source also contains an optional collector-aspect factor with a TODO, which has no effect at alpha=1. Its diagnostic height-kernel routine uses effective-area radii, unlike the production geometric kernel, and its Wang–Ji driver assigns `ty='oblate'` while checking `'oblates'` in two geometry branches. Those source differences prevent assuming that every published curve follows directly from the pinned master revision. Retrieval failures, reading extents and the selected dataset rows are recorded in `boehm_sources.json`.

Lunar gravity and local air density and viscosity enter the continuum equations explicitly. Their use provides a physically specified transfer calculation, with inherited roughness, orientation and flow assumptions. A future comparison of evolving geometry needs at least monomer count and a stated rule for deposition, sticking and fragmentation. This adapter changes none of those closure assumptions, and its raw collision efficiency remains an experimental diagnostic pending a defensible applicability and normalization assessment.
