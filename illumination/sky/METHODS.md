# Clear-sky solver: methods

This scalar spherical solver computes spectral sky radiance and surface
irradiance for Earth, the Open Moon with a stretched 300 DU ozone layer, and the
Open Moon without ozone, for Sun elevations from −90° to +90°. It was written for
the "month of light" viewer (now in
[visualization/month-of-light](../../visualization/month-of-light/)) and moved here
because it is illumination science, not display code. The display models the
viewer adds (clouds, fog, water, exposure) are documented with the viewer.

## Atmospheric calculation

The underlying scalar spherical optical solver is the first viewer's solver,
with its ground scattering-order increment criterion extended from the interval
−24° to +12° through **+90°**. The original algorithm and full profile definitions
are documented in `validation/legacy_v1/METHODS_v1.md`. This document records the
extension.

Earth radius: 6,371 km; molecular density scale height: 8 km; surface density
multiplier: 1. Lunar radius: 1,737.4 km; scale height: 48.428 km; density multiplier:
1.2. The exponential profiles are chosen optical proxies. They are not the
project's solved variable-gravity thermal/chemical atmospheric columns. The
neutral ground albedo is 0.1; the source is the unfiltered reference solar
spectrum. The Moon has either a stretched 300 DU ozone profile or zero ozone.

Fields were recomputed on the standard spatial/angular source grid and 22
wavelengths, 380–800 nm at 20 nm spacing. Earth stopped after 14 scattering orders;
both lunar cases after 55. The last-order increment remains a convergence
indicator, rather than a certified bound on the sum of all omitted orders.
The calculated panorama grid is now **174 Sun elevations from −90° to +90°**, with
81 quadratically spaced view elevations and 65 azimuths over the symmetric half
sky. Above +12°, the output angle spacing is 1° through 30°, then 2° through 90°.
The underlying diffuse-source grid remains coarser away from twilight, especially
below −35°. Raw signed linear RGB and its metadata remain in the NPZ archives.

Cloud illumination inputs are sampled from the same optical source field at a
prescribed altitude of 2.5 km. They provide a colour and intensity reference for
the visual cloud layer; this does not turn the weather renderer into a coupled
cloud/molecular transport solution.

## Precision and open checks

The atlas archives are float32 RGB with the original source metadata. The earlier
study's global lunar energy-accounting residual of up to 4.6% remains open; its
records are in `validation/legacy_v1/`. `test_solver.py` (21 tests) covers the
geometry, transport and quadrature routines, `test_atlas.py` (10 tests) checks
invariants of the generated atlases, and `validation/noon_monte_carlo.json`
records six independent 100,000-path noon Monte Carlo spot checks
(`check_noon.py`). These address selected numerical properties, not a whole-sky
error bound or empirical lunar-atmosphere validation. Deep-night outputs inherit
the coarser diffuse-source grid below −35° and are exploratory.
