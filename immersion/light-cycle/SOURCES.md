# Sources and provenance

Inspected through the web and GitHub connector on 20–21 September 2026. The
container had no working external DNS, so no external package is claimed to have
been downloaded, installed or executed. This implementation is a newly written
scalar solver, informed by the methods below. It is not libRadtran/MYSTIC and
not the Bruneton GPU implementation.

1. Eric Bruneton. **Precomputed Atmospheric Scattering: a New Implementation**
   (2017).
   https://ebruneton.github.io/precomputed_atmospheric_scattering/
   Method and implementation reference for spherical transmittance, successive
   scattering orders, surface reflection and spectral display conversion.
   The annotated equations and code in
   https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/functions.glsl.html
   were inspected, particularly the transmittance, Rayleigh phase, scattering
   density, multiple scattering and irradiance sections. The upstream validation
   described there belongs to that implementation, not to this new solver.

2. Bruneton reference demonstration, `atmosphere/demo/demo.cc`.
   https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/demo/demo.cc.html
   https://github.com/ebruneton/precomputed_atmospheric_scattering/blob/master/atmosphere/demo/demo.cc
   Inspected blob: `edd1ae134a8a79c6ff7254682a26cc99351a75ea`.
   Source of the 48 binned extraterrestrial irradiance values, 360–830 nm, and
   the 48 ozone absorption cross sections at 233 K. The upstream comments identify
   the solar data as ASTM G173 and the ozone data as University of Bremen data.
   Those original datasets were not independently retrieved for this deliverable.
   The reference Rayleigh coefficient, effective 8 km scale height and triangular
   300 DU ozone profile are also taken from this transparent reference setup.
   Its BSD-3-Clause notice is reproduced in LICENSES.txt.

3. Chris Wyman, Peter-Pike Sloan and Peter Shirley (2013).
   **Simple Analytic Approximations to the CIE XYZ Color Matching Functions**.
   Journal of Computer Graphics Techniques 2(2), 1–11.
   https://jcgt.org/published/0002/02/01/paper.pdf
   The complete paper was read; the page containing Table 1 was visually checked.
   The asymmetric multi-lobe fits provide the CIE 1931 2-degree x,y,z functions
   used in `render_data.py`. Coefficients are implemented as equations, not as
   a copied image or reproduction of the paper. The shortwave truncation and
   20 nm transport sampling remain separate approximation choices.

4. libRadtran official **Basic usage** documentation.
   https://www.libradtran.org/doku.php?id=basic_usage
   Establishes the relevance of fully spherical twilight calculations, and the
   limitations of plane-parallel and pseudospherical treatments below the horizon.
   MYSTIC examples were inspected as a method lead. None of its published
   benchmark numbers is reported as a test performed on this implementation.

5. Terluna repository, `illumination/README.md` on main, inspected before work.
   https://github.com/dgoldman0/Terluna/blob/main/illumination/README.md
   Blob: `d905884590279c4457cd4c8aa1cc511b18d28930`.
   Establishes that the pre-existing illumination work consisted of an angular
   sweep diagnostic and that scattering, colour and refraction were still open.
   Lunar period 29.53059 days and the intended 1.2-pressure-scale atmosphere are
   inherited project inputs. The prescribed optical model in this deliverable
   is explicitly separate from the repository's upper thermal-column solver.

All source URLs are provenance links, not runtime dependencies of the HTML.

## Full-cycle extension and weather display

- Pharr, Jakob & Humphreys, *Physically Based Rendering*, fourth edition,
  "Transmittance": https://pbr-book.org/4ed/Volume_Scattering/Transmittance
  Used as the primary reference for exponential path transmittance. The weather
  display is an independently implemented reduced approximation, not PBRT code.
- libRadtran, "Basic Usage", especially cloud inputs and solver geometry:
  https://www.libradtran.org/doku.php?id=basic_usage
  This describes the inputs and more complete scattering tools that future
  cloud validation would require. The viewer does not run libRadtran or MYSTIC.
- Original solver and new numeric data were rerun locally for this version.
  `validation/noon_monte_carlo.json` records six independent noon spot checks.

No outside landscape imagery, tree models, font files or third-party JavaScript
libraries are redistributed. The landscape is procedurally constructed. The
third-party solar/ozone-array notice is retained in the source and bundled HTML.
