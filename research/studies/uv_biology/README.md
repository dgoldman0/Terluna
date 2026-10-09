# Biology under the shield and a restricted ultraviolet option

The working choice is to engineer the biosphere for the existing filtered spectrum and use local ultraviolet
where it serves a defined function. A narrow, attenuated UV-B option remains worth investigating. The EUV and
short-wave rejection requirements remain controlling; no new transmission band is adopted here.

The [runner](run.py), `python -m research.studies.uv_biology.run`, reweights the committed
[surface-light product](../../../illumination/surface_light/results/surface_light.json). Its paired filtered and
unfiltered solar cases use the same almost ozone-free lunar column, making monochromatic transfer linear in
incoming intensity. Direct and diffuse light, atmospheric multiple scattering and ground albedo 0.1 are retained.
The [results](results/uv_biology.json) pin that product, the code and the shared constants.
The accompanying [primary-source register](../aerial_ecology/sources.json) records the photobiological evidence
and access limits. This is an independently checked transfer calculation, rather than empirical validation of
the prescribed UV windows.

## The initial numerical comparison

The proposed windows replace transmission only inside the named band. All other wavelengths retain the reference
filter. These are prescribed spectra, without a manufacturable film design. Percentages are absolute incoming
transmission. The 310–400 nm case opens a wider comparison band; it is not the old 310-nm equilibrium experiment.

| Prescribed spectrum | Overhead surface UV index | UV index with Sun at 30° | Added overhead ground UV irradiance |
|---|---:|---:|---:|
| Reference film and dimmer | 0.109 | 0.039 | 0 |
| 295–315 nm at 1% | 0.371 | 0.133 | 0.0208 W/m² |
| 295–315 nm at 5% | 1.418 | 0.507 | 0.1042 W/m² |
| 295–315 nm at 10% | 2.727 | 0.975 | 0.2084 W/m² |
| 310–315 nm at 10% | 0.235 | 0.084 | 0.0698 W/m² |
| 310–400 nm at 95% | 2.261 | 0.809 | 14.976 W/m² |

The values are conditional instantaneous surface exposures. The energy column is downward ground irradiance,
not a global absorbed forcing or upper-atmospheric heat budget. The UV index is erythemally weighted; it does not
measure vitamin-D yield, plant signaling, DNA damage for each organism or the accumulated dose over a long lunar
day. The product's unweighted photon bands are clearly labelled proxies. Low-Sun two-stream limitations in the
parent product carry through.

These results establish that a restricted band can produce a substantial surface dose even in the deeply
scattering atmosphere. **Strong scattering alone does not require global UV-B transmission.** Biological benefit,
damage and atmospheric consequences remain the decision criteria. Narrower 310–315 nm transmission supplies none
of the 295–305 nm proxy band; its benefit cannot be inferred from an equal UV index.

## Corrections to the discussion

The present film suppresses much of **315–350 nm short UV-A** as well as UV-B. The stored spectrum puts the
overhead short-UV-A photon flux at about 0.194% of the Earth control. Almost all retained UV-A lies farther toward
350–400 nm. The claim that the current design already preserves short UV-A was inaccurate. Rai et al. (2020,
[doi:10.1111/pce.13752](https://doi.org/10.1111/pce.13752)) demonstrate UVR8 responses extending toward 350 nm;
receptor adaptation and the available wavelength distribution need evaluation together.

The historical 310-nm cutoff gave a surface UV index near 1.1 and 175 K at 1 Pa in an LTE column. It assumed 0.1%
leakage below the edge and developed about 138 DU of ozone. Neither its ozone nor its temperature belongs to the
new clean-window cases. In particular, 1 Pa is not the exobase. The current shield's full solar-cycle and aperture
feedback requirements must still be tested.

Preserving rejection below about 240 nm maintains the principal O2-photolysis restriction. Admitted UV-B can
nevertheless photolyse existing ozone and alter OH chemistry; biological NOx/VOC sources also matter. Recent
oxygen-transport calculations already produce a small ozone column under the shield. Holding an almost ozone-free
column fixed is therefore a first optical comparison, not a chemically self-consistent atmosphere. Additional
UV-B cannot be assumed either to restore terrestrial atmospheric cleansing or to leave trace gases unchanged.

## What is locked in

1. Preserve the current spectral shield as the reference, including its EUV/far-UV rejection and cold-exobase
   acceptance criteria. Model wavelengths below the product's 202-nm boundary as explicitly unchanged; the
   numerical equality check covers the available 202–295 nm bins only.
2. Develop biological provision of essential compounds, photoreceptors and ecological signals under the available
   spectrum. Local controlled UV remains an option for processes that benefit from it. Dietary vitamin D is a
   supply route; the route that produces it without a global natural UV-B field still needs to be specified.
3. Retain a narrow attenuated UV-B window as a candidate to test. The new surface screen justifies that test and
   gives exposure scales; no transmission percentage is selected.

Before adopting a window, use the current ring geometry, leakage and solar-cycle spectra in a coupled
radiative/photochemical calculation with biological emissions; calculate actinic doses through the aerial biomes;
couple the non-LTE middle air to the current cooled exobase and shadow-expansion feedback; and include film
absorption, infrared re-emission and compensating climate dimming. Vitamin-D and taxon-specific action spectra
must replace band proxies for those decisions. Wet biological aerosols introduce additional opacity and chemistry.

The full equilibrium inputs and correlated-k caches are absent from this checkout. Rebuilding them would require
external data restoration and hours of table generation, and the present chemistry still lacks the complete
living-atmosphere reactions. This bounded calculation uses the stored transfer product and leaves that
larger evaluation explicitly open. No new climate run or relaxed protection requirement is implied.


Executed numerical and repository checks, including unsuccessful full-suite gates and their causes, are in the
[check record](../aerial_ecology/checks.json). The new aerial and UV checks pass; the whole repository is not reported green.
