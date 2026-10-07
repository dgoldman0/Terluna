# Atmospheric protection

The [primary magnetic architecture comparison](../research/studies/solar_shield_array/magnetic_architecture.md)
adds exact finite-loop fields in `magnetic_fields.py` and compares larger
regional circuits, conditional weaker moments and upstream sources. It finds
a strong regional material advantage while leaving plasma transport and final
hardware selection open. The imported September implementation remains intact.

`model.py`, `verify.py` and the compact result tables are verbatim imports from the September 2026 protection package. `results/results.json` retains selected original constants, optical summaries, reference/storage cases and input hashes for the verifier; duplicate table arrays and unrelated summaries are omitted, with values unchanged. They calculate candidate optics, aperture geometry, lunar-phase holding forces, propulsion/power/storage feedback, magnetic structures, particle-rigidity diagnostics and renewal budgets.

The protected atmosphere and transport comparison are inherited from the [September feasibility baseline](../research/baselines/feasibility/), with fixed inventory constants in the original code. Their use is not an independent atmospheric validation.

## Inputs and reproduction

Five external optical/solar input files are pinned in [inputs.json](inputs.json). They are deliberately outside this source commit; data attribution, retrieval URLs, sizes and exact checksums are retained. The original files are available in the recovered project ZIP in the current session. Restore them from that package or fetch the recorded upstream bytes:

```sh
python protection/fetch_inputs.py --archive /path/to/Lunar_Protection_Model.zip
# Alternatively, when network access is available:
python protection/fetch_inputs.py --download
python protection/verify.py
```

The verifier imports `model.py`, which immediately loads the TiO2 table; restoring inputs is required even for that original verifier. A SHA-256 mismatch stops retrieval. Upstream changes require an explicit reviewed update; hashes are never silently rewritten.

After installing `research/requirements.txt`, run `python protection/model.py` to regenerate the optical optimization and all reference outputs. It writes into `protection/results/`; preserve the selected reference snapshot with Git or run in a copied working directory. The fixed optimization seed aids reproducibility but floating-point/optimizer versions can change results. `python research/check.py` performs the safe checks in a temporary copy.

## Spectral transmission product

[spectra/transmission.py](spectra/transmission.py) writes
[spectra/shield_transmission.json](spectra/shield_transmission.json) (schema
`terluna.protection.shield-transmission/1`). It holds T(λ) from 200 nm to 5 µm
for the stored titania–silica design and for idealised edge filters. The design
is evaluated with `model.py`'s own thin-film code; the edge filters are
scenarios, not designs. The atmosphere domain reads this product to filter the
sunlight in its climate and photochemistry calculations. Regenerate it with
`python -m protection.spectra.transmission` after restoring the inputs.

## Short-wave transmission (2026-09-25)

[spectra/short_wave.py](spectra/short_wave.py) writes
[spectra/stack_short_wave.json](spectra/stack_short_wave.json) (schema
`terluna.protection.stack-short-wave/1`): the stored design's T(λ) from 0.1 to
210 nm, with the design's own layers and thin-film code. The optical constants
are published:
- CXRO atomic scattering factors (Henke, Gullikson and Davis 1993) for both
  oxides below 24.8 nm;
- fused silica from 24.8 nm (Franta et al. 2016, pinned in
  [inputs.json](inputs.json));
- titania from 120 nm (the design's Siefke data).

Between 24.8 and 120 nm no measured titania set was found, so the CXRO
estimate stands in; the 10 µm of silica alone is opaque there (below 10⁻¹⁴⁶),
so this does not affect the result. The rebuilt layer sequence reproduces
`model.optical_stack` to 10⁻¹², and the X-ray absorption reproduces the design's
[xray_absorption.csv](results/xray_absorption.csv).

The film passes hard X-rays only:

| Wavelength | Transmission |
|---|---|
| 0.1 nm | 96% |
| 0.5 nm | 9% |
| 1 nm | 0.5% |
| 1.5 nm | 2×10⁻⁷ |
| Anything beyond 2.5 nm | 4×10⁻⁸ at most |
| 5–200 nm | below 10⁻²⁰ |

The atmosphere domain's escape calculation reads this. The film lets at most
10⁻⁹ W/m² of heat into the upper air for the quiet Sun and 10⁻⁷ W/m² at solar
maximum, so the exobase stays within 1.5 K of its base. How much extreme
ultraviolet reaches the Moon is therefore set by light that bypasses the film:
gaps, pinholes and edges of the aperture, and off-normal incidence, which
[transmission/](transmission/README.md) now screens (below). The sky adds heat of its own:
interplanetary hydrogen glows in Lyman-alpha at about 1,000 rayleigh and reaches
the upper air from every direction, where no Sun-facing shield can block it
([report.md](report.md), section 9; a first estimate is in the atmosphere
domain's [middle-atmosphere README](../atmosphere/middle_atmosphere/README.md)).
Regenerate the product with `python -m protection.spectra.short_wave` after
restoring the inputs.

## Annulus films (2026-10-07)

[spectra/annulus_film.py](spectra/annulus_film.py) writes
[spectra/annulus_film.json](spectra/annulus_film.json) (schema
`terluna.protection.annulus-film/1`). Beyond the climate window, the shield's
aperture out to four lunar radii has to stop the ultraviolet and pass visible
light. These films use the same two oxides: titania 0.02–1 µm thick on silica
0.5–10 µm thick, bare and with the stored stack's anti-reflection and matching
layers. Each is evaluated from 0.1 to 2,500 nm with `model.py`'s thin-film code
and short_wave's optical constants below 210 nm. It is weighted by the
quiet-Sun WHI 2008 spectrum below 202 nm and by TSIS-1 above.

| Coated film | Areal mass | Worst, 10–175 nm | Sunlight below 175 nm | Visible passed |
|---|---:|---:|---:|---:|
| 0.1 µm titania on 1 µm silica | 3.6 g/m² | 4×10⁻⁷ | 1.1×10⁻⁴ | 98.6% |
| 0.1 µm titania on 2 µm silica | 5.8 g/m² | 1.4×10⁻⁷ | 3.7×10⁻⁵ | 98.6% |
| 0.1 µm titania on 4 µm silica | 10.2 g/m² | 1.3×10⁻⁷ | 1.0×10⁻⁵ | 98.7% |
| Stored stack, 1 µm on 10 µm | 26.9 g/m² | 2×10⁻⁴³ | 6.7×10⁻⁷ | 97.8% |

- **Ultraviolet.** A tenth of a micrometre of titania stops the extreme and far
  ultraviolet.
- **Soft X-rays.** Between 2.5 and 10 nm, soft X-rays set the silica.
- **Hard X-rays.** Below 2.5 nm, hard X-rays pass every film lighter than the
  stored stack, so they decide each light film's share of the sunlight below
  175 nm.
- **Upper-air heat.** The product also gives each film's heat of the upper air.
  It uses the atmosphere domain's film-heat count, with every transmitted X-ray
  absorbed above the base, at quiet Sun and at solar maximum. Beside it is the
  same model's heat for light that passes gaps at O1's swarm levels.
- **Light pressure.** The reflected and absorbed shares give the light pressure,
  (2R + A) times S/c: 0.12–0.16 for coated films and about 0.5 for bare ones.
- **Not modelled.** Normal incidence only. Gaps, pinholes, the support's
  strength and ageing are not included.

Regenerate it with `python -m protection.spectra.annulus_film` after restoring
both this domain's inputs and the WHI spectrum
(`python -m atmosphere.middle_atmosphere.fetch_inputs --download`).

## The September design report (imported 2026-09-26)

[report.md](report.md) is the package's full design report of 9 September 2026,
verbatim except for two private conversation links removed from its first
source. It gives the reasoning behind the tables: the Moon-following screen and
its roughly 976,000 replaceable 10×10 km cells, the power and propellant
closure, exhaust isolation (section 8), residual heating including the sky's
Lyman-alpha (section 9), the regional magnets, corridors and cosmic rays,
manufacturing and billion-year renewal, a staged schedule, decision gates, and
the alternatives near Sun–Earth L1 (section 18).

Its architecture keeps the held screen light and puts mass elsewhere: a
solar-filter complex (optical swarm, power and propulsion, metrology, stores and
servicing); an industrial hub for fabrication, propellant conditioning, freight
and repair; lunar orbital depots; heavy inhabited and industrial facilities on
economical trajectories that service the moving screen; and the original
Earth–Sun L1/L2 hubs for interplanetary freight, power and industry.

The author's requirements, recorded in the September feasibility and protection
reports:

- construction within 500 years of 2026, and operation for at least 10⁹ years;
- an open atmosphere with no pressure dome, and no superconducting planetary ring;
- an exobase near 250 K preferred, with 260 K a candidate limit;
- a surface radiation dose of at most 0.027 mSv/day;
- protection that expands around occupied destinations and Earth–Moon traffic
  corridors, secondary to protecting the lunar atmosphere and surface.

[reference/historical/](reference/historical/) holds reconstructions of older
concepts whose original files are unrecovered. The recovery pass in the
author's other project space wrote them on 2026-09-26 from the 2025
conversation records, so their values are as recalled there and none is an
original file. They cover Lunashield-L1 (2025), a 23–40 t plasma-inflated
mini-magnetosphere for the solar wind, with a consolidated record of its system
design dossier; the Earth–Sun L1 super hub EL1-SH (2025), with an EUV metascreen
of tiles, rafts and veils, transparent photovoltaics capped at 1–2% of sunlight,
power beaming and a mature hub near 300 TW, also with a consolidated record; and
the lineage from the 2025 plasma ideas to this design. The Lunashield record
keeps the staged idea for charged-particle protection: seed stations with local
cavities, then more nodes whose cavities grow and merge until they fill a
protected volume, then redundancy. Paths inside the lineage refer to the
recovery dump's layout ([research/archive_status.md](../research/archive_status.md)).

The loss response's exosphere step
([atmosphere/loss_response](../atmosphere/loss_response/README.md)) sizes the
magnets from the air's side. Without a magnetosphere the solar wind carries off
the ions the sunlit exosphere makes, and its charge exchange with the dense
exosphere near the exobase costs about 0.7–3 kg/s that no optical shadow
removes, so a 1 kg/s budget needs magnetic protection. A lunar dipole of about
3×10¹⁹ A·m², whose stand-off clears the dense exosphere of cool upper air,
removes that loss in the screening; the four regional installations give
1.5×10²¹ A·m² and a stand-off near 10 lunar radii. Neither holds the neutral
fragments of molecules that sunlight breaks up outside the shadow, which the
optical shield's reach has to cover.

## Module catalogue

[modules/catalogue.py](modules/catalogue.py) writes
[modules/catalogue.json](modules/catalogue.json) (schema
`terluna.protection.module-catalogue/1`): each replaceable unit of the September
reference design with its count, mass, materials, power, consumables and
lifetime scenarios. There are about 976,000 optical cells of 10 × 10 km at
5×10⁶ kg each, their share of the holding hardware (1.1×10⁶ kg, 183 MW mean power
and 0.29 kg/s of propellant per cell) and of the propellant buffer, and the four
regional magnet installations. The values come only from this model's stored
outputs; the engineering domain's supply ledger
([engineering/network](../engineering/network/README.md)) reads the catalogue as
demand. Regenerate it with `python -m protection.modules.catalogue`.

## UV transmission of the swarm (2026-09-26)

[transmission/model.py](transmission/model.py) writes
[transmission/results/uv_transmission.json](transmission/results/uv_transmission.json)
(schema `terluna.protection.uv-transmission/1`): the share of sunlight below
175 nm that reaches the protected region through a formation of 10 × 10 km
cells. It counts seams, micrometeoroid holes, missing cells and manufacturing
defects for three levels of design choices. Tight choices hold 3–4×10⁻⁵,
standard ones about 2×10⁻⁴ and relaxed ones 2.5×10⁻³. Failed cells waiting to
be covered make up most of it, and overlapping seams and patching every 10–20
years keep the rest small. The design levels are assumptions to test;
[transmission/README.md](transmission/README.md) has the results and limits.

## Condition and remaining questions

- Optical layers use measured constituent constants plus effective-medium and normal-incidence approximations. This is not a measured complete coating, irradiated lifetime test, or gap-free aperture.
- The holding-force calculation is an instantaneous circular-phase sweep, not full ephemeris propagation or a passive lunar orbit.
- Specific electrical power, material properties, fuel buffer, exhaust cant and storage capability are assumed design parameters. Clean plume operation remains unestablished. Report section 8 sets the condition: for a 10 kg/s loss allowance, less than about 94 MW of plume power (7.5×10⁻⁷ of the 125 TW exhaust) may reach the upper air, and it lists fallbacks from thrusters outside the footprint to externally delivered momentum. Whether the solar wind, which blows from the screen toward the Moon, carries ionized exhaust into the upper air is open.
- Holding the screen takes 2.8×10⁵ kg/s of propellant. The feasibility report's energy-limited loss with no filter at all is about 5,300–195,000 kg/s at 10% efficiency ([baseline report](../research/baselines/feasibility/report.md), section 6), so the reference screen spends more mass than it saves. It saves mass only if it is much lighter, uses much faster exhaust or is held without propellant.
- The hub architecture above (industrial hub, depots, Earth–Sun L1/L2 hubs) has no design or budget here yet.
- Magnetic moment, pressure-balance and rigidity estimates are component diagnostics. Global plasma performance, reconnection, trapped particles and human radiation dose remain unmodelled.
- Thermospheric chemistry, lower climate, water loss and species escape remain to be coupled. The original model does not establish that its filter produces a 250 K exobase. The film's own short-wave transmission is now computed (above), and the atmosphere domain finds that it holds the exobase at 142–193 K before the sky's Lyman-alpha glow is counted. Light bypassing the film through the aperture is screened for assumed design choices (above); formation control, cell failure rates and the oxide film's hole size are untested.
- Hardware replacement may be much smaller than construction throughput while accumulated propellant/resource use remains substantial over geological time.

The preserved results include successful and failed propulsion closures. Large sampled optical/phase grids and the rendered PDF are omitted from this curated import and recorded in [provenance](../research/provenance.json). Their numerical tables can be regenerated from the original code once inputs are restored. The older [shield geometry table](reference/legacy_shield_geometry.csv) remains explicitly separate from this later design.

[Checks](../research/checks.json) report precisely which implementation tests were run, without granting environmental or engineering validation.

## Full-ephemeris array study (2026-10-04)

[The solar shield/habitat study](../research/studies/solar_shield_array/report.md)
adds [DE440 dynamics](dynamics/README.md), conserved photon-momentum bounds,
a bounded monthly trajectory search, area quadrature and 19-year eclipse
checks. The 50-g/m² reference reduces mean holding power from about 309 to
238 TW on the selected moving path, with about 371,000 kg/s of propellant.
A 10-g/m² replacement film gives about 40 TW on that path; the stored oxides
alone already weigh about 26 g/m². The imported September model stays intact.
Optical material performance, safe plumes, tile formation, passive solutions
and lifetime material closure remain open.

The [cycling follow-up](../research/studies/solar_shield_array/cycling.md)
keeps 50 g/m² and lets gravity return tiles between shadow-service passes.
A 15,000-km seed stays bounded for three years in the ideal point-sail model,
with about 15.7% useful projected area and no electric thrust. A separate
optimized solar-sail return arc is independently replayed. These results
support a fleet-inventory alternative; spatial coverage, real spectral and
attitude response, seams, collision avoidance and payload mass remain open.

The [simultaneous fleet follow-up](../research/studies/solar_shield_array/fleet.md)
tests finite squares in several planes and radial families over a common month.
Solar-disk ray unions expose substantial gaps even when summed areas exceed
the protected aperture. Swept collision and mutual-shadow exclusions reduce
the candidate populations further. The 50 m placement gate fails on independent
replay, and large-square probes need extended-body dynamics. Four lunar radii
is the primary target; continuous coverage and payload power remain unresolved.

The [active-assistance tests](../research/studies/solar_shield_array/active.md)
use separate 10 km tiles. Bounded feedback closes local finite-Sun gaps while
preserving physical separation, and independently replayed phase transfers
move a tile's projected arrival by about 1,400 km with 5.45–6.79 m/s of correction.
Their energy and propellant budgets are explicit. A global assignment and its
recurring handover/return costs remain to be solved.
