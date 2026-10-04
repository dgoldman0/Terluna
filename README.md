# Terluna

Research toward constructing and sustaining an open, living Moon.

The project combines shared numerical research with a five-paper ensemble. The next full manuscript is **Constructing and Sustaining an Open Moon**: a grounded vision of a nearby world with diverse environments, living communities and ways of inhabiting space. Four companions develop the physical, biological, human and industrial arguments independently.

## Start here

- [Research inventory](research/README.md), [condition register](research/status.json) and [decisions register](research/decisions.md)
- [The environment screens' equations, bounds and limits](research/findings.md)
- [Five-paper ensemble and accepted seeds](ensemble/README.md)
- [Core-first research plan](research/plan.md)

## Layout

| Lane | Folders | Purpose |
|---|---|---|
| Domains | `atmosphere/`, `climate/`, `biosphere/`, `geography/`, `illumination/`, `protection/`, `engineering/`, `habitation/` | Models, simulations and reference calculations, each with its tests |
| Research hub | [research/](research/README.md) | Findings, plan, status register, provenance, and cross-domain studies in `research/studies/` |
| Visualization | [visualization/](visualization/README.md) | Scientific and engineering rendering of domain results |
| Immersion | [immersion/](immersion/README.md) | The explorable experience: engine, world content, experiences; reads domain results as baked products |
| Shared | [shared/](shared/README.md) | Constants, scenarios, the data-product convention and the layer check |
| Papers | [ensemble/](ensemble/README.md) | The five-paper manuscript ensemble |
| Archive | [archive/](archive/README.md) | Deprecated material kept for provenance |

## Available research

| Area | Available work | Evidence boundary |
|---|---|---|
| [Atmosphere](atmosphere/) | A solved thermal column; line-by-line radiative–convective and photochemical columns under the shield (ozone, surface ultraviolet); the exobase, its heating and escape; the loss response to protection | Three-dimensional middle-atmosphere transport; light bypassing the shield film through aperture gaps |
| [Climate](climate/) | An ExoPlaSim GCM design climate near 299.6 K with lunar-gravity corrections; CM1 cloud-resolving rings and a highland box; [the seas' waves](climate/waves/) through two lunar cycles, from the open sea to run-up on the shore | Paused since 2026-09-30 while CM1 and the GCM disagree on warmth and humidity over land; the waves rest on one sampled month and a 1-degree sea grid |
| [Biosphere](biosphere/) | The carbon reserve theorem and model; an oxygen box; tree and forest-patch mechanics; canopy photosynthesis by wavelength; a plant's carbon through the lunar night; fruit designed for the lunar day | Measured traits, complete life cycles and ecological interactions |
| [Geography](geography/) | LOLA topography above the GRAIL geoid; the atlas at 28% water, named from the IAU gazetteer; rivers and rain-fed lakes; [every sea's monthly tide](geography/README.md#the-monthly-tide) | How the water divides between seas and lakes; groundwater; Mare Fecunditatis's strait |
| [Illumination](illumination/) | A spectral clear-sky solver and light-transport references; Sun, Earth and star geometry with earthlight; surface spectra; a shielded spherical sky; evening cloud radiance | Regional profiles, three-dimensional clouds, aerosols, polarization and the earthlight spectrum |
| [Protection](protection/) | The component model and its September design report; the film's transmission from hard X-rays to the far ultraviolet | Aperture leakage, clean operation, particle transport and lifetime resources |
| [Engineering](engineering/) | Resource transport, growth and renewal accounts; plume-heat sensitivity | A complete industrial network and safe source-to-use routes |
| [Habitation](habitation/) | Design concepts; a first optical-comfort screen | Practical inhabited capacity; visual comfort in resolved scenes and human validation |

Cross-domain studies, from conservation and the protection architecture to
lunar-cycle ecology, optical comfort and evening clouds, are listed in the
[research inventory](research/README.md#studies).

The environment screens' [results](research/studies/environment_screens/results/) contain numerical cases and approximation flags. Solving the selected equations establishes their conditional consequences. Environmental compatibility, biological persistence and engineering performance require additional evidence.

## Run

```sh
python -m pip install -r research/requirements.txt -r immersion/bake/requirements.txt
(cd immersion && npm install)
make check                                   # layer check, Python and JS tests, provenance, ensemble
OPENBLAS_NUM_THREADS=1 python -m research.studies.environment_screens.run
```

The new screens run with NumPy and SciPy. To inspect the original optical input ranges, restore the protection inputs or supply `--protection-archive /path/to/Lunar_Protection_Model.zip` to the screen runner. Without those inputs the audit explicitly records that their ranges were not inspected in that run. The band-heating interface rejects unspecified responses and spectral gaps.

Use ordinary filenames and Git history for ongoing research. The imported baselines, locked planning snapshot and accepted seeds retain their provenance. See [AGENTS.md](AGENTS.md) before making changes.

## Historical material

[Training data](archive/training_data/README.md) are deprecated synthetic worldbuilding records, retained for provenance. They are outside the scientific evidence and model-calibration workflow. Their historical log is preserved separately; no training records have been deleted.
