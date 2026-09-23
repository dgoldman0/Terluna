# Terluna

Research toward constructing and sustaining an open, living Moon.

The project combines shared numerical research with a five-paper ensemble. The next full manuscript is **Constructing and Sustaining an Open Moon**: a grounded vision of a nearby world with diverse environments, living communities and ways of inhabiting space. Four companions develop the physical, biological, human and industrial arguments independently.

## Start here

- [Current findings, equations and limits](research/findings.md)
- [Research inventory](research/README.md) and [condition register](research/status.json)
- [Five-paper ensemble and accepted seeds](ensemble/README.md)
- [Core-first research plan](research/plan.md)

## Available research

| Area | Available work | Evidence boundary |
|---|---|---|
| [Atmosphere](atmosphere/) | Historical baselines; a solved molecular thermal column; band-energy interface | Prescribed lower boundary; full spectral chemistry and kinetic escape remain open |
| [Climate](climate/) | Conservative, periodic latitude–longitude energy-balance screen | Synthetic geography and uncalibrated radiation/transport; no weather or ice model |
| [Biosphere](biosphere/) | Periodic carbon-storage theorem and model; oxygen budget with explicit shortfalls | Hypothetical functional traits; no life-cycle or ecosystem validation |
| [Protection](protection/) | Optical, positioning, magnetic and renewal estimates | Component calculations; important EUV response/input gap documented |
| [Engineering](engineering/) | Resource transport, growth and maintenance accounting | Conditional budgets; complete industrial and safety closure remains open |
| [Geography](geography/), [illumination](illumination/), [habitation](habitation/) | Data leads, angular geometry and design requirements | Actual-terrain climates, visual appearance and settlement capacity still require work |

The [results](research/studies/environment_screens/results/) contain numerical cases and approximation flags. Solving the selected equations establishes their conditional consequences. Environmental compatibility, biological persistence and engineering performance require additional evidence.

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
