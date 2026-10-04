# Shared foundation

What every lane may use, and the rules that keep the lanes apart.

## Constants

[constants.json](constants.json) is the one source for physical and astronomical
constants (Boltzmann, Avogadro, lunar and terrestrial radius and GM, synodic and
sidereal months, the solar constant). [constants.py](constants.py) and
[constants.js](constants.js) read it, so Python models and the browser use
identical numbers, and derived values (lunar surface gravity GM/R², the gas
constant) are computed rather than copied.

One known inconsistency is recorded instead of hidden: the sky solver, the
atmospheric column model and the immersion use a rounded lunar gravity of
1.62 m/s², while the other lunar models use GM/R² = 1.6242 m/s². Moving them over
means regenerating the sky atlas, the column set and the A1–A3 references, so
`LEGACY_MOON_GRAVITY` names the rounded value until that decision is made.

Historical, byte-pinned models (the April simulator, the September feasibility
and protection imports) keep their own literals; [research/check.py](../research/check.py)
reproduces them as they were.

The `model_closures` section records CM1 r22's dry-air and water-vapour gas
constants for interpreting its surface thermodynamics in the wave forcing
export. These values apply to CM1's thermodynamic closure; the other atmospheric
models use gas constants for their specified compositions.

## Scenarios

[scenarios/sites.json](scenarios/sites.json) holds sites where models and experiences
evaluate the sky. The `development` site (0° N, 65° W) is a placeholder chosen so
Earth is visible about 25° up in the east while the Sun still crosses the zenith;
it is not a proposed Open Moon location.

[scenarios/water.json](scenarios/water.json) holds the Open Moon's standing-water cover: 28% of the
surface under seas and lakes, selected on 25 September 2026 within the author's 25–35% range. It
records the basis for the choice. The geography [atlas](../geography/README.md#atlas-at-the-selected-water-share)
reads the share, and the water to deliver follows from the domains' hydrology.

[scenarios/waves.json](scenarios/waves.json) records the 1,025 kg/m³ seawater
density used by the wind-wave studies and their shoreline continuation.
It is a calculation scenario; ocean composition and density variations remain
separate research inputs.

## Lanes and the layer check

| Lane | Folders | May import |
|---|---|---|
| shared | `shared/` | nothing else |
| domain | `atmosphere/`, `climate/`, `biosphere/`, `protection/`, `engineering/`, `geography/`, `habitation/`, `illumination/` | shared, other domains |
| research | `research/` (hub and studies) | shared, domains |
| visualization | `visualization/` | shared, domains, research, immersion (to measure it) |
| immersion engine | `immersion/engine/` | shared, the engine (never a world) |
| immersion world | `immersion/world/` | shared, the engine, worlds |
| immersion | `immersion/experiences/`, tests | shared, engine, worlds |
| immersion bake | `immersion/bake/` | shared, all of immersion, domains |

Nothing imports `archive/`. [check_layers.py](check_layers.py) enforces the table
for Python imports and JavaScript import/require/URL specifiers; `make check`
runs it.

## Data products

When one lane's results feed another, they travel as a data product rather than
an import. A product is a file whose header states:

- `schema`: a versioned name such as `terluna.atmosphere.column-set/1`;
- `producer`: the domain, the model file and its hash, and the exporter;
- `evidence`: what the numbers are and are not (for example "selected column
  experiments, not forecasts");
- units and the rule for reading it (for example the interpolation rule).

Consumers check the schema and record the hash of the product they used. The
[column set](../atmosphere/column/export_columns.cjs) is the current example; the
immersion's [bake step](../immersion/bake/) packs products into runtime assets.

A producer records the shared constants its code reads by name and value
instead of pinning `constants.json` or `constants.py` by hash:
[provenance.py](provenance.py) finds them from the `shared.constants` imports of
the producer's pinned Python files. The product's tests compare the recorded
values with the current ones, so adding a constant anywhere leaves products
current while a change to a value one of them used shows up. This follows
workflow tools such as DVC and Snakemake, which track the parameters a step
uses. The optical, sky and cloud-twilight producers follow it; the wave, tide
and coastal products still record the constants file's hash as a record of
their run, without a test of it.
