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

The Earth's visual geometric albedo is 0.242, from Robinson et al.'s (2025)
synthesis of the observations, with their analytic fit to its phase curve
(`visual_phase_*` in [constants.json](constants.json)). It replaced the long-quoted
0.367 on 2026-10-04; [illumination/earthlight](../illumination/earthlight/README.md),
illumination/ephemeris.py and the immersion follow the same curve.

Historical, byte-pinned models (the April simulator, the September feasibility
and protection imports) keep their own literals; [research/check.py](../research/check.py)
reproduces them as they were.

The `model_closures` section records CM1 r22's dry-air and water-vapour gas
constants for interpreting its surface thermodynamics in the wave forcing
export. These values apply to CM1's thermodynamic closure; the other atmospheric
models use gas constants for their specified compositions. It also records the
rain fall-speed law of CM1's Morrison microphysics, which
[atmosphere/rain](../atmosphere/rain/README.md) checks against its single-drop
solver, and the scheme's Earth limits on rain's size, which
[atmosphere/rain/breakup.py](../atmosphere/rain/breakup.py) carries to lunar
gravity.

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
| domain | `atmosphere/`, `climate/`, `biosphere/`, `protection/`, `engineering/`, `geography/`, `geology/`, `habitation/`, `humanities/`, `illumination/`, `provisioning/`, `resources/` | shared, other domains |
| research | `research/` (hub and studies) | shared, domains |
| visualization | `visualization/` | shared, domains, research, immersion (to measure it) |
| immersion engine | `immersion/engine/` | shared, the engine (never a world) |
| immersion world | `immersion/world/` | shared, the engine, worlds |
| immersion | `immersion/experiences/`, tests | shared, engine, worlds |
| immersion bake | `immersion/bake/` | shared, all of immersion, domains |

Nothing imports `archive/`. [check_layers.py](check_layers.py) enforces the table
for Python imports and JavaScript import/require/URL specifiers. Local checks
run it when source files change; the full suite always runs it.

## Checks

`make check` verifies the affected work before a local commit. It compares the
working tree with `HEAD`, including staged, unstaged and untracked files. It
checks changed JSON, TOML and Python syntax and the tracked diff's whitespace,
then runs the selected tests and validators. `make check-plan` prints the
commands and dependency paths without running them. A clean tree produces no
test selection; use a base revision when checking work already committed:

```sh
make check-plan CHECK_ARGS='--base main'
make check CHECK_ARGS='--base <revision-before-the-work>'
make check CHECK_ARGS='--include research/runs/example/forcing.json'
make check-full
```

[check_changed.py](check_changed.py) follows imports and literal file references
through the current source and the paths in product metadata, including paths
used as JSON keys. It follows downstream consumers transitively. A reference
inside a single Python test can select that test alone; a module-level input
or import selects its test file. Test fixtures in `conftest.py` affect the tests
below them. Uncovered code or data changes fall back to their domain or study's
configured tests. Changes to shared scientific foundations, Python collection
configuration or unrecognized repository-wide inputs request the full suite.

The selector reads source without importing models or regenerating products.
It caches dependency extraction under `.pytest_cache/check_changed/` and hashes
the source bytes again on each invocation; changes to the inventory or selector
invalidate the resolved graph. Explanatory citations in prose and JSON records
do not become input dependencies; actual file reads and pinned identities do.
It cannot establish every dependency assembled at runtime, supplied by an
external program or kept outside Git. Review its printed scope and include
known affected inputs or run the relevant explicit test command when those
dependencies change. A reported lack of configured coverage requires the
artifact's documented checks and inspection. A passing selected check reports
only its selected scope. Product freshness tests can require a separate
product rebuild; selecting tests never performs that rebuild automatically.

`make check-full` preserves the complete configured suite. Use it for substantial
integrations, changes with unbounded effects and explicit full verification.
Its component targets also remain available: `test-python`, `test-js`,
`test-immersion`, `test-column`, `test-light-transport`, `test-column-viewer`,
`test-calendar`, `check-layers`, `check-provenance` and `check-ensemble`.
The runners use one numerical thread by default and launch checks sequentially.

Successful results can cover later local commits if their relevant code, tests,
inputs, installed dependencies and environment remain unchanged. Record what
ran and reuse that evidence when only an unrelated record changes. The tool
does not cache test successes: ignored scientific inputs can change independently
of Git, so an unchanged commit identifier is insufficient evidence for reuse.
For selection audits, `python3 shared/check_changed.py --json --files <paths>`
prints a plan for explicitly named paths instead of the Git diff.

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
