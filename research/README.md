# Terluna research

Shared computational and theoretical support for the five-paper [Open Moon ensemble](../ensemble/). **Constructing and Sustaining an Open Moon is the next full manuscript.** Companion analyses develop alongside it; the core's opening landscape does not define the limits of the research.

| Topic | Available material | Current condition |
|---|---|---|
| [Atmosphere](../atmosphere/) | April screening code; September hydrostatic/Jeans code and cases | Conditional models with different assumptions; compatibility with a realizable spectrum remains open |
| [Climate](../climate/) | 18-case linear thermal-response diagnostic | Executable scalar sensitivity; no spatial weather model |
| [Geography](../geography/) | Terrain-data leads and basin/hydrology plan | No elevation grid or simulation loaded |
| [Illumination](../illumination/) | Existing angular-twilight diagnostic | Runnable geometry arithmetic; no sky-brightness or spectral model |
| [Biosphere](../biosphere/) | Mechanisms, candidate biomes and empirical leads | No executable life-cycle/ecosystem model recovered |
| [Habitation](../habitation/) | Spatial/time-use concepts and mechanical relationships | No complete settlement, vehicle or services simulation |
| [Engineering](../engineering/) | Transfer, propulsion, growth and maintenance tables | Conditional accounting; complete fleet/manufacturing network missing |
| [Protection](../protection/) | Optical, holding-force, magnetic and renewal code/results | Component calculations; five pinned external inputs required to run |

Read [status.json](status.json) for condition and next-task detail, [plan.md](plan.md) for the core-first work order, and [archive_status.md](archive_status.md) for unresolved recovery. [provenance.json](provenance.json) lists every original archive member and its disposition, including hashes for retained files and omitted material. Model code and CSV tables are verbatim; the two JSON reference files explicitly select original fields with unchanged values.

## Reproduction

```sh
python -m pip install -r research/requirements.txt
python research/check.py
# Restore the exact five protection inputs from the original project ZIP:
python protection/fetch_inputs.py --archive /path/to/Lunar_Protection_Model.zip
python research/check.py --require-inputs
```

Upstream retrieval is also available through `python protection/fetch_inputs.py --download`; the exact expected hash and size must match. Direct network retrieval was unavailable in the preparation container and was not newly verified. With missing protection data the checker reports `BLOCKED` explicitly; it never substitutes synthetic inputs.

The [checked snapshot](checks.json) records an actual run using the original archived inputs. Future checks write to ignored `research/runs/` by default. Numerical reproduction and file integrity do not grant physical validation, full-source admission, or manuscript clearance.

## Organization

The September feasibility implementation is a historical multi-domain script. Its single intact copy lives in [baselines/feasibility](baselines/feasibility/); topic folders hold its original compact reference tables. Protection has its own intact implementation. Large generated grids, rendered figures/PDFs, duplicate baseline copies and external spectral bytes are omitted from the curated commit; their source records and reproduction/restoration paths remain explicit.

Ordinary file names and Git history manage ongoing changes. The earlier dated planning snapshot and accepted seeds stay intact. New calculations should state assumptions, track conservation/residuals, test convergence and identify which conclusion their outputs can change.
