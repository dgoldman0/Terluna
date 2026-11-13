# Extraction protection and scientific records

This study applies the project's existing preservation requirements to candidate lunar and near-Earth sources. The reader is a resource-model reviewer deciding what an inventory calculation establishes and what evidence must precede disturbance. The evaluator records whether each required condition has supporting attestations. An inventory, cost allowance or complete set of attestations confers no extraction authorization.

The scope is enclosed lunar communities and their actual supporting source and transport facilities. Open Moon transformations, the habitat array and outer-Solar-System development remain outside this calculation. The conservation framework supplies applicable principles and records; its atmospheric and flooding milestones are not a timetable for the enclosed phase.

## Authority and application

The conservation framework is dated 25 September 2026 and retains decided items, proposals and open questions separately. Its heritage evaluations are a first pass awaiting review, apart from individually decided treatments. The implementation preserves those distinctions in each gate's `authority_status`. It applies requirements of proposed origin as conservative conditions of this study; it does not change their standing in the project's decisions register.

| Condition applied to this study | Existing authority and exact locator | Implementation and limit |
|---|---|---|
| Ordinary sources; protected major bodies remain untouched | Decisions register, line 53, 9 October 2026; conservation principles 4–5, lines 63–66 | `ordinary_source_selection`; protected categories are Saturn's rings, Titan's atmosphere and ocean worlds. These category exclusions do not expand the enclosed project's source region. |
| Separate scientific, heritage and functional values | Conservation decisions D2–D3, lines 158–159; explanation, lines 72–91; register criteria, lines 179–212 | `scientific_heritage_review`, `scientific_function_continuity`. Samples and measurements cannot replace heritage whose significance depends on place. Any lost scientific function needs an at-least-as-useful replacement before loss. |
| Survey and preserve information before irreversible disturbance | Conservation principle 1, line 59; science programme, lines 283–292; gates, lines 301–313 | `pre_disturbance_survey`, `representative_samples`, `archive_before_disturbance`. The proposed principle is a conservative release condition. Landing, prospecting, dust and exhaust can trigger it before any atmospheric transformation. |
| Preserve samples under original conditions and keep material unopened | Conservation principle 2, line 60; methods, lines 138–141; sample requirements, lines 294–299 | `native_sample_conditions`, `unopened_samples`. Sealed vacuum, dry inert storage and native-temperature preservation apply where relevant. The protocol specifies an unopened part, without a numerical fraction. |
| Maintain samples and records across locations and institutions | Conservation principle 3, line 62; archive requirements, lines 294–299; stewardship, lines 343–347 | `multiple_archive_locations`, `institutional_continuity`. Copies in multiple places and maintenance by multiple organizations are independent conditions. |
| Preserve provenance, scientific review and access | Conservation lines 76–78, 294–299 and 343–347; open-record principle, line 67 | `recorded_provenance`, `open_scientific_review`, `open_record_access`. A stored file alone does not satisfy all three. |
| Reserve accessible future resources and recover scarce elements | Decisions register, line 599, 31 October 2026; world-economy development study, lines 64–80 | `future_accessible_stocks`, `scarce_element_recovery`. Future use, emergency stocks, protected material, uncertainty and recovery losses precede a surplus claim. The lunar direction is applied conservatively to source extraction. The development study's synthetic masses and percentages are not imported as approved reserves. |
| Preserve affected polar records and control contamination | Conservation D5, line 164; scientific programme, lines 277–285; local-disturbance gates, lines 303–313 | `polar_record_protection` applies only to lunar cases; `dust_exhaust_contamination` applies to both regions. An unaffected finding still needs a record. No fixed buffer distance or mining-area fraction is invented. |
| Maintain stewardship and resolve participation conditions | Conservation D8–D9, lines 170–171; proposed consent principle, line 68; register and stewardship, lines 205–212 and 343–348 | `consent_and_stewardship`. Stewardship is decided; the consent wording and equal standing of affected communities and relevant expertise retain their open/proposed status. |
| Preserve all six standing cargo-safety requirements | Decisions register, lines 84–95, confirmed 26 September 2026 | Independent gates cover passive failure, Earth-clear catcher geometry, remote braking, packet energy, tracking/interception/corridors and verified capture before a planetary intercept. The historical 10⁶–10⁸ kg mature-traffic range is illustrative; the exact risk-qualified energy cap remains open. |

The protection architecture's S4 invokes these traffic rules (requirements, lines 196–216). S7 also requires tracing outflows to their endpoints, including later Earth passages, but its atmospheric escape-region budget share is open. The enclosed study applies outflow/contamination assessment without transferring that numerical atmospheric budget. S1–S3 and S5–S6 describe a planetary protection architecture, construction sequence, duration and atmospheric tradeoff; this evaluator does not impose that architecture on enclosed communities.

The heritage-site register and the 31 mapped scientific targets are supporting records for review. Neither establishes a complete surveyed set of exclusion polygons. Proposed sealed reference areas, vault locations and register additions remain unresolved design decisions. A comprehensive registry and site-specific reviews must therefore accompany any eventual source qualification.

## Consistent source exclusions

The existing enclosed-source screen lists Bennu, Ryugu, Itokawa and Eros as illustrative scientific exclusions. Their JPL source IDs are retained here so a display-name change cannot bypass the screen. Matching either an identifier or the named object creates a hold. The same rule applies to every pool before union or subtraction. Removing a scientific object from a narrow pool must not release it into a broader remainder.

| Retained illustrative exclusion | JPL source ID |
|---|---|
| Bennu | 20101955 |
| Ryugu | 20162173 |
| Itokawa | 20025143 |
| Eros | 20000433 |

This four-object list is not a comprehensive or approved protected-source registry. Every other catalogue object remains unqualified until its physical, scientific, archive and governance conditions have supporting records. Explicit protection records and the adopted protected categories also create holds. Conflicting identifiers require reconciliation before set subtraction; any identifier matching a held object preserves the hold.

## Evaluator interface and reading rule

`protection.load_inputs()` reads the versioned input file and verifies every recorded source SHA-256. A changed source causes an error rather than silently changing the interpretation. Each gate records its exact source line range, requirement, scope and authority status. The evaluator has no numerical extraction quota and creates no legal permission.

`evaluate_release(subject, attestations=None, config=None)` requires `subject.region` to be `lunar` or `near_earth`. Attestations are keyed by gate ID, with an explicit status and a list of evidence-record references. `passed`, `failed` and `hold` require supporting references. Missing evidence is `unknown`; it is neither a demonstrated failure nor a pass. A condition outside the selected region is `not_applicable`. Inapplicability to a particular activity within a region needs an attested review, rather than omission of the gate.

| Returned `release_status` | Meaning |
|---|---|
| `unresolved` | At least one applicable condition lacks evidence and no explicit failure or hold has been recorded. |
| `blocked` | At least one applicable condition failed or is held, or the source has a protection hold. Other conditions may still be unknown. |
| `conditions_satisfied` | All applicable conditions have supplied pass attestations and no hold applies. The truth and sufficiency of those records are not independently verified by this calculation. |

`extraction_authorized` is always `False`, including the complete synthetic controls. The name `release_status` describes the status of this study's evidence conditions, not an operational permit. The current catalogue products supply no site-specific attestations for these gates, so their candidates remain unresolved or held.

`filter_source_sets(mapping, config=None)` returns `source_sets` and excluded records for each named pool. `source_difference(broad, narrow, config=None)` screens both pools, subtracts unambiguous stable source IDs, and returns `sources` plus exclusion records. Input IDs may be supplied as `source_id`, `spkid` or `id`; when more than one is supplied they must identify the same object for set subtraction.

```python
from resources.source_access import protection

config = protection.load_inputs()
status = protection.evaluate_release(
    {"source_id": "candidate-id", "region": "near_earth"},
    config=config,
)
assert status["release_status"] == "unresolved"
assert status["extraction_authorized"] is False
```

## Source identity and checks

The following files are local project records. Their bytes, authority descriptions and exact line locators are pinned in [protection_inputs.json](protection_inputs.json). No external paper or catalogue is copied into this module. The existing science-screen input also records the identity of its archived JPL catalogue.

| Source record | SHA-256 |
|---|---|
| `research/studies/conservation/README.md` | `850534d9f0615a95e8dd20618c7fda7480b2e96852f64c821e75fc33c2863631` |
| `research/decisions.md` | `f82f959838eebaeb1b86e560a48467f6c2e23e0782e969268b409992d420afa2` |
| `research/studies/world_economy/development.md` | `ec9c2034760f6b6f7b029b9f06eae03f448817c719d2ad7fc226a61c49a82e8f` |
| `research/studies/protection_architecture/requirements.md` | `b5c6d83181bd771fe21ac1a14b32dd0bd4118b6c272bb63c08fa25ca19acba25` |
| `research/studies/conservation/sites.json` | `cdfa3fcd37a32edccd0678cf7f698cde5f617de9188d258da62624a62b48ba26` |
| `research/studies/conservation/targets.json` | `0d12f8d6e570df0e2b316a11a500239d82c4c34d46d241bccff8d91d05bb5387` |
| `resources/enclosed_mining/inputs.json` | `b0ee1d80e96e1e91b764eae05d7db4d31fd733a834b0f2f339e07d7b5eb7b9d3` |

The focused test suite checks unknown versus failed/held conditions, every independent sample/archive requirement, synthetic complete attestations, explicit protection holds, identical exclusions across pools, protected identifier aliases, conflicting identifiers and changed source hashes. It does not validate a mine or inspect the truth of hypothetical attestation records.

Reproduce the focused checks from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest resources/source_access/tests/test_protection.py -q
```
