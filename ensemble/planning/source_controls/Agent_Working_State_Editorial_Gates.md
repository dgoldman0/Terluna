# Editorial gates for scholarly systems papers

**Control revision:** control-01, 19 September 2026  
**Current manuscript status:** editorial clearance blocked; substantive revision deferred.

These controls consolidate the author's existing paragraph-inheritance, reorderability, audience, and affirmative-development requirements with the failures identified in the preceding review. They govern future revisions of this project and can be carried into later paper projects with an updated audience and technical brief. They establish release criteria; this control pass performs no manuscript revision.

## Scope and precedence

The intended audience works in computer architecture, systems, managed runtimes, or agent infrastructure. Assume familiarity with memory hierarchies, allocation and garbage collection, asynchronous execution, and host/device offload. Explain a concept where its precise use introduces a new contract, ambiguity, or consequence for this proposal. Calibrate any future field adaptation through the same audience test.

Keep the architecture generic, the argument affirmative, the benefits conditional on evidence, and the paper within its agreed short-paper scope. Preserve the current author, linked ORCID, and scholarly preamble unless a later request changes them. A technical change opens the affected semantic and editorial checks; a formatting change opens the affected rendered-page checks. Later user instructions take precedence within their stated scope.

## Gate G1 — An argument instead of an expanded inventory

**Required development.** Each prose paragraph advances a question, operation, causal relationship, evidentiary comparison, or design consequence. Its governing relation determines the order of its sentences and the point at which the next paragraph begins. A section develops its question through such paragraphs, carrying a consequence into the following section.

**Test.** Reduce the paragraph to its main assertions. State why that ordering matters, then try swapping its middle assertions or moving its last assertion earlier. For an expository paragraph with independently meaningful items, identify the reason for keeping them in prose. Put genuine catalogs in a table or concise list when that format better serves the argument.

**Blocking failure.** A feature tour, definition series, or evaluation checklist survives with substantially arbitrary order. Repeated references to a shared object supply topical continuity while leaving the reasoning undeveloped. Adding *while*, *as*, semicolons, or relative clauses leaves the inventory intact.

**Clearance evidence.** Identify the paragraph, record the relationship that controls its order, and explain what logical work a swap would lose. For a repaired passage, include a short before/after excerpt and identify the changed inference. Simply restating its topic or listing the edits leaves the gate open.

## Gate G2 — Contribution appropriate to the readership

**Required development.** Established systems concepts enter through a proposal-specific question: a placement trade-off, execution cost, lifetime dependency, failure case, or interface obligation. A paragraph earns its space by investigating that question. Definitions remain available where the chosen meaning changes the design or its evaluation.

**Test.** Delete a familiar explanatory sentence mentally and ask what a systems reader loses: a necessary assumption, a distinctive mechanism, an inferential step, or merely a reminder of established knowledge. Apply the same test to sentences that explain the significance of the immediately preceding passage. Also test headings and captions for unnecessary reader instruction.

**Blocking failure.** The text explains that backing storage extends fast memory, that completion precedes safe buffer reuse, or that end-to-end measurements include communication, then presents the explanation itself as a contribution. A paragraph ends by announcing an implication already fully expressed. Reader-management language tells the audience what to recognize, remember, or conclude.

**Clearance evidence.** Locate the suspect passage and name the proposal-specific work it now performs, or record its removal with the continuity preserved. Retained foundational explanation requires a concrete ambiguity or dependency, not a generic appeal to accessibility. Technical precision remains compatible with an expert audience.

## Gate G3 — Inheritance, syntactic hierarchy, and cadence

**Required development.** Except at a purposeful turn, each sentence receives an actor, state, condition, unresolved question, or result from its predecessor and transforms it. Syntax expresses causal, temporal, constitutive, and evidentiary relations where they actually exist; coordination serves propositions of equal force. Paragraph endings complete a development or generate the next question.

**Test.** Trace what each sentence inherits and changes. Check clause-to-clause relationships within long sentences as carefully as the stops between short ones. Repeated subject–verb resets, overloaded lists inside subordinate clauses, stock transitions, and abstract recap landings trigger inspection. Read the passage continuously for cadence after its argument is sound.

**Blocking failure.** A chain of equal-weight declaratives persists, including inside longer sentences. Length or connective density increases while the dependency structure remains unchanged. Every paragraph finishes with a generalized lesson. A conclusion merely repeats the section sequence.

**Clearance evidence.** Record the inherited element and its changed role at the previously failing boundary, together with the landing's contribution. Describe the cadence check actually performed; a silent reread must be recorded as silent, and a read-aloud claim requires that activity. Short sentences remain legitimate hinges and landings, while sentence-length targets and compulsory subordination have no role as pass criteria.

## Gate G4 — Affirmative scholarly stance and economy

**Required development.** Claims unfold through the proposal, the evidence, and the conditions under which the proposal could work. Qualifications attach to the claim they qualify. The prose uses precise operations and physical or logical consequences, with self-reference reserved for useful methodological disclosure.

**Test.** Inspect anticipatory rebuttals, objection catalogs, self-praise, generic importance claims, surplus explanation, and inflated contribution verbs. Check the manuscript for the specifically excluded phrases `not` and `rather than`; inspect other negations and imperative reader instructions in context. Captions and the abstract receive the same review as body prose.

**Blocking failure.** Defensive framing organizes the argument; a lexical substitution preserves its defensive shape; elementary implications are presented as discoveries; a conditional benefit becomes a performance claim during stylistic revision. Promotional framing or elaborate connective prose obscures the substantive question.

**Clearance evidence.** Record actual findings and their disposition. Lexical scans are screening evidence only. Preserve exact bibliographic titles, faithful source quotations, code, and mathematically necessary expressions; log any literal exception by location and reason. When accurate technical meaning requires an exception in authorial prose, surface it for resolution instead of concealing it with a misleading paraphrase.

## Gate G5 — Technical meaning and evaluation discrimination

**Required development.** An editorial change preserves the architecture's ownership, authority, lifetime, recovery, and concurrency conditions. Hardware mechanisms remain distinct from software policy, and a conceptual separation receives only the capacity or performance implications that its implementation and evidence support.

**Test.** Check each substantive edit against its originating claim and downstream uses, including captions and equations. For an evaluation item, name the hypothesis it discriminates, the competing explanation or baseline, the held-constant semantics, and the observation that would change the architectural judgment. Shared acceptance tests can cover several mechanisms where that relationship is explicit.

**Blocking failure.** Smoother prose silently adds concurrency, isolation, unlimited capacity, automatic semantic validity, or exactly-once effects. Generic metrics and familiar benchmarks substitute for a discriminating experiment. Tags are treated as authority, collector reachability as semantic relevance, or fast device kernels as established whole-system gains.

**Clearance evidence.** Preserve a compact semantic-diff record and a hypothesis-to-test mapping. Record changed assumptions explicitly. Retain rejection regions and costs within the analysis; their technical role supports an affirmative account of where specialization earns its place.

## Gate G6 — Full-document source admission

**Required evidence.** Every cited source has an accessible complete local copy, its exact version and checksum, a completed full-document reading record, and claim-level page or section locators. The reading record covers appendices and the supplied reference list as well as the main text, and identifies actual quality or scope concerns. Reading a source qualifies access; claim fit and evidentiary quality qualify its use.

**Blocking failure.** A discovery snippet, abstract, selected passage, or browser-only access report substitutes for the requested full download and reading. A prior review label silently becomes fresh verification. Access failures become an implicit waiver of the source requirement.

**Clearance evidence.** A source manifest identifies the archived artifact and exact version, with the reading record and claim ledger linked to it. Record proposals and deductions as such, separately from sourced findings. Unresolved access holds source clearance until the source is obtained or the affected citation and claim are revised under an authorized content pass. This control pass leaves the existing source-access problem open.

## Gate G7 — Exact-build visual and production review

**Required evidence.** The final compiled PDF is identified by its hash and compared with the source version being released. Inspect every page at readable scale, then inspect each figure, table, and equation closely. Check borders, labels, arrows, captions, pagination, links, embedded fonts, and author metadata. Previously problematic vector diagrams receive a second-renderer check.

**Blocking failure.** A warning-free compile, text extraction, a contact sheet alone, or a prior build's images substitute for final-page inspection. Label backgrounds hide borders; opposed transfers share illegible labels; a final edit changes wrapping after review. Render claims refer to an unspecified build.

**Clearance evidence.** Save page-specific observations and the rendered image artifacts, recording the exact PDF hash and renderer. A changed build reopens this gate. Preamble quality, citation resolution, and visual correctness remain separate from editorial quality.

## Conduct of an authorized revision

Work one section at a time, beginning with its question and reading adjacent sections before changing prose. Establish the governing relation of a failing paragraph, repair its argument, then address audience calibration and sentence movement. Read the revised section continuously, check its technical dependencies, and update its evidence record before advancing. Revisit the abstract, transitions, captions, and conclusion after the body has stabilized.

A section record can be brief. It needs the exact source hash, locations reviewed, the paragraph movement, the relevant deletion/reorderability findings, repairs to previously failed passages, and dispositions under G1–G4. Explicitly mark any paragraph or caption intentionally functioning as a catalog or method statement. Coverage concerns the entire section, while supporting excerpts can focus on its difficult passages. A new substantive section requires its own record.

After local review, use a separate review pass to challenge the manuscript's hardest passages against G1–G5. A separate pass by the same reviewer remains a separate pass; it carries no claim of independent peer review. A later change reopens the affected section and adjacent dependencies, with a final continuous read checking the complete argument.

## Release authority and record integrity

The active record is `review/release-status.json`. Each gate has a status of `PASS`, `BLOCKED`, or `DEFERRED`, a reviewer, substantive evidence files, and the artifact hash to which the decision applies. G1–G4 additionally require clearance for the abstract/front matter and each body section. A `PASS` is an evidenced judgment, not a description of intended work.

The release check requires every gate and section result to pass, requires current source hashes and accessible evidence files, and requires the final PDF hash for visual clearance. Automated validation establishes record completeness and currency; it cannot establish the truth of a prose judgment or the adequacy of source reading. A mechanical success supports the review record and never replaces substantive examination.

A failure holds an **edited/final** delivery. A requested draft may remain available with accurate status. When the user defers editing, retain the issue and stop at the control record. Neither a wording detector nor a polished layout can clear an argument or audience failure, and a disclaimer leaves an unresolved source gate unresolved.

## Calibration cases carried forward from the preceding assessment

These are existing findings, retained as future regression checks. The manuscript remains unchanged during this control pass.

- **Section 2, memory-tier paragraph:** the assertion that retention/residency separation supplies an enterprise-scale hierarchy culminates in a familiar backing-store implication. G2 requires concrete placement consequences to justify the space; G5 preserves the distinction between conceptual separation and achieved scale.
- **Section 3, opening hardware paragraph:** subordinate clauses join a tour of tags, allocation, barriers, and collection engines. G1 requires a real ordering relation; G3 checks that the revised syntax carries that relation rather than concealing another inventory.
- **Section 3.2, bridge/reclamation landing:** the final sentence repeats the preceding transfer-lifetime implication. G2 and G3 require a consequential landing or an economical ending after the result already established.
- **Section 6, evaluation sequence:** familiar baselines, replay modes, ablations, and metrics require a clearer relationship to competing architectural predictions. G1 and G5 assess that relationship.
- **Section 4, frozen-view passage:** serialization, approval, transmission, and response binding offer an existing example of an order grounded in operation dependencies. Use its relation as a calibration example, without treating the paragraph or section as pre-cleared.
