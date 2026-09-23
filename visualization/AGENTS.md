# Working on visualization

Visualization is scientific and engineering rendering: tools that show what a
model computed, faithfully and with its evidence boundary attached. Examples are
the month-of-light viewer of the clear-sky atlases, the A1–A3 accuracy labs and
the B1 reference path tracer.

- **Display, don't own.** Physics belongs to its domain. A tool here reads domain
  products or calls domain code; it does not become the source of truth. A
  display approximation (a fog formula, a tone curve) says that it is one.
- **Nothing depends on this lane.** Visualization may read any domain, research
  study or the immersion (to measure it), but no other lane imports from here.
- **Keep builds reproducible** from the products they read. Record product and
  source hashes in what you publish, and keep generated pages and images out of
  Git.
- **Faithfulness is the test.** Check that what is drawn matches the numbers,
  with round-trip values, reference images or stated tolerances.
