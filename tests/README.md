# Environment model tests

Run `python -m unittest discover -s tests -v` from the repository root. The tests use NumPy and SciPy and cover supplied-input validation, conservation, independent analytic comparisons, the reserve proof's discrete implementation, periodic boundaries and numerical resolution. No network access or external spectral archive is required for unit tests.

The recorded pass covers 39 tests. Model-domain warnings remain visible even when a test passes. Whole-repository historical checks are still performed separately by `research/check.py` on a complete checkout.
