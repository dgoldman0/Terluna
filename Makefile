# Every repository check: `make check`. Each target also runs on its own.
PYTHON ?= python3
export OPENBLAS_NUM_THREADS ?= 1

.PHONY: check check-layers test-python test-js check-provenance check-ensemble build-immersion

check: check-layers test-python test-js check-provenance check-ensemble

# Imports respect the lanes (shared, domains, research, visualization, immersion).
check-layers:
	$(PYTHON) shared/check_layers.py

# Domain, study and shared tests (the sky solver tests run when numba is installed).
test-python:
	$(PYTHON) -m pytest

# The immersion (after its bake), the atmospheric column model and the light-transport references.
test-js:
	cd immersion && npm test
	node --test atmosphere/column/tests/ illumination/references/tests/

# Byte-pinned historical imports and the reproduced feasibility baseline.
check-provenance:
	$(PYTHON) research/check.py

check-ensemble:
	$(PYTHON) ensemble/tools/validate_setup.py

build-immersion:
	cd immersion && npm run build
