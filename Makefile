# Affected work: `make check`. Every repository check: `make check-full`.
PYTHON ?= python3
CHECK_ARGS ?=
export OPENBLAS_NUM_THREADS ?= 1
export OMP_NUM_THREADS ?= 1
export MKL_NUM_THREADS ?= 1

.PHONY: check check-plan check-full check-layers test-python test-js test-immersion test-column test-light-transport test-column-viewer test-calendar check-provenance check-ensemble build-immersion

check:
	$(PYTHON) shared/check_changed.py $(CHECK_ARGS)

check-plan:
	$(PYTHON) shared/check_changed.py --plan $(CHECK_ARGS)

check-full: check-layers test-python test-js test-calendar check-provenance check-ensemble

# Imports respect the lanes (shared, domains, research, visualization, immersion).
check-layers:
	$(PYTHON) shared/check_layers.py

# Domain, study and shared tests, and the Jura Island flyover's (the sky solver tests run when numba is installed).
test-python:
	$(PYTHON) -m pytest

# The immersion (after its bake), the atmospheric column model, the light-transport references
# and the column viewer.
test-js: test-immersion test-column test-light-transport test-column-viewer

test-immersion:
	cd immersion && npm test

test-column:
	node --test atmosphere/column/tests/

test-light-transport:
	node --test illumination/references/tests/

test-column-viewer:
	node --test visualization/atmospheric-columns/tests/

test-calendar:
	$(PYTHON) visualization/light-calendar/build.py
	node --test illumination/calendar/tests/*.test.mjs visualization/light-calendar/tests/*.test.mjs

# Byte-pinned historical imports and the reproduced feasibility baseline.
check-provenance:
	$(PYTHON) research/check.py

check-ensemble:
	$(PYTHON) ensemble/tools/validate_setup.py

build-immersion:
	cd immersion && npm run build
