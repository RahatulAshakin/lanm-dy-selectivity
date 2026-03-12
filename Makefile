PYTHON ?= python
PIP ?= $(PYTHON) -m pip

.PHONY: setup audit fetch annotate test

setup:
	$(PIP) install -e .[phase1,dev]

audit:
	$(PYTHON) -m lanm.cli.audit_data

fetch:
	$(PYTHON) -m lanm.cli.fetch_public_structures

annotate:
	$(PYTHON) -m lanm.cli.annotate_metal_sites

test:
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py'
