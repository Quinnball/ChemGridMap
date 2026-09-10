PYTHON ?= python3
export MPLCONFIGDIR := $(CURDIR)/.cache/matplotlib
VENV_PYTHON := .venv/bin/python

.PHONY: setup check test paper paper-legacy paper-current paper-application paper-record-tasks
setup:
	$(PYTHON) -c "import sys; assert sys.version_info >= (3, 10), 'Use Python 3.10 or newer, not the legacy chemicalmap environment'"
	$(PYTHON) -m venv .venv
	$(VENV_PYTHON) -m pip install -e '.[dev,umap,paper]'

check:
	$(VENV_PYTHON) scripts/check_environment.py

test: check
	$(VENV_PYTHON) -m pytest -q

paper: paper-current

paper-current: check
	$(VENV_PYTHON) scripts/restore_paper_data.py --download
	$(VENV_PYTHON) paper/validate_current.py
	$(VENV_PYTHON) paper/verify_source_context.py
	$(VENV_PYTHON) paper/prepare_tasks.py
	$(VENV_PYTHON) -m pytest -q --junitxml=paper/output/current/pytest.xml

paper-legacy: check
	$(VENV_PYTHON) paper/run_revision_validation.py
	$(VENV_PYTHON) paper/validate_adaptive_assignment.py
	$(VENV_PYTHON) paper/build_revision_figures.py
	$(VENV_PYTHON) paper/audit_revision_outputs.py

paper-application: check
	$(VENV_PYTHON) paper/run_application_validation.py
	$(VENV_PYTHON) paper/build_revision_figures.py --source paper/output/revision_v5 --output paper/output/revision_v5/figures
	$(VENV_PYTHON) paper/audit_revision_outputs.py --source paper/output/revision_v5

paper-record-tasks: check
	$(VENV_PYTHON) paper/run_record_tasks.py --source paper/output/revision_v5
	$(VENV_PYTHON) paper/build_revision_figures.py --source paper/output/revision_v5 --output paper/output/revision_v5/figures
