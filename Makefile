PYTHON ?= python3
export MPLCONFIGDIR := $(CURDIR)/.cache/matplotlib
VENV_PYTHON := .venv/bin/python

.PHONY: setup check test paper paper-legacy paper-current paper-application paper-record-tasks paper-structure \
	paper-f0 paper-e5 paper-e5b paper-e5c paper-e5d paper-e6-prepare \
	paper-e2c paper-e6-regression paper-e6 paper-e7 paper-e9 \
	paper-e10 paper-e11 paper-e12a paper-e12b paper-posthoc paper-t1 \
	paper-fig2-data paper-figures-data paper-figures paper-figure6-capture \
	paper-reproduce-tests
setup:
	$(PYTHON) -c "import sys; assert sys.version_info >= (3, 10), 'Use Python 3.10 or newer, not the legacy chemicalmap environment'"
	$(PYTHON) -m venv .venv
	$(VENV_PYTHON) -m pip install -e '.[dev,umap,paper]'

check:
	$(VENV_PYTHON) scripts/check_environment.py

test: check
	$(VENV_PYTHON) -m pytest -q

paper: paper-current

paper-structure: check
	$(VENV_PYTHON) paper/validate_structural_comparison.py
	$(VENV_PYTHON) paper/build_structure_figures.py

paper-current: check
	$(VENV_PYTHON) scripts/restore_paper_data.py --download
	$(VENV_PYTHON) paper/validate_current.py
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

PAPER_SOURCE := paper/strengthening/source
PAPER_FIG2_INPUT := validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv
PAPER_FIG2 := $(VENV_PYTHON) -m chemgridmap.cli $(PAPER_FIG2_INPUT) \
	--input-format chembl --target-id CHEMBL205 --activity-type IC50 \
	--lower-threshold 6 --upper-threshold 7 --representation morgan \
	--projection umap --umap-metric jaccard --radius 2 --n-bits 2048 \
	--random-state 42 -k 10 --grid-occupancy 0.4 \
	--assignment-method dense --render-detail overview \
	--grid-refinement neighbors --core-min-coverage 0.8 \
	--support-level 0.8 --offset-display uncertain \
	--offset-min-shared 3 --offset-min-difference 1 \
	--offset-sign-agreement 0.8 \
	--map-color-modes activity class_support neighbor_recall structure_group \
	--output-formats svg png pdf --name caii_seed42

paper-f0: check
	$(VENV_PYTHON) $(PAPER_SOURCE)/f0_build.py

paper-e5: check
	$(VENV_PYTHON) $(PAPER_SOURCE)/e5_analysis.py

paper-e5b: paper-e5 paper-fig2-data
	$(VENV_PYTHON) $(PAPER_SOURCE)/e5b_analysis.py

paper-e5c: paper-e5b
	$(VENV_PYTHON) $(PAPER_SOURCE)/e5c_analysis.py

paper-e5d: paper-e5c
	$(VENV_PYTHON) $(PAPER_SOURCE)/e5d_analysis.py

paper-e2c: check
	$(VENV_PYTHON) $(PAPER_SOURCE)/e2c_analysis.py

paper-e6-prepare: paper-e2c
	$(VENV_PYTHON) $(PAPER_SOURCE)/e6_prepare.py

paper-e6-regression: paper-e6-prepare
	$(VENV_PYTHON) $(PAPER_SOURCE)/e6_regression_diag.py

paper-e6: paper-e6-prepare paper-e6-regression
	$(VENV_PYTHON) $(PAPER_SOURCE)/e6_analysis.py

paper-e7: paper-e2c
	$(VENV_PYTHON) $(PAPER_SOURCE)/e7_analysis.py

paper-e9: paper-e5d
	$(VENV_PYTHON) $(PAPER_SOURCE)/e9_analysis.py

paper-e10: paper-e6
	$(VENV_PYTHON) $(PAPER_SOURCE)/e10_analysis.py

paper-e11: paper-e5c
	$(VENV_PYTHON) $(PAPER_SOURCE)/e11_analysis.py

paper-e12a: paper-e11 paper-e9
	$(VENV_PYTHON) $(PAPER_SOURCE)/e12a_rule_calibration.py

paper-e12b: paper-e7 paper-e9
	$(VENV_PYTHON) $(PAPER_SOURCE)/e12b_analysis.py

paper-posthoc: paper-e5 paper-e5b paper-e5c paper-e5d paper-e6 paper-e7 paper-e9 paper-e10 paper-e11 paper-e12a paper-e12b

paper-t1: check
	$(VENV_PYTHON) $(PAPER_SOURCE)/t1_scale.py --mode all

paper-fig2-data: check
	$(PAPER_FIG2) --core-rule coverage --fmcs-timeout 1 \
		--output-dir paper/figures/data_fig2/default
	$(PAPER_FIG2) --core-rule fmcs --fmcs-timeout 10 \
		--output-dir paper/figures/data_fig2/fmcs

paper-figures-data: paper-posthoc paper-f0
	$(VENV_PYTHON) paper/figures/source/data_final.py

paper-figure6-capture: check
	node paper/figures/source/fig6_capture.cjs \
		--output-dir paper/figures/data_fig6

paper-figures: paper-figures-data paper-fig2-data
	$(VENV_PYTHON) paper/figures/source/build_figures.py

paper-reproduce-tests: check
	$(VENV_PYTHON) -m pytest paper/strengthening/source paper/figures/source
