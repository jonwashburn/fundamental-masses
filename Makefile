# Makefile for Particle Masses Empirical Regularity Validation
# This reproduces all results in the manuscript

.PHONY: all clean test rg residues certs ablations figures paper

# Python interpreter
PYTHON := python3

# Directories
CODE_DIR := code
SCRIPTS_DIR := $(CODE_DIR)/scripts
CORE_DIR := $(CODE_DIR)/core
OUT_DIR := out
CSV_DIR := $(OUT_DIR)/csv
TEX_DIR := $(OUT_DIR)/tex
FIG_DIR := $(OUT_DIR)/figs
DATA_DIR := data

# Create output directories
$(shell mkdir -p $(CSV_DIR) $(TEX_DIR) $(FIG_DIR))

# Default target - run everything
all: clean rg residues certs test ablations figures paper
	@echo "=========================================="
	@echo "COMPLETE PIPELINE EXECUTED SUCCESSFULLY"
	@echo "=========================================="
	@echo "Results in: $(OUT_DIR)/"
	@echo "Paper: Fundamental-Masses-Revised.pdf"

# Clean previous outputs
clean:
	@echo "Cleaning previous outputs..."
	@rm -f $(CSV_DIR)/*.csv
	@rm -f $(TEX_DIR)/*.tex
	@rm -f $(FIG_DIR)/*.pdf
	@rm -f $(OUT_DIR)/*.log

# Step 1: RG transport (PDG → μ*)
rg:
	@echo "=========================================="
	@echo "Step 1: RG Transport PDG → μ*"
	@echo "=========================================="
	$(PYTHON) $(CORE_DIR)/quark_rg.py \
		--mustar 182.201 \
		--output $(CSV_DIR)/transported_masses.csv | tee $(OUT_DIR)/rg.log
	@echo "✓ Transport complete"

# Step 2: Compute residues
residues:
	@echo "=========================================="
	@echo "Step 2: Computing Residues at μ*"
	@echo "=========================================="
	$(PYTHON) $(CORE_DIR)/pm_rs_native_full.py \
		--input $(CSV_DIR)/transported_masses.csv \
		--output $(CSV_DIR)/residues_at_mustar.csv | tee $(OUT_DIR)/residues.log
	@echo "✓ Residues computed"

# Step 3: Certificate validation
certs:
	@echo "=========================================="
	@echo "Step 3: Certificate Validation"
	@echo "=========================================="
	$(PYTHON) $(SCRIPTS_DIR)/validate_certificates.py \
		--input $(CSV_DIR)/residues_at_mustar.csv \
		--output $(CSV_DIR)/certificates.csv | tee $(OUT_DIR)/certificates.log
	@echo "✓ Certificates validated"

# Step 4: Main test (pass/fail decision)
test: certs
	@echo "=========================================="
	@echo "Step 4: Main Validation Test"
	@echo "=========================================="
	@$(PYTHON) -c "import pandas as pd; \
		df = pd.read_csv('$(CSV_DIR)/certificates.csv'); \
		all_pass = df['passes'].all(); \
		print('OVERALL RESULT:', 'PASS ✓' if all_pass else 'FAIL ✗'); \
		exit(0 if all_pass else 1)"

# Step 5: Ablation studies
ablations:
	@echo "=========================================="
	@echo "Step 5: Ablation Studies"
	@echo "=========================================="
	$(PYTHON) $(SCRIPTS_DIR)/validate_certificates.py --ablations \
		--output $(CSV_DIR)/ablations.csv | tee $(OUT_DIR)/ablations.log
	@echo "✓ Ablations complete"

# Step 6: Generate figures
figures:
	@echo "=========================================="
	@echo "Step 6: Generating Figures"
	@echo "=========================================="
	$(PYTHON) $(SCRIPTS_DIR)/make_figures.py \
		--residues $(CSV_DIR)/residues_at_mustar.csv \
		--ablations $(CSV_DIR)/ablations.csv \
		--output $(FIG_DIR)/ | tee $(OUT_DIR)/figures.log
	@echo "✓ Figures generated"

# Step 7: Generate LaTeX tables
tables:
	@echo "=========================================="
	@echo "Step 7: Generating LaTeX Tables"  
	@echo "=========================================="
	$(PYTHON) $(SCRIPTS_DIR)/make_all_masses_table.py \
		--input $(CSV_DIR)/residues_at_mustar.csv \
		--output $(TEX_DIR)/validation_table.tex
	@echo "✓ Tables generated"

# Step 8: Compile paper
paper: figures tables
	@echo "=========================================="
	@echo "Step 8: Compiling Paper"
	@echo "=========================================="
	pdflatex -interaction=nonstopmode Fundamental-Masses-Revised.tex
	bibtex Fundamental-Masses-Revised || true
	pdflatex -interaction=nonstopmode Fundamental-Masses-Revised.tex
	pdflatex -interaction=nonstopmode Fundamental-Masses-Revised.tex
	@echo "✓ Paper compiled: Fundamental-Masses-Revised.pdf"

# Verification targets
verify-lean:
	@echo "Verifying Lean formalization..."
	cd . && lake build
	@echo "✓ Lean verification complete"

verify-python:
	@echo "Verifying Python environment..."
	$(PYTHON) -c "import numpy, pandas, scipy; print('✓ Python packages OK')"

# Help target
help:
	@echo "Particle Masses Validation Pipeline"
	@echo "===================================="
	@echo "Targets:"
	@echo "  all       - Run complete pipeline (default)"
	@echo "  clean     - Remove all generated files"
	@echo "  rg        - Run RG transport PDG → μ*"
	@echo "  residues  - Compute residues at anchor"
	@echo "  certs     - Validate certificates"
	@echo "  test      - Run main validation test"
	@echo "  ablations - Run ablation studies"
	@echo "  figures   - Generate all figures"
	@echo "  tables    - Generate LaTeX tables"
	@echo "  paper     - Compile the paper"
	@echo ""
	@echo "Verification:"
	@echo "  verify-lean   - Check Lean formalization"
	@echo "  verify-python - Check Python environment"

# Installation helpers
install-deps:
	pip install -r requirements.txt
	lake exe cache get  # For Lean/Mathlib

# Archive for Zenodo
archive:
	@echo "Creating archive for Zenodo..."
	git archive --format=zip --prefix=particle-masses/ HEAD > particle-masses.zip
	@echo "✓ Archive created: particle-masses.zip"
