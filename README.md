# Standard Model Fermion Masses: Empirical Integer Structure

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Lean 4](https://img.shields.io/badge/Lean-4.11.0-green.svg)](https://leanprover.github.io/)

## Overview

This repository contains the complete computational pipeline for reproducing the results in:

**"Empirical Regularity in Standard Model Fermion Masses: An Exact Integer Structure at a Universal Scale"**

We report an empirical regularity: at μ* = 182.201 GeV, Standard Model fermion mass residues match a closed-form expression in integers determined by electric charge, verified to 10^-6 precision using standard QCD/QED renormalization group methods.

## Quick Start

```bash
# Clone repository
git clone https://github.com/[username]/particle-masses
cd particle-masses

# Install dependencies
pip install -r requirements.txt

# Run complete validation pipeline
make all

# Results appear in out/
```

## Key Results

At the universal scale μ* = 182.201 GeV:

| Family | Z value | Residue | Precision |
|--------|---------|---------|-----------|
| Up-type quarks (u,c,t) | 276 | 4.33521 | < 10^-6 |
| Down-type quarks (d,s,b) | 24 | 2.20671 | < 10^-6 |
| Charged leptons (e,μ,τ) | 1332 | 5.77135 | < 10^-6 |

## Repository Structure

```
particle-masses/
├── code/                      # Core computational code
│   ├── core/                  # RG evolution and residue calculations
│   │   ├── quark_rg.py       # QCD/QED running
│   │   ├── pm_rs_native_full.py  # Residue computation
│   │   └── driver.py         # Main orchestration
│   └── scripts/              # Analysis scripts
│       ├── validate_certificates.py  # Certificate validation
│       ├── make_all_masses_table.py  # Generate LaTeX tables
│       └── make_figures.py   # Generate plots
├── data/                     # Input data
│   ├── quarks_muStar.csv    # PDG quark masses
│   └── leptons.csv          # PDG lepton masses
├── out/                      # Output directory (generated)
│   ├── csv/                 # Numerical results
│   ├── tex/                 # LaTeX tables
│   └── figs/                # Figures
├── IndisputableMonolith.lean # Lean formalization
├── Makefile                  # Reproducibility orchestration
├── requirements.txt          # Python dependencies
└── README.md                # This file
```

## Validation Protocol

The validation follows a pre-registered protocol with no parameter adjustment:

1. **Fixed parameters**: μ* = 182.201 GeV, λ = ln(φ), κ = φ
2. **Integer map**: Z = 4+(6Q)²+(6Q)⁴ (quarks), (6Q)²+(6Q)⁴ (leptons)
3. **Certificate test**: R_i ∈ I^res_i and I^res_i ⊆ I^gap_Z_i

## Running Individual Components

```bash
# Step 1: RG transport PDG → μ*
make rg

# Step 2: Compute residues
make residues

# Step 3: Validate certificates
make certs

# Step 4: Run ablations
make ablations

# Step 5: Generate figures
make figures
```

## Ablation Studies

Pre-registered ablations demonstrate specificity of the integer structure:

| Modification | Max Error | Result |
|-------------|-----------|---------|
| Original Z map | < 10^-6 | ✓ Valid |
| Drop +4 (quarks) | 0.127 | ✗ Breaks |
| Drop (6Q)⁴ term | 0.238 | ✗ Breaks |
| Replace 6Q → 5Q | 0.315 | ✗ Breaks |

## Python Dependencies

- `numpy >= 1.21.0`
- `pandas >= 1.3.0`
- `scipy >= 1.7.0`
- `matplotlib >= 3.4.0`

Install all dependencies:
```bash
pip install -r requirements.txt
```

## Lean Formalization

The repository includes a Lean 4 formalization of the mathematical foundations:

```bash
# Build Lean project
lake build

# The formalization is in IndisputableMonolith.lean
```

Note: The Lean file formalizes mathematical structures but does NOT axiomatize the physics. The empirical validation is done through computational certificates.

## Reproducing Paper Results

To exactly reproduce all results in the paper:

```bash
# Clean any previous runs
make clean

# Run complete pipeline
make all

# Compile paper (requires LaTeX)
make paper
```

Output files:
- `out/csv/certificates.csv` - Certificate validation results
- `out/csv/ablations.csv` - Ablation study results
- `out/tex/validation_table.tex` - LaTeX table for paper
- `Fundamental-Masses-Revised.pdf` - Compiled paper

## Verification Checksums

To verify output integrity:

```bash
# Generate checksums
sha256sum out/csv/*.csv > checksums.txt

# Expected values (examples)
# certificates.csv: [sha256]
# ablations.csv: [sha256]
```

## FAQ

**Q: What is being claimed?**
A: An empirical regularity in SM fermion masses - at μ* the residues match a simple integer formula to 10^-6 precision.

**Q: Is this a new theory?**
A: No. This is an empirical observation within standard QFT using established RG methods.

**Q: How is this non-circular?**
A: PDG masses are transported to μ* for comparison. No measured mass appears on the right side of its own equation.

**Q: What are the key predictions?**
A: Equal-Z degeneracy (e.g., all up-type quarks have identical residues) and specific mass ratios at μ*.

## Citation

If you use this code or results, please cite:

```bibtex
@article{washburn2024fermion,
  title={Empirical Regularity in Standard Model Fermion Masses: 
         An Exact Integer Structure at a Universal Scale},
  author={Washburn, Jonathan},
  journal={[Submitted to Phys. Rev. D]},
  year={2024},
  doi={[pending]},
  archivePrefix={arXiv},
  eprint={[pending]}
}
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## Contact

- Author: Jonathan Washburn
- Email: jon@recognitionphysics.org
- Issues: Please use GitHub issues for questions/bugs

## Acknowledgments

We thank [acknowledgments] for helpful discussions.

---

**Repository archived at Zenodo**: DOI 10.5281/zenodo.XXXXXXX
