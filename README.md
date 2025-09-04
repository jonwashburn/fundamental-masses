# Empirical Regularity in Standard Model Fermion Masses

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the complete computational pipeline for the paper:

**"Empirical Regularity in Standard Model Fermion Masses: An Exact Integer Structure at a Universal Scale"**  
*Jonathan Washburn*

We report an empirical regularity: at μ* = 182.201 GeV, Standard Model fermion mass residues match a closed-form expression in integers determined by electric charge, verified to 10^-6 precision.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run complete validation
make all
```

## Key Result

At μ* = 182.201 GeV, the dimensionless residue R_i equals F(Z_i) where:
- F(Z) = ln(1 + Z/φ) / ln(φ)  
- Z is an integer determined by electric charge
- Verified to < 10^-6 for all quarks and charged leptons

## Repository Contents

- `Fundamental-Masses-Revised.pdf` - Full manuscript
- `Fundamental-Masses-Revised.tex` - LaTeX source
- `code/` - Python implementation
- `data/` - Input PDG masses  
- `Makefile` - Reproducibility pipeline
- `IndisputableMonolith.lean` - Lean formalization

## Reproducing Results

```bash
make clean    # Clean previous outputs
make all      # Run complete pipeline
make test     # Verify certificates
```

## Citation

```bibtex
@article{washburn2024fermion,
  title={Empirical Regularity in Standard Model Fermion Masses},
  author={Washburn, Jonathan},
  journal={Physical Review D},
  year={2024},
  note={Submitted}
}
```

## License

MIT License - See [LICENSE](LICENSE) file for details.