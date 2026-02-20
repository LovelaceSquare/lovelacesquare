# LAsLS: Local Asymmetric Least Squares Baseline Correction (MATLAB)

The Local Asymmetric Least Squares (LAsLS) algorithm extends asymmetric least squares (AsLS) smoothing to one-dimensional data vectors by allowing different asymmetry parameters (`pVals`) and second-derivative smoothing penalties (`lambdasAsym`) within specified intervals. Outside these intervals, a uniform smoothing penalty (`lambdaWhit`) and a global first-derivative penalty (`mu`) are enforced.

Through an Iteratively Reweighted Least Squares (IRLS) procedure, the method updates a weight matrix based on residuals to adaptively estimate a smooth, locally tailored baseline.

**Reference of the original AsLS**:
Eilers, P. H. C., & Boelens, H. F. M. (2005). *Baseline correction with asymmetric least squares smoothing*. Leiden University Medical Centre Report, 1(1), 5.

---

## Contents

- **`LASLS_CL.m`** — Command-line function for LAsLS baseline correction
- **`test_LASLS_CL.m`** — Test script with synthetic data for the command-line function
- **`LASLS_GUI/`** — Interactive graphical interface (MATLAB AppBase + uihtml)
  - `LASLS.m` — Main app class
  - `LASLS_test.m` — Synthetic data generator for the GUI
  - `business_logic/` — Core algorithm (`@LASLSCorrector`) and validation (`@DataValidator`)
  - `ui/` — HTML/JS frontend

---

## Overview

The `LASLS_CL` function implements a Local Asymmetric Least Squares (LAsLS) baseline correction algorithm tailored for one-dimensional data vectors, such as spectra or chromatograms.

Key features:
- **Interval-specific asymmetry (`pVals`)**: Different sensitivity to positive/negative residuals across intervals.
- **Interval-specific smoothing (`lambdasAsym`)**: Finer curvature control in defined segments.
- **Global penalties**:
  - `lambdaWhit`: uniform second-derivative penalty outside intervals
  - `mu`: global first-derivative penalty discouraging abrupt slope changes

The method solves the following weighted system in each IRLS iteration:

    (W + Dᵀ * diag(lambdaVec) * D + μ * Lᵀ * L) * b = W * y

Where:
- `y`: observed data vector
- `W`: diagonal weight matrix updated per iteration
- `D`, `L`: second- and first-order difference operators
- `lambdaVec`: vector combining local and global smoothing terms

---

## Usage Example

    [estimatedBaseline, weights] = LASLS_CL(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);

See `test_LASLS_CL.m` for a full working example with synthetic data.

---

## Installation

### Prerequisites
- MATLAB R2016a or later (uses sparse matrices and solvers)

### Setup

1. Add the folder to your MATLAB path:

        addpath('path/to/LAsLS');

2. For the GUI, also add subfolders:

        addpath(genpath('path/to/LAsLS/LASLS_GUI'));

### Dependencies
- Built-in MATLAB functions only: `spdiags`, matrix division (`\`), `norm`, `ones`, `length`

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez** — Universitat de Barcelona & Lovelace's Square, Barcelona, Spain
- **Berta Torres-Cobos** — University of Copenhagen, Denmark & Lovelace's Square, Barcelona, Spain
- **Rodrigo Rocha de Oliveira** — Universitat de Barcelona & Lovelace's Square, Barcelona, Spain

---

## Changelog

- **v1.1 (2025)**: Renamed to LASLS_CL/LASLS_GUI; added co-authors; added GUI
- **v1.0 (2024-12-16)**: Initial implementation of per-interval LAsLS baseline correction

---

## Keywords

- Baseline correction
- Asymmetric least squares
- AsLS
- Local asymmetric least squares
- LAsLS
- Spectral preprocessing
- MATLAB
- Detrending
