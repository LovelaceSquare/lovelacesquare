# LAsLS: Local Asymmetric Least Squares Baseline Correction (MATLAB)

The Local Asymmetric Least Squares (LAsLS) algorithm extends asymmetric least squares (AsLS) smoothing to one-dimensional data vectors by allowing different asymmetry parameters (`pVals`) and second-derivative smoothing penalties (`lambdasAsym`) within specified intervals. Outside these intervals, a uniform smoothing penalty (`lambdaWhit`) and a global first-derivative penalty (`mu`) are enforced.

Through an Iteratively Reweighted Least Squares (IRLS) procedure, the method updates a weight matrix based on residuals to adaptively estimate a smooth, locally tailored baseline.

**Reference of the original AsLS**:
Eilers, P. H. C., & Boelens, H. F. M. (2005). *Baseline correction with asymmetric least squares smoothing*. Leiden University Medical Centre Report, 1(1), 5.

---

## Contents

Two packages are available under `Codes/`:

- **`LASLS_CL/`** — Command-line function for scripted/batch use
  - `LASLS_CL.m` — Standalone baseline correction function
  - `test_LASLS_CL.m` — Test script with synthetic data
- **`LASLS_GUI/`** — Interactive graphical interface
  - `LASLS.m` — Main app class (MATLAB AppBase + uihtml)
  - `LASLS_test.m` — Synthetic data generator for the GUI
  - `business_logic/` — Core algorithm (`@LASLSCorrector`) and validation (`@DataValidator`)
  - `ui/` — HTML/JS frontend

---

## Command-Line Usage (`LASLS_CL`)

    [estimatedBaseline, weights] = LASLS_CL(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);

See `test_LASLS_CL.m` for a full working example with synthetic data.

---

## GUI Usage (`LASLS_GUI`)

Launch the graphical interface from MATLAB:

    LASLS

### Quick start

1. **Load data** — Click *Load* to import a matrix variable from the workspace (rows = samples, columns = channels). To try the tool without your own data, use the *Generate demo data* option in the *+* menu.
2. **Draw intervals** — Click *Draw* and drag on the preview chart to define peak regions. Each interval gets its own local smoothing (λ) and asymmetry (p) parameters.
3. **Tune parameters** — Use the *Global* tab to adjust baseline smoothness (λ), asymmetry (p), and tension (μ). Switch to the *Interval* tab to fine-tune parameters for individual intervals.
4. **Preview** — The baseline estimate updates in real time as you adjust parameters. The corrected signal is shown in the bottom chart.
5. **Apply & Export** — Click *Apply* to compute baselines for all samples, then *Export* to save corrected data and baselines to the workspace.

### Features

- **Signal-by-signal navigation** — Browse individual samples using the badge controls in the chart header.
- **Per-signal parameters** — Assign different global parameters to each sample for heterogeneous datasets.
- **Peak detection** — Automatic peak finder to help define intervals.
- **Batch processing** — Apply correction to all samples at once.
- **Import/Export intervals** — Save and reload interval definitions.
- **Interactive zoom & pan** — Scroll to zoom, drag to pan on both charts.
- **Dark mode** — Toggle via the *+* menu.
- **Built-in tutorial** — Click *Tutorial* in the *+* menu for a step-by-step guided tour of all interface features.

---

## Installation

### Prerequisites
- MATLAB R2016a or later (uses sparse matrices and solvers)

### Setup

1. Add the command-line function to your MATLAB path:

        addpath('path/to/LASLS_CL');

2. For the GUI, add with subfolders:

        addpath(genpath('path/to/LASLS_GUI'));

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
