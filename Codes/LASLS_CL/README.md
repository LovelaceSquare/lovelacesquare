# LAsLS Command-Line: Local Asymmetric Least Squares Baseline Correction

Standalone MATLAB function for Local Asymmetric Least Squares (LAsLS) baseline correction. LAsLS extends Asymmetric Least Squares (AsLS) by allowing different asymmetry parameters (`pVals`) and smoothing penalties (`lambdasAsym`) within user-defined intervals. Outside these intervals, symmetric weighting and a uniform smoothing penalty (`lambdaWhit`) are applied. A global first-derivative penalty (`mu`) enforces baseline tension across the full signal.

The baseline is estimated via Iteratively Reweighted Least Squares (IRLS):

```
(W + D' * diag(lambdaVec) * D + mu * L' * L) * baseline = W * y
```

**Reference**: Eilers, P. H. C., & Boelens, H. F. M. (2005). *Baseline correction with asymmetric least squares smoothing*. Leiden University Medical Centre Report, 1(1), 5.

---

## Usage

```matlab
[baseline, weights] = LASLS_CL(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);
```

### Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `y` | (n x 1) vector | Input signal |
| `intervals` | (m x 2) matrix or cell | Each row defines `[startIdx, endIdx]` of a peak region |
| `pVals` | (m x 1) vector | Asymmetry parameter for each interval (small values ignore peaks) |
| `lambdasAsym` | (m x 1) vector | Local smoothing penalty for each interval |
| `lambdaWhit` | scalar | Smoothing penalty outside intervals |
| `mu` | scalar | Global first-derivative penalty (baseline tension) |
| `maxIter` | scalar | Maximum IRLS iterations |
| `tol` | scalar | Convergence tolerance (relative baseline change) |

### Outputs

| Output | Type | Description |
|--------|------|-------------|
| `baseline` | (n x 1) vector | Estimated baseline |
| `weights` | (n x 1) vector | Final IRLS weights |

### Example

```matlab
y = rand(100,1)*10;
intervals = [10 20; 40 50];
pVals = [0.01; 0.001];
lambdasAsym = [1e4; 2e5];
lambdaWhit = 100;
mu = 1e3;
maxIter = 50;
tol = 1e-6;
[baseline, weights] = LASLS_CL(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);
```

See `test_LASLS_CL.m` for a complete working example with synthetic data and visualization.

---

## Contents

| File | Description |
|------|-------------|
| `LASLS_CL.m` | LAsLS baseline correction function |
| `test_LASLS_CL.m` | Test script with synthetic data, plots, and MSE evaluation |

---

## Installation

Add to your MATLAB path:

```matlab
addpath('path/to/LASLS_CL');
```

### Dependencies

Built-in MATLAB functions only: `spdiags`, matrix division (`\`), `norm`.

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrian Gomez-Sanchez** -- Universitat de Barcelona & Lovelace's Square, Barcelona, Spain
- **Berta Torres-Cobos** -- University of Copenhagen, Denmark & Lovelace's Square, Barcelona, Spain
- **Rodrigo Rocha de Oliveira** -- Universitat de Barcelona & Lovelace's Square, Barcelona, Spain
