# Lovelace's Square Code Library

This repository collects MATLAB implementations of common preprocessing and analysis routines used in chemometrics. Every algorithm resides under `Codes/` in its own folder containing the main function, a test script, and a README.

## Quick Reference

The table below lists the available algorithms and a short description of each.

| Folder | Purpose |
|-------|---------|
| AsLs | Asymmetric Least Squares baseline correction |
| BIRFI LS | IRF estimation with Tikhonov regularization |
| BIRFI LS smoothing | IRF estimation with additional smoothing penalty |
| EMSC | Extended Multiplicative Scatter Correction |
| Fourier filter | Frequency-domain filtering of time-series data |
| I-SVD | Iterative SVD-based PCA imputation |
| Kernelize | Convolve signals with predefined kernels |
| LALS | Local Asymmetric Least Squares baseline correction |
| MSC | Multiplicative Scatter Correction |
| O ALS (2) | Orthogonalized Alternating Least Squares |
| PCA filter | Low-rank filtering via PCA |
| SavGol | Savitzky–Golay smoothing and derivatives |
| WienerFiltering | Wiener filtering for spectra |
| autoscale | Autoscale by rows or columns |
| binning v 2.0 | Bin N‑D arrays by summing or averaging |
| corr_map | Plot correlation heatmap |
| cosmicpeakcorrection | Remove cosmic spikes from spectra |
| cropBackground | Remove low/high intensity pixels from image cubes |
| normMatrix | Normalize matrices using several norms |
| simulate_spectra | Generate synthetic Gaussian spectra |
| unfoldImage | Flatten a 3‑D image cube into 2‑D |
| whittakerSmoother | Whittaker smoothing |
| whittakerSmoother imputation | Whittaker smoothing with imputation |

For usage details consult the `README.txt` inside each folder.
