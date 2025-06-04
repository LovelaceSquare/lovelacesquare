# Lovelace's Square Code Library

This repository collects MATLAB implementations of common preprocessing and analysis methods used in chemometrics. Each algorithm lives in its own folder under `Codes/` and includes:

- The main `.m` function implementing the algorithm.
- A `test.m` script that demonstrates basic usage.
- A `README.txt` describing the purpose, usage and authorship information.

If you are exploring the project with an LLM, this file provides a quick overview. See `lovelace_square_readme.md` for complete contribution guidelines and coding standards.

## Available Algorithms

| Folder | Description |
|-------|-------------|
| AsLs | AsLS. Perform Asymmetric Least Squares (AsLS) baseline correction on |
| BIRFI LS | BIRFI_LS. Perform IRF estimation using Tikhonov regularization. |
| BIRFI LS smoothing | BIRFI_LS_SMOOTHING. Perform IRF estimation with smoothing using Tikhonov regularization. |
| EMSC | EMSC.  Perform Extended Multiplicative Scatter Correction on a 2D data matrix. |
| Fourier filter | fourierFilter. Apply frequency-domain filtering to time-series signals. |
| I-SVD | I_SVD  Perform Iterative SVD-based PCA Imputation on the input data. |
| Kernelize | kernelize Apply kernelization to the input matrix D |
| LALS | LALS. Perform Local Asymmetric Least Squares (LALS) baseline correction with per-interval parameters. |
| MSC | MSC.  Perform Multiplicative Scatter Correction on a 2D data matrix. |
| O ALS (2) | OALS  Perform Orthogonalized Alternating Least Squares (OALS) on the input data. |
| PCA filter | pcaFilter Performs PCA filtering on the input matrix. |
| SavGol | SavGol.  Apply Savitzky-Golay filtering to a 2D data matrix. |
| WienerFiltering | WienerFiltering. Apply Wiener filtering to spectral data. |
| autoscale | AUTOSCALE Perform autoscaling on a 2D data matrix. |
| binning v 2.0 | binning Bins the data array according to specified bin sizes and mode. |
| corr_map | Plot the correlation map |
| cosmicpeakcorrection | COSMICPEAKCORRECTION  Removes cosmic spikes from spectral data by |
| cropBackground | CROPBACKGROUND Crop background pixels from a 3D image cube by intensity range. |
| normMatrix | normMatrix. Normalizes the input matrix using specified norms. |
| simulate_spectra | Simulate spectra with n gaussian without noise. |
| unfoldImage | UNFOLDIMAGE. Flatten a 3D image cube into a 2D matrix – because chemometrics |
| whittakerSmoother | WhittakerSmoother. Apply the Whittaker smoother for signal smoothing. |
| whittakerSmoother imputation | WHITTAKERSMOOTHERIMPUTE  Apply the Whittaker smoother for signal smoothing, |

For more details on each algorithm, see the README.txt inside each folder.
