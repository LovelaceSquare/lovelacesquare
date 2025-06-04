# AsLS

Asymmetric Least Squares (AsLS) baseline correction for spectral data.

## Inputs
- `data` (matrix): spectral data with one spectrum per row.
- `lambda` (double): smoothing parameter controlling baseline smoothness.
- `p` (double): asymmetry parameter (0 < p < 1) weighting positive/negative residuals.

## Outputs
- `correctedData` (matrix): baseline-corrected spectra.
- `baseline` (matrix): estimated baseline for each row.

## Example
```matlab
lambdaVal = 1e6;
pVal = 0.001;
[corrected, base] = AsLS(X, lambdaVal, pVal);
```

## Author
Adrián Gómez-Sánchez

Licensed under the MIT License.
