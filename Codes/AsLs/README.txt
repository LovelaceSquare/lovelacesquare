# AsLS

Asymmetric Least Squares (AsLS) baseline correction for spectral data. The
method models a slowly varying background beneath each spectrum by solving a
weighted least squares problem. Weights are updated iteratively so that large
positive residuals (peaks) are penalized more heavily than negative residuals.
The smoothing parameter `lambda` controls how stiff the estimated baseline is
while the asymmetry parameter `p` tunes the emphasis on negative excursions.
Typical values range from 1e4 to 1e8 for `lambda` and around 0.001 for `p`.

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

## What the code does
Asymmetric Least Squares (AsLS) baseline correction for spectral data. The

## How to use it
Run `AsLS.m` in MATLAB. A basic demonstration is provided in `test_AsLS.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_AsLS.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
