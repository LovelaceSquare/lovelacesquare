# MSC

Classic Multiplicative Scatter Correction using a reference spectrum. Each
spectrum is linearly regressed onto the reference and then rescaled by the
resulting slope and intercept, thereby removing additive and multiplicative
scatter effects. The reference can be the mean spectrum or a chosen sample and
should represent the ideal shape of the data.

## Inputs
- `data` (matrix): spectra with one row per sample.
- `refType` ("Mean Spectrum" or "Reference Index"): how to choose the reference.
- `refIndex` (int): row index when `refType` is "Reference Index".

## Outputs
- `correctedData` (matrix): MSC‑corrected spectra.
- `referenceSpec` (row vector): reference spectrum used.

## Example
```matlab
[Xcorr, ref] = MSC(X,'Mean Spectrum',1);
```

## What the code does
Classic Multiplicative Scatter Correction using a reference spectrum. Each

## How to use it
Run `MSC.m` in MATLAB. A basic demonstration is provided in `test_MSC.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_MSC.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
