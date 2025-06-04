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

## Author
Adrián Gómez-Sánchez

MIT License.
