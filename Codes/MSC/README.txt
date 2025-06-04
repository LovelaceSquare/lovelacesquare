# MSC

Classic Multiplicative Scatter Correction using a reference spectrum.

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
