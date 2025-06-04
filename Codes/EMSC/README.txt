# EMSC

Extended Multiplicative Scatter Correction.
Removes additive and multiplicative effects by aligning each spectrum to a reference and optional polynomial baseline.

## Inputs
- `data` (matrix): spectral data with one spectrum per row.
- `refType` ("Mean" or "Median"): how to compute the reference spectrum.
- `polyOrder` (int): polynomial order for baseline modeling (0 for none).

## Outputs
- `correctedData` (matrix): EMSC-corrected data.
- `referenceSpec` (row vector): reference spectrum used.

## Example
```matlab
[XCorr, ref] = EMSC(X, 'Mean', 1);
```

## Author
Adrián Gómez-Sánchez

MIT License.
