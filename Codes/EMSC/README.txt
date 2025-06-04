# EMSC

Extended Multiplicative Scatter Correction. The procedure fits each spectrum to
a common reference while optionally modeling low‑order polynomial trends. By
separating multiplicative scaling from additive effects, EMSC is able to reduce
scattering variability in Raman and NIR data. The reference can be computed as
the mean/median spectrum or supplied explicitly, and the polynomial order
controls how much baseline curvature is removed.

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
