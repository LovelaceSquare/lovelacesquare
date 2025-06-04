# SavGol

Savitzky–Golay smoothing and differentiation for spectral matrices. Within each
sliding window, a polynomial is fitted and evaluated at the center point. The
result preserves peak heights and widths while suppressing noise. Derivative
orders up to the polynomial order can be computed to enhance subtle features or
estimate slopes. Window size should be odd and large enough to capture the
desired peak width.

## Inputs
- `data` (matrix): spectra arranged row-wise.
- `windowSize` (odd int): filter window length.
- `polyOrder` (int): polynomial order of the fit.
- `derivOrder` (int): derivative order (0=smoothing).
- `edgeMethod` (string, optional): edge handling mode.

## Output
- `filteredData` (matrix): smoothed or differentiated spectra.

## Example
```matlab
Xf = SavGol(X,11,3,0,'Reflection');
```

## Author
Adrián Gómez-Sánchez

MIT License.
