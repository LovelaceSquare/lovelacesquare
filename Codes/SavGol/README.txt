# SavGol

Savitzky–Golay smoothing and differentiation for spectral matrices.

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
