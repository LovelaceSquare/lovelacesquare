# WienerFiltering

Row-wise Wiener filter using Welch-estimated power spectra.

## Inputs
- `data` (matrix): spectral data matrix.
- `noiseVar` (double): estimated noise variance.
- `segmentLength` (int): segment length for PSD estimation.
- `overlap` (int): overlap between segments.
- `nfft` (int): FFT length.

## Output
- `filteredData` (matrix): Wiener-filtered signals.

## Example
```matlab
Xfilt = WienerFiltering(X,0.01,256,128,512);
```

## Author
Adrián Gómez-Sánchez

MIT License.
