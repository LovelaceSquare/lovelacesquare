# WienerFiltering

Row‑wise Wiener filter using Welch-estimated power spectra. For each signal the
power spectral density of both the noisy data and the estimated noise are
computed with Welch's method. A Wiener gain is formed to attenuate frequency
regions dominated by noise, leading to an improved signal-to-noise ratio when
the noise can be considered additive and stationary.

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
