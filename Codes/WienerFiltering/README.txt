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

## Outputs
- `filteredData` (matrix): Wiener-filtered signals.

## Example
```matlab
Xfilt = WienerFiltering(X,0.01,256,128,512);
```

## What the code does
Row‑wise Wiener filter using Welch-estimated power spectra. For each signal the

## How to use it
Run `WienerFiltering.m` in MATLAB. A basic demonstration is provided in `test.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
