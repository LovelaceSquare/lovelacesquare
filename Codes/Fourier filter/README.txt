# fourierFilter

Apply frequency-domain filtering to each row of a time‑series matrix using the
fast Fourier transform. The routine converts each signal to the frequency
domain, zeroes out frequencies outside the provided pass bands, and then performs
an inverse transform. This is useful for removing electrical noise or other
unwanted frequency components from densely sampled data.

## Inputs
- `data` (matrix): signals arranged as rows.
- `freqIntervals` (Nx2 array): frequency bands [low high] in Hz to keep.
- `dt` (double, optional): sampling interval in seconds.

## Outputs
- `filteredData` (matrix): signals after filtering.
- `freqVector` (vector): frequencies corresponding to the FFT bins.

## Example
```matlab
intervals = [0 15; 40 60];
[Xfilt, f] = fourierFilter(X, intervals, 0.01);
```

## What the code does
Apply frequency-domain filtering to each row of a time‑series matrix using the

## How to use it
Run `fourierFilter.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
