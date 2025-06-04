# simulate_spectra

Generate synthetic spectra by summing Gaussian peaks. Peak centers are chosen at
random and each peak has a Gaussian shape with width controlled by the
`width` parameter. When `type` is `'same'` all spectra share the same set of
peak positions, otherwise each spectrum is generated independently. This
function is handy for creating test data or evaluating preprocessing routines.
## Inputs
- `num_spectra` (int): number of spectra.
- `num_variables` (int): number of wavelength points.
- `n_peaks` (int): number of Gaussians per spectrum.
- `width` (double): scale factor for Gaussian width.
- `type` ("same" or "different"): whether peak positions are shared.
## Outputs
- `spectra` (matrix): generated spectra.
## Example
```matlab
spec = simulate_spectra(10,200,3,5,'same');
```

## What the code does
Generate synthetic spectra by summing Gaussian peaks. Peak centers are chosen at

## How to use it
Run `simulate_spectra.m` in MATLAB. A basic demonstration is provided in `test_simulate_spectra.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_simulate_spectra.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
