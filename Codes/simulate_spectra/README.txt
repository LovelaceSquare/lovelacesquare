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
## Output
- `spectra` (matrix): generated spectra.
## Example
```matlab
spec = simulate_spectra(10,200,3,5,'same');
```

## Author
MIT License.
=======
## What the code does
Simulate spectra with n gaussian without noise.
Run `simulate_spectra.m` in MATLAB. A basic demonstration is provided in `test.m`.
MATLAB R2018b or later. Add this folder to your MATLAB path.
See `test.m` for a usage example.

contact@lovelacesquare.org
MIT

## Authors
Ruggero Guerrini

## License
MIT

## Version
1.0

## Date Created


## Reviewed by Lovelace's Square team
No
