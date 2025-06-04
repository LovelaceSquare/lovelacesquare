# birfi_ls_smoothing

Variant of `birfi_ls` that enforces extra smoothness on the recovered IRF. It
solves the same Hankel-based least squares problem but includes a second‑order
derivative penalty directly in the objective. Adjust the `penalty` parameter to
control how strongly rapid fluctuations are suppressed. This is especially
useful when the measured decay is noisy and the true IRF is believed to be
smooth.

## Inputs
- `decay` (vector): measured decay signal.
- `irf_size` (int): desired IRF length.
- `penalty` (double): smoothing regularization parameter.

## Outputs
- `irf` (vector): reconstructed instrument response function.

## Example
```matlab
pen = 1e8;
irf = birfi_ls_smoothing(decay, 100, pen);
```

## What the code does
Variant of `birfi_ls` that enforces extra smoothness on the recovered IRF. It

## How to use it
Run `birfi_ls_smoothing.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
