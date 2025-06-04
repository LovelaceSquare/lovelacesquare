# birfi_ls

Estimate the instrument response function (IRF) of a fluorescence decay. The
procedure builds a structured Hankel matrix from the decay trace and solves a
regularized least squares system. A Tikhonov term penalizes large second
derivatives of the IRF, yielding a smooth estimate even when the signal is
noisy. This routine is useful when an experimental IRF measurement is not
available.

## Inputs
- `decay` (vector): measured decay signal.
- `irf_size` (int): size of the IRF to recover.
- `lambda` (double): regularization strength controlling smoothness.

## Outputs
- `irf` (vector): estimated instrument response function.

## Example
```matlab
irf_len = 50;
lambdaVal = 1e5;
irf = birfi_ls(decay_signal, irf_len, lambdaVal);
```

## What the code does
Estimate the instrument response function (IRF) of a fluorescence decay. The

## How to use it
Run `birfi_ls.m` in MATLAB. A basic demonstration is provided in `test_BIRFI_LS.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_BIRFI_LS.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
