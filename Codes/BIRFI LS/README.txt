# birfi_ls

Estimates the instrument response function (IRF) of a fluorescence decay using a Hankel matrix approach with Tikhonov regularization.

## Inputs
- `decay` (vector): measured decay signal.
- `irf_size` (int): size of the IRF to recover.
- `lambda` (double): regularization strength controlling smoothness.

## Output
- `irf` (vector): estimated instrument response function.

## Example
```matlab
irf_len = 50;
lambdaVal = 1e5;
irf = birfi_ls(decay_signal, irf_len, lambdaVal);
```

## Author
Adrián Gómez-Sánchez

MIT License.
