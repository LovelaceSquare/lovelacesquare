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

## Output
- `irf` (vector): reconstructed instrument response function.

## Example
```matlab
pen = 1e8;
irf = birfi_ls_smoothing(decay, 100, pen);
```

## Author
Adrián Gómez-Sánchez

MIT License.
