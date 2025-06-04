# birfi_ls_smoothing

Variant of BIRFI_LS that adds a second-derivative smoothing penalty when estimating the instrument response function.

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
