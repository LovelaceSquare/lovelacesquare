# whittakerSmoother

Apply the Whittaker penalized least squares smoother to each row of a matrix.
The method fits a smooth curve that balances fidelity to the data with a
penalty on the dth-order differences, effectively acting as a low-pass filter.
Larger `lambda` values yield smoother curves while the difference order `d`
controls which derivative is penalized.

## Inputs
- `inputMatrix` (matrix): data to smooth.
- `lambda` (double): smoothing parameter.
- `d` (int): order of the difference operator.

## Output
- `smoothedMatrix` (matrix): smoothed signals.

## Example
```matlab
S = whittakerSmoother(X,1e3,2);
```

## Author
Adrián Gómez-Sánchez

MIT License.
