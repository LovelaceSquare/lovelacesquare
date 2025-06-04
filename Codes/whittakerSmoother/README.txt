# whittakerSmoother

Apply the Whittaker penalized least squares smoother to each row of a matrix.

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
