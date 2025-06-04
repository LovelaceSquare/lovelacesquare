# whittakerSmootherImpute

Iteratively apply the Whittaker smoother to fill missing values (NaNs) in a matrix.

## Inputs
- `inputMatrix` (matrix): data containing NaNs.
- `lambda` (double): smoothing parameter.
- `d` (int): difference order.
- `maxIter` (int): maximum iterations.
- `tol` (double): convergence tolerance.

## Outputs
- `smoothedMatrix` (matrix): smoothed result.
- `imputedMatrix` (matrix): data with NaNs replaced.

## Example
```matlab
[S, I] = whittakerSmootherImpute(X,10,2,5,1e-4);
```

## Author
Adrián Gómez-Sánchez

MIT License.
