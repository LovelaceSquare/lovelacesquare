# whittakerSmootherImpute

Iteratively apply the Whittaker smoother to fill missing values (NaNs) in a
matrix. At each iteration the current estimate is smoothed, and the original
observed points are reinserted. The routine stops when successive estimates
change by less than the specified tolerance. Use this when the underlying signal
is expected to vary smoothly over the missing regions.

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

## What the code does
Iteratively apply the Whittaker smoother to fill missing values (NaNs) in a

## How to use it
Run `whittakerSmootherImpute.m` in MATLAB. A basic demonstration is provided in `test_whittakerSmootherImpute.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_whittakerSmootherImpute.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
