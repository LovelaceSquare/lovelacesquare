# LALS

Local Asymmetric Least Squares baseline correction with interval-specific asymmetry and smoothing parameters.

## Inputs
- `y` (vector): 1‑D signal to correct.
- `intervals` (Mx2 array or cell): index ranges where local parameters apply.
- `pVals` (vector): asymmetry parameter for each interval.
- `lambdasAsym` (vector): local smoothing penalties.
- `lambdaWhit` (scalar): smoothing outside intervals.
- `mu` (scalar): global first derivative penalty.
- `maxIter` (int): maximum IRLS iterations.
- `tol` (double): convergence tolerance.

## Outputs
- `baseline` (vector): estimated baseline.
- `weights` (vector): final IRLS weights.

## Example
```matlab
[b, w] = LALS(y,[10 20;40 50],[0.01;0.001],[1e4;2e5],100,1e3,50,1e-6);
```

## Author
Adrián Gómez-Sánchez

MIT License.
