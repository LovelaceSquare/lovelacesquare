# LALS

Local Asymmetric Least Squares baseline correction with interval‑dependent
parameters. The signal is divided into user supplied regions and each region can
have its own asymmetry level and smoothing strength. Outside those regions a
global Whittaker smoother is applied. This flexibility lets you aggressively
remove background around broad peaks while preserving narrow features elsewhere.

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
=======
## What the code does
LALS. Perform Local Asymmetric Least Squares (LALS) baseline correction with per-interval parameters.

## How to use it
Run `LALS.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
