# normMatrix

Normalize a matrix using a variety of norms. Depending on the `normType`
parameter, each row or column (or the whole matrix) is divided by its chosen
norm, such as L1, L2, infinity, or Frobenius. This allows easy comparison of
spectra or variables on a common scale.
## Inputs
- `data` (matrix): matrix to normalize.
- `normType` (string): one of `'max'`, `'euclidean'`, `'l1'`, `'l2'`, `'linf'`, `'frobenius'`.
- `dimension` (string): `'all'`, `'row'`, or `'column'`.
## Outputs
- `normalizedMatrix` (matrix): normalized matrix.
## Example
```matlab
N = normMatrix(X,'euclidean','row');
```

## What the code does
Normalize a matrix using a variety of norms. Depending on the `normType`

## How to use it
Run `normMatrix.m` in MATLAB. A basic demonstration is provided in `test_normMatrix.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_normMatrix.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
