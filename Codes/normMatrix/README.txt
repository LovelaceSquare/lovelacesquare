# normMatrix

Normalize a matrix using a variety of norms. Depending on the `normType`
parameter, each row or column (or the whole matrix) is divided by its chosen
norm, such as L1, L2, infinity, or Frobenius. This allows easy comparison of
spectra or variables on a common scale.
## Inputs
- `data` (matrix): matrix to normalize.
- `normType` (string): one of `'max'`, `'euclidean'`, `'l1'`, `'l2'`, `'linf'`, `'frobenius'`.
- `dimension` (string): `'all'`, `'row'`, or `'column'`.
## Output
- `normalizedMatrix` (matrix): normalized matrix.
## Example
```matlab
N = normMatrix(X,'euclidean','row');
```
## Author
MIT License.
 'l2'         - L2 norm normalization (entire matrix only).
 'linf'       - L-infinity norm normalization.
 'frobenius'  - Frobenius norm normalization (entire matrix only).
 dimension (string): The dimension along which to perform normalization:
 'all'    - Normalization over the entire matrix.
 'row'    - Row-wise normalization.
 'column' - Column-wise normalization.
 
 OUTPUT:
 normalizedMatrix (array): The normalized matrix based on the specified norm.
 
 USE EXAMPLE:
 normalizedMatrix = normMatrix(data, 'euclidean', 'row');
 
 DISCLAIMER:
 Authors and Lovelace's Square are not responsible for any issues, inaccuracies, 
 or data loss arising from the use of this function.
=======
## What the code does
normMatrix. Normalizes the input matrix using specified norms.
Run `normMatrix.m` in MATLAB. A basic demonstration is provided in `test.m`.
MATLAB R2018b or later. Add this folder to your MATLAB path.
See `test.m` for a usage example.
contact@lovelacesquare.org
MIT
 
 INPUTS:
 data (array)      : The input matrix to be normalized.
 normType (string) : The type of normalization:
 'max'        - Maximum value normalization.
 'euclidean'  - Euclidean norm normalization (row/column).
 'l1'         - L1 norm normalization.
 'l2'         - L2 norm normalization (entire matrix only).
 'linf'       - L-infinity norm normalization.
 'frobenius'  - Frobenius norm normalization (entire matrix only).
 dimension (string): The dimension along which to perform normalization:
 'all'    - Normalization over the entire matrix.
 'row'    - Row-wise normalization.
 'column' - Column-wise normalization.
 
 OUTPUT:
 normalizedMatrix (array): The normalized matrix based on the specified norm.
 
 USE EXAMPLE:
 normalizedMatrix = normMatrix(data, 'euclidean', 'row');
 
 DISCLAIMER:
 Authors and Lovelace's Square are not responsible for any issues, inaccuracies, 
 or data loss arising from the use of this function.

## License
MIT

## Version
1.0

## Date Created
2024-12-14

## Reviewed by Lovelace's Square team
No
