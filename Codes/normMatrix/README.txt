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
 'l2'         - L2 norm normalization (entire matrix only). 'linf'       - L-infinity norm normalization. 'frobenius'  - Frobenius norm normalization (entire matrix only). dimension (string): The dimension along which to perform normalization: 'all'    - Normalization over the entire matrix. 'row'    - Row-wise normalization. 'column' - Column-wise normalization.  OUTPUT: normalizedMatrix (array): The normalized matrix based on the specified norm.  USE EXAMPLE: normalizedMatrix = normMatrix(data, 'euclidean', 'row');  DISCLAIMER: Authors and Lovelace's Square are not responsible for any issues, inaccuracies,  or data loss arising from the use of this function.

## How to use it
[Placeholder: Explain how to use the code, e.g., main function and arguments]

## Installation/setup instructions
(Please list any installation or setup steps required)

## Usage examples
(Please provide one or more examples of how to use the code)

## Contact information
(Please provide contact information for questions or support)

## Authors
Adrián Gómez-Sánchez

## License
(Refer to lovelace_square_readme.md for acceptable licenses. Please fill this manually.)

## Version
1.0

## Date Created
2024-12-14

## Reviewed by Lovelace's Square team
No
