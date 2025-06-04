# normMatrix

## What the code does
normMatrix. Normalizes the input matrix using specified norms.
Run `normMatrix.m` in MATLAB. A basic demonstration is provided in `test.m`.
MATLAB R2018b or later. Add this folder to your MATLAB path.
See `test.m` for a usage example.
contact@lovelacesquare.org
MIT
  INPUTS: data (array)      : The input matrix to be normalized. normType (string) : The type of normalization: 'max'        - Maximum value normalization. 'euclidean'  - Euclidean norm normalization (row/column). 'l1'         - L1 norm normalization. 'l2'         - L2 norm normalization (entire matrix only). 'linf'       - L-infinity norm normalization. 'frobenius'  - Frobenius norm normalization (entire matrix only). dimension (string): The dimension along which to perform normalization: 'all'    - Normalization over the entire matrix. 'row'    - Row-wise normalization. 'column' - Column-wise normalization.  OUTPUT: normalizedMatrix (array): The normalized matrix based on the specified norm.  USE EXAMPLE: normalizedMatrix = normMatrix(data, 'euclidean', 'row');  DISCLAIMER: Authors and Lovelace's Square are not responsible for any issues, inaccuracies,  or data loss arising from the use of this function.

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
