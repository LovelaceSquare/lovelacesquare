# I_SVD

Iterative SVD‑based PCA imputation. The algorithm repeatedly fills missing
values by projecting the partially observed data onto a specified number of
principal components. After each SVD reconstruction the updated matrix replaces
the missing entries and the process continues until the change between
iterations falls below a tolerance. This works well when the underlying data are
approximately low rank.

## Inputs
- `D` (matrix): data matrix possibly containing NaNs.
- `nComp` (int): number of principal components to retain.
- `maxIter` (int): maximum number of global iterations.
- `tol` (double, optional): convergence tolerance.

## Outputs
- `Dimp` (matrix): imputed data matrix.
- `T` (matrix): score matrix.
- `P` (matrix): loading matrix.
- `r2` (vector): explained variance per iteration.
- `lofc` (vector): lack of fit per iteration.

## Example
```matlab
[Dimp,T,P,r2,lof] = I_SVD(D,3,50);
```

## Author
Adrián Gómez-Sánchez

MIT License.
=======
## What the code does
I_SVD  Perform Iterative SVD-based PCA Imputation on the input data.

## How to use it
Run `I_SVD.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
