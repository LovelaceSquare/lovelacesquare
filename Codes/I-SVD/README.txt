# I_SVD

Iterative SVD-based PCA imputation. Missing values are replaced by estimates from a low-rank SVD reconstruction until convergence.

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
