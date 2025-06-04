# OALS

Orthogonalized Alternating Least Squares for PCA. Scores and loadings are
updated in an alternating fashion, and the loadings are orthogonalized at each
step to maintain principal component properties. If a loading matrix is
provided, it serves as the starting guess for refinement. Useful for PCA on
large data sets where SVD may be impractical.

## Inputs
- `D` (matrix): data matrix.
- `iter` (int): maximum number of iterations.
- `nPC` (int): number of components.
- `P` (matrix, optional): initial loadings.

## Outputs
- `Dr` (matrix): reconstructed data matrix.
- `T` (matrix): scores.
- `P` (matrix): final loadings.
- `r2` (vector): explained variance per iteration.
- `lofc` (vector): lack of fit per iteration.

## Example
```matlab
[Dr,T,P,r2,lof] = OALS(D,100,3);
```

## Author
Adrián Gómez-Sánchez

MIT License.
