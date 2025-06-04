# OALS

Orthogonalized Alternating Least Squares for PCA with optional user-supplied loadings.

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
