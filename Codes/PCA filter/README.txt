# pcaFilter

Performs dimensionality reduction by projecting the data matrix onto a specified number of principal components and reconstructing it.

## Inputs
- `inputMatrix` (matrix): data to filter.
- `numComponents` (int): number of principal components to retain.

## Outputs
- `filteredMatrix` (matrix): reconstructed data using the selected components.

## Example
```matlab
[Xfilt] = pcaFilter(X,5);
```

## Author
Adrián Gómez-Sánchez

MIT License.
