# pcaFilter

Performs dimensionality reduction by projecting the data matrix onto a selected
number of principal components and then reconstructing it from that subspace.
Components capturing only noise are discarded, so the filtered matrix retains
the dominant structures of the original data. This method can be used for
baseline removal or signal denoising.

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
=======
## What the code does
pcaFilter Performs PCA filtering on the input matrix.

## How to use it
Run `pcaFilter.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
