# unfoldImage

Flatten a 3D image cube into a 2D matrix so that each pixel becomes a row. This
reshaping step is typically performed prior to multivariate analysis so that
common matrix algorithms can be applied to the spectral dimension while keeping
track of pixel order for later refolding.
## Inputs
- `Cube` (x × y × z array): image cube.
## Outputs
- `Matrix` ( (x*y) × z array ): unfolded matrix.
## Example
```matlab
M = unfoldImage(imgCube);
```

## What the code does
Flatten a 3D image cube into a 2D matrix so that each pixel becomes a row. This

## How to use it
Run `unfoldImage.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
