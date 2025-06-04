# unfoldImage

Flatten a 3D image cube into a 2D matrix so that each pixel becomes a row. This
reshaping step is typically performed prior to multivariate analysis so that
common matrix algorithms can be applied to the spectral dimension while keeping
track of pixel order for later refolding.
## Input
- `Cube` (x × y × z array): image cube.
## Output
- `Matrix` ( (x*y) × z array ): unfolded matrix.
## Example
```matlab
M = unfoldImage(imgCube);
```
## Author
MIT License.
 % Suppose 'imgCube' is a 3D image array of dimensions 100x100x3:
 unfolded = unfoldImage(imgCube);
 % 'unfolded' will be a 10000x3 matrix – nice and flat, just like all our chemometric models like it.
 
 DISCLAIMER:
 Authors and Lovelace's Square are not responsible for any existential crises 
 that might arise from the difficulty of this code.
=======
## What the code does
UNFOLDIMAGE. Flatten a 3D image cube into a 2D matrix – because chemometrics
Run `unfoldImage.m` in MATLAB. A basic demonstration is provided in `test.m`.
MATLAB R2018b or later. Add this folder to your MATLAB path.
See `test.m` for a usage example.
contact@lovelacesquare.org
MIT. Do we need a License for this code? Haha
 INPUT:
 Cube (array): A 3D numeric array of dimensions [x, y, z] representing an image cube.
 
 OUTPUT:
 Matrix (array): A 2D numeric matrix of size (x*y) x z, with the spatial dimensions flattened.
 
 EXAMPLE:
 % Suppose 'imgCube' is a 3D image array of dimensions 100x100x3:
 unfolded = unfoldImage(imgCube);
 % 'unfolded' will be a 10000x3 matrix – nice and flat, just like all our chemometric models like it.
 
 DISCLAIMER:
 Authors and Lovelace's Square are not responsible for any existential crises 
 that might arise from the difficulty of this code.


## Authors
Adrián Gómez-Sánchez

## License
MIT

## Version
1.0

## Date Created
models only date flat data, and Lovelace's Square is here to make your life easier.

## Reviewed by Lovelace's Square team
No
