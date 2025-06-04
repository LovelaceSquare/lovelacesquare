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
 % Suppose 'imgCube' is a 3D image array of dimensions 100x100x3: unfolded = unfoldImage(imgCube); % 'unfolded' will be a 10000x3 matrix – nice and flat, just like all our chemometric models like it.  DISCLAIMER: Authors and Lovelace's Square are not responsible for any existential crises  that might arise from the difficulty of this code.

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
models only date flat data, and Lovelace's Square is here to make your life easier.

## Reviewed by Lovelace's Square team
No
