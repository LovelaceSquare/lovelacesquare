# autoscale

Center and scale a matrix to zero mean and unit variance by rows or columns. The
operation standardizes each row or column separately so that subsequent analysis
is not biased by absolute intensities. The accompanying structure returns the
mean and standard deviation for possible inverse transformation.

## Inputs
- `inputMatrix` (matrix): data to scale.
- `direction` ("column" or "row", default "column"): orientation of scaling.

## Outputs
- `scaledMatrix` (matrix): autoscaled data.
- `scalingParams` (struct): means and standard deviations used.

## Example
```matlab
[Xs, params] = autoscale(X,'column');
```

## What the code does
Center and scale a matrix to zero mean and unit variance by rows or columns. The

## How to use it
Run `autoscale.m` in MATLAB. A basic demonstration is provided in `test_autoscale.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_autoscale.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
