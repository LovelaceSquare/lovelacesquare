# whittakerSmoother

Apply the Whittaker penalized least squares smoother to each row of a matrix.
The method fits a smooth curve that balances fidelity to the data with a
penalty on the dth-order differences, effectively acting as a low-pass filter.
Larger `lambda` values yield smoother curves while the difference order `d`
controls which derivative is penalized.

## Inputs
- `inputMatrix` (matrix): data to smooth.
- `lambda` (double): smoothing parameter.
- `d` (int): order of the difference operator.

## Outputs
- `smoothedMatrix` (matrix): smoothed signals.

## Example
```matlab
S = whittakerSmoother(X,1e3,2);
```

## What the code does
Apply the Whittaker penalized least squares smoother to each row of a matrix.

## How to use it
Run `whittakerSmoother.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
