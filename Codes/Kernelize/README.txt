# kernelize

Convolve each signal in the input matrix with a bank of normalized kernels. Each
kernel is typically Gaussian shaped and shifted across the signal, creating a
third dimension representing kernel responses. The resulting tensor can be used
for matched filtering or as features for classification tasks. The width and
number of kernels determine the time resolution and coverage.

## Inputs
- `D` (matrix): input data (samples × timepoints).
- `num_kernels` (int): number of kernels to generate.
- `kernel_width` (int): width of each kernel in points.

## Outputs
- `D_kernelized` (3D array): kernelized data [samples × kernels × time].

## Example
```matlab
Dk = kernelize(D, 5, 20);
```

## What the code does
Convolve each signal in the input matrix with a bank of normalized kernels. Each

## How to use it
Run `kernelize.m` in MATLAB. A basic demonstration is provided in `test_kernelize.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_kernelize.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
