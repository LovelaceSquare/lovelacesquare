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

## Output
- `D_kernelized` (3D array): kernelized data [samples × kernels × time].

## Example
```matlab
Dk = kernelize(D, 5, 20);
```

## Author
Adrián Gómez-Sánchez

MIT License.
