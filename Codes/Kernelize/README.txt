# kernelize

Convolve each signal in the input matrix with a set of normalized kernels, producing a 3‑way array of kernelized data.

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
