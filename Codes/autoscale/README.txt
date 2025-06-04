# autoscale

Center and scale a matrix to zero mean and unit variance by rows or columns.

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

## Author
Adrián Gómez-Sánchez

MIT License.
