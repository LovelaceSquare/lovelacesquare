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

## Author
Adrián Gómez-Sánchez

MIT License.
