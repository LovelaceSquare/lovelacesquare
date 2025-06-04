# binning

Group adjacent elements of an N‑D array and sum or average them according to the specified bin sizes.

## Inputs
- `data` (array): array to bin.
- `binVector` (vector): bin size for each dimension.
- `mode` ("sum" or "mean"): how to combine elements.

## Output
- `dataBinned` (array): binned array with reduced size.

## Example
```matlab
B = binning(rand(100,200),[4 5],'sum');
```

## Author
Adrián Gómez-Sánchez

MIT License.
