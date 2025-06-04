# cropBackground

Remove pixels from a hyperspectral image cube whose summed intensity lies
outside a user-defined range. Each pixel is integrated across all wavelengths
and compared against the `minThresh` and `maxThresh` values. Pixels that fall
outside this interval are discarded, allowing you to crop away low-signal
background areas before further analysis.

## Inputs
- `imageCube` (rows × cols × channels): hyperspectral image.
- `minThresh` (double): lower intensity bound.
- `maxThresh` (double): upper intensity bound.

## Outputs
- `croppedMatrix` (matrix): retained pixel spectra.
- `retainedPixelsIdx` (vector): indices of kept pixels.
- `discardedPixelsIdx` (vector): indices of removed pixels.

## Example
```matlab
[cropData, keepIdx, dropIdx] = cropBackground(img,100,500);
```

## Author
Adrián Gómez-Sánchez

MIT License.
