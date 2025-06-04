# cropBackground

Remove pixels from an image cube whose summed intensity is outside a given range.

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
