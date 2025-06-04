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
=======
## What the code does
CROPBACKGROUND Crop background pixels from a 3D image cube by intensity range.

## How to use it
Run `cropBackground.m` in MATLAB. A basic demonstration is provided in `test.m`.

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
