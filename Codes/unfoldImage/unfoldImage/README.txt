# unfoldImage: Flatten 3D Image Cubes for Chemometric Use (MATLAB)

Because flattening a 3D image cube by hand is so last century.  
Turns `[x × y × z]` chaos into `(x*y) × z` neatness—your chemometric models will thank you.

**Reference**:  
Gómez-Sánchez, Adrián. (2024). *unfoldImage* function.  
Lovelace’s Square. https://lovelacesquare.org/

---

## Overview

Got a fancy hyperspectral cube that’s mocking you every time you try to reshape it?  
`unfoldImage` is here to rescue you from the tyranny of dimensions.  
It smashes your `[x, y]` spatial grid into one long list of pixels, while keeping the spectral channels intact. Voilà: a tidy table of pixel spectra, minus the sweat and tears.

Why you need it:

- **Tedium exterminator**: No more typing `reshape(Cube, x*y, z)` and praying it works  
- **No more squinting** at the workspace trying to find the dimensions  
- **RSI prevention**: Protect your precious fingers from endless resizing  
- **Instant gratification**: Feed it directly to PCA, PLS, or whatever monster you’re wrestling with

---

## Inputs

- `imgCube` (3D array): A numeric array of size `[x × y × z]`, where `x` and `y` are spatial dimensions and `z` is the spectral/channel dimension

## Outputs

- `flat` (2D matrix): A reshaped array of size `[x*y × z]`, suitable for model input

---

## Usage Example

Paste into MATLAB:

```matlab
% Imagine you have an image cube called imgCube (size 100×100×50)
imgCube = rand(100,100,50); 

% Flatten that beast
flat = unfoldImage(imgCube);  % now 10000×50
```

---

## Installation

### Compatibility

- Works with any MATLAB version that supports `reshape` (MATLAB 4.0 or later)

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: December 16, 2024  
- **Reviewed by**: Lovelace’s Square (this one was a hard one)

---

## Changelog

- **v1.0 (2024-12-16)**:  
  Initial release

---

## Keywords

- flatten  
- image  
- cube  
- chemometrics  
- reshape  
- hyperspectral  
- MATLAB  
- pixel spectra  
- unfold
