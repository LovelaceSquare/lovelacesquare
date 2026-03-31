# Lovelace Square Chemometric Codes

This repository collects the chemometric algorithms and utilities developed by the Lovelace Square team.
It complements the material published at [Lovelace Square](https://lovelacesquare.org), an open hub for code, data, and learning resources in chemometrics.

## Structure

Algorithms are organized under `Codes/`, usually one folder per method or tool.

Current top-level folders include:

```text
Codes/
|-- AsLS_GUI
|-- AsLs
|-- Autoscale
|-- BIRFI LS
|-- BIRFI LS smoothing
|-- Binning
|-- ColorLayer
|-- Cosmic Peak Correction
|-- Cosmic Peak Correction GUI
|-- Crop Background
|-- Crop Background GUI
|-- EMSC
|-- Fourier filter
|-- I-SVD
|-- Image registration
|-- Kernelize
|-- LASLS_CL
|-- LASLS_GUI
|-- MCR-ALS Lite
|-- MSC
|-- NIPALS
|-- Normalize Matrix
|-- O ALS
|-- PARAFAC-ALS Lite
|-- PCA ALS-QR
|-- PCA filter
|-- Pure
|-- Saturation O-ALS
|-- SavGol
|-- UnfoldImage
|-- Whittaker Smoother
|-- Whittaker Smoother with missing values
`-- Wiener Filtering
```

Many methods are provided in more than one form:
- plain MATLAB function or script implementations
- interactive GUI variants for exploratory use
- command-line or workflow-oriented variants when appropriate

Each folder typically contains:
- `*.m` MATLAB implementation files
- test or example scripts
- a local `README.md` with method-specific notes

## Featured Lite Implementations

### MCR-ALS Lite
Multivariate Curve Resolution - Alternating Least Squares.

- Bilinear decomposition: `D = C * S + E`
- Non-negativity constraints on concentration and spectral profiles
- Real-time convergence visualization
- Good starting point for learning the MCR-ALS workflow
- Documentation: [Codes/MCR-ALS Lite/README.md](Codes/MCR-ALS%20Lite/README.md)

### PARAFAC-ALS Lite
Parallel Factor Analysis - Alternating Least Squares.

- Trilinear decomposition with Khatri-Rao products
- Non-negativity constraints on all factor matrices
- Real-time convergence visualization
- Good starting point for learning PARAFAC
- Documentation: [Codes/PARAFAC-ALS Lite/README.md](Codes/PARAFAC-ALS%20Lite/README.md)

## Getting Started

1. Clone the repository:

   ```bash
   git clone https://github.com/LovelaceSquare/lovelacesquare.git
   cd lovelacesquare
   ```

2. Add the folder you want to use to the MATLAB path:

   ```matlab
   addpath(genpath('Codes/LASLS_GUI'))
   ```

3. Open the relevant local README and run the example or entry-point file for that tool.

## Documentation

Each algorithm folder includes its own `README.md` with:
- a short explanation of the method
- usage instructions
- references when available

Additional learning resources are published through [The Library](https://library.lovelacesquare.org).

## Contributing

1. Fork this repository.
2. Create a branch.
3. Commit and push your changes.
4. Open a pull request.

## License

Unless otherwise stated, code in this repository is released under the [MIT License](LICENSE).
