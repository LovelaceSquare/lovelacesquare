# Lovelace Square's – Chemometric Codes

This repository collects the **chemometric algorithms and utilities** written by our team.
It is part of the contribution effort to [Lovelace's Square](https://lovelacesquare.org), an open and collaborative hub for code,
data, and learning resources in chemometrics.

## 📂 Structure

The repository is organized by algorithms, each under its own subfolder in `Codes/`:

```
Codes/
├── AsLs/                         # Asymmetric Least Squares baseline correction
├── BIRFI LS/                     # Baseline Iterative Reweighted Fitting
├── BIRFI LS smoothing/           # Smoothed Baseline Iterative Reweighted Fitting
├── EMSC/                         # Extended Multiplicative Scatter Correction
├── Fourier filter/               # Fourier-based spectral filtering
├── I-SVD/                        # Iterative SVD
├── Image registration/           # Image registration utilities
├── Kernelize/                    # Kernel preprocessing
├── LALS/                         # Localized ALS
├── MCR-ALS Classic/              # Multivariate Curve Resolution - Alternating Least Squares
├── MCR-ALS Lite/                 # MCR-ALS Lite: Lightweight foundational implementation
├── MSC/                          # Multiplicative Scatter Correction
├── O ALS/                        # Orthogonal ALS
├── PARAFAC-ALS Lite/             # PARAFAC-ALS Lite: Lightweight foundational implementation
├── PCA ALS-QR/                   # PCA-based ALS with QR decomposition
├── PCA filter/                   # PCA-based filtering
├── Pure/                         # Pure variable selection
├── Saturation O-ALS/             # Orthogonal ALS with saturation peak recovery
├── SavGol/                       # Savitzky–Golay filtering
├── WienerFiltering/              # Wiener filtering
├── autoscale/                    # Autoscaling methods
├── binning/                      # Spectral binning
├── cosmicpeakcorrection/         # Cosmic ray peak correction
├── cropBackground/               # Background cropping
├── normMatrix/                   # Normalization utilities
├── unfoldImage/                  # Image unfolding utilities
├── whittakerSmoother/            # Whittaker smoothing
└── whittakerSmoother_imputation/ # Whittaker smoothing with imputation
```

Each folder typically contains:
- `*.m` — MATLAB implementation(s)
- `test_*.m` — example/test script
- `README.md` — short description and usage notes

## 🌟 Featured: Lite Implementations

The **Lite** series provides lightweight, foundational implementations of core chemometric algorithms:

### MCR-ALS Lite
**Multivariate Curve Resolution - Alternating Least Squares**
- Bilinear decomposition: `D = C × S + E`
- Non-negativity constraints on concentration (C) and spectral (S) profiles
- Real-time convergence visualization
- Perfect for learning the fundamentals of MCR-ALS
- [📖 Documentation](Codes/MCR-ALS%20Lite/README.md)

### PARAFAC-ALS Lite
**Parallel Factor Analysis - Alternating Least Squares**
- Trilinear decomposition: `X ≈ ∑ aᵣ ⊗ bᵣ ⊗ cᵣ`
- Three-way tensor decomposition with Khatri-Rao products
- Non-negativity constraints on all factor matrices
- Real-time convergence visualization
- Perfect for learning the fundamentals of PARAFAC
- [📖 Documentation](Codes/PARAFAC-ALS%20Lite/README.md)

**What makes Lite implementations special:**
- 🎯 **Consistent structure** across implementations for easy learning
- 📚 **Educational focus** with detailed comments and documentation
- 🚀 **No external dependencies** (custom FNNLS solver included)
- 🔬 **Foundation for advanced variants** that may be added in the future
- 📊 **Real-time visualization** to understand algorithm behavior

---

## If you want to download all codes:

1. Clone this repository:
   ```bash
   git clone https://github.com/LovelaceSquare/lovelacesquare.git
   cd lovelacesquare/Codes
   ```
2. Open your code interpreter
3. Add the desired subfolder to your path:
   ```matlab
   addpath('Codes/AsLs');
   ```
4. Run the test script (if available):
   ```matlab
   test_AsLS
   ```

## 📖 Documentation

Each algorithm folder includes a `README.md` with:
- A short explanation of the method
- Basic usage instructions
- References to the original publication(s)

Additional learning resources will be available through [The Library](https://library.lovelacesquare.org).

To contribute:
1. Fork this repository.
2. Create a new branch: `git checkout -b feature/my-algorithm`.
3. Commit your changes and push: `git push origin feature/my-algorithm`.
4. Open a Pull Request.

## 📜 License

Unless otherwise stated, code in this repository is released under the [MIT License](LICENSE).
