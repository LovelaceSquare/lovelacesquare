# MCR-ALS Lite: Multivariate Curve Resolution - Alternating Least Squares

**Lite** implementation of MCR-ALS for resolving spectral mixtures into pure component profiles.

## 📖 Overview

**MCR-ALS Lite** is the fundamental, lightweight implementation of Multivariate Curve Resolution using Alternating Least Squares. This Lite version provides the core algorithm with essential non-negativity constraints and spectral normalization. It serves as the foundation for more advanced MCR variants that may be added to this repository in the future.

MCR-ALS is a powerful soft-modeling technique for decomposing complex multivariate data into chemically meaningful components. It resolves mixtures by iteratively estimating concentration profiles and pure spectral profiles from experimental data.

### The Bilinear Model

MCR-ALS decomposes a data matrix according to the bilinear model:

**D = C × S + E**

where:
- **D** (n × m): Experimental data matrix (e.g., spectra collected over time)
- **C** (n × k): Concentration profiles of k components across n samples
- **S** (k × m): Pure spectral profiles of k components across m variables
- **E** (n × m): Residual matrix (experimental noise and model error)

### Key Features of the Lite Version

- ✅ **Flexible initialization** with either C_init OR S_init
- ✅ **Non-negativity constraints** on both C and S
- ✅ **Spectral normalization** (unit Euclidean norm) to prevent scale ambiguity
- ✅ **Real-time visualization** of convergence progress
- ✅ **Lack of Fit (LOF)** monitoring at each iteration
- ✅ **Simple, lightweight implementation** with minimal dependencies
- ✅ **Foundation for future MCR variants** in this repository

> **Note**: This is the **Lite** version. Future additions to this repository may include MCR-ALS with additional constraints (unimodality, closure, selectivity), MCR-BANDS, weighted MCR-ALS, and other advanced variants.

---

## Applications

MCR-ALS is widely used in:

- **Spectroscopy**: UV-Vis, NIR, Raman, FTIR, NMR  
- **Chromatography**: HPLC-DAD, GC-MS, LC-MS  
- **Process monitoring**: Reaction kinetics, process analytical technology (PAT)  
- **Environmental analysis**: Mixture quantification  
- **Quality control**: Pharmaceutical, food, and chemical industries

---

## 📊 Algorithm

### Alternating Least Squares (ALS)

1. **Initialize** with either concentration profiles **C** OR spectral profiles **S**
2. **Iterate** until convergence:

   **If C_init provided:**
   - **Step 1**: Fix C, solve for S with non-negativity: `min ||D - C×S||²  s.t.  S ≥ 0`
   - **Step 2**: Normalize each row of S to unit norm, compensate in C
   - **Step 3**: Fix S, solve for C with non-negativity: `min ||D - C×S||²  s.t.  C ≥ 0`

   **If S_init provided:**
   - **Step 1**: Fix S, solve for C with non-negativity: `min ||D - C×S||²  s.t.  C ≥ 0`
   - **Step 2**: Fix C, solve for S with non-negativity: `min ||D - C×S||²  s.t.  S ≥ 0`
   - **Step 3**: Normalize each row of S to unit norm, compensate in C

   - **Step 4**: Calculate Lack of Fit: `LOF = 100 × ||D - C×S||_F / ||D||_F`
   - **Step 5**: Check convergence: if `|LOF(i-1) - LOF(i)| < tol`, stop

3. **Return** optimized C, S, and LOF history

### Normalization Strategy

Each spectral component (row of S) is normalized by its Euclidean norm:

% Normalize spectral row i and compensate in C
S(i,:) = S(i,:) / norm(S(i,:), 2);
C(:,i) = C(:,i) * norm(S(i,:), 2);

This ensures:
- Unique scaling between C and S  
- Spectral profiles have unit intensity  
- Concentration profiles carry the scaling information

---

## 🚀 Installation

### Prerequisites
- MATLAB R2016a or later  
- No additional toolboxes required (uses custom `fnnls` implementation)

### Setup

1. Clone or download the repository
2. Add the **MCR-ALS Lite** folder to your MATLAB path:

addpath('path/to/Codes/MCR-ALS Lite');

3. Verify installation:

which MCR_ALS_Lite

---

## 📝 Usage

### Basic Example - Initialize with Concentration Profiles

% Load your data matrix D (samples × variables)
% For example: D could be 50 spectra × 200 wavelengths

% Initialize concentration profiles (e.g., 3 components)
nComponents = 3;
C_init = rand(size(D,1), nComponents);  % Random initialization

% Run MCR-ALS Lite (C_init provided, S_init is [])
[C, S, lof] = MCR_ALS_Lite(D, C_init, [], 100, 1e-6);

% C: Concentration profiles (samples × components)
% S: Spectral profiles (components × variables)
% lof: Lack of fit per iteration (%)

### Basic Example - Initialize with Spectral Profiles

% Load your data matrix D (samples × variables)
% For example: D could be 50 spectra × 200 wavelengths

% Initialize spectral profiles (e.g., 3 components)
nComponents = 3;
S_init = rand(nComponents, size(D,2));  % Random initialization

% Run MCR-ALS Lite (C_init is [], S_init provided)
[C, S, lof] = MCR_ALS_Lite(D, [], S_init, 100, 1e-6);

% C: Concentration profiles (samples × components)
% S: Spectral profiles (components × variables)
% lof: Lack of fit per iteration (%)

### Advanced Initialization

For better results, use informed initialization. These methods can provide either concentration (C) or spectral (S) profiles:

Method | Description | Output | When to Use
---|---|---|---
**PUREST** | PURity-based Evolving Self-modeling Technique | S profiles | Pure variable detection
**ESI** | Essential spectra | S profiles | Spectral libraries available
**Pure windows** | Known regions where components are isolated | C or S | Domain knowledge of mixture
**Random** | Random positive values | C or S | Last resort, may need multiple runs

---


### Running the Test Script

% Run the provided test on synthetic data
test_MCR_ALS_Lite

This will:
- Generate synthetic spectral data with 3 known components  
- Apply MCR-ALS Lite to recover concentration and spectral profiles  
- Display convergence plots in real-time  
- Compare recovered vs. true profiles

---

## 📋 Function Reference

### `MCR_ALS_Lite`

[C, S, lof] = MCR_ALS_Lite(D, C_init, S_init, maxIter, tol)

**Inputs:**
- `D` — Data matrix (n × m)
- `C_init` — Initial concentration profiles (n × k) OR `[]` if using S_init
- `S_init` — Initial spectral profiles (k × m) OR `[]` if using C_init
- `maxIter` — Maximum iterations (default: 100)
- `tol` — Convergence tolerance for LOF change (default: 1e-6)

**Note:** Provide either `C_init` OR `S_init` (set the other to `[]`). Exactly one must be non-empty.

**Outputs:**
- `C` — Final concentration profiles (n × k)
- `S` — Final spectral profiles (k × m), each row normalized to unit norm
- `lof` — Lack of fit per iteration (%)

**Dependencies:**
- `fnnls.m` — Fast Non-Negative Least Squares solver (included)

### `fnnls`

X = fnnls(A, B, tol, maxIter)

Fast non-negative least squares solver for multiple right-hand sides.  
Solves `min_X ||A*X - B||_F^2` subject to `X ≥ 0`, one column of `B` at a time.

**Inputs:**
- `A` — Design matrix (n × p)  
- `B` — Right-hand sides (n × q)  
- `tol` — Stationarity/zero tolerance (default: `1e-12 * ||A||_F`)  
- `maxIter` — Max active-set expansions per RHS (default: `5*p`)

**Outputs:**
- `X` — Solution matrix (p × q) with non-negative entries

**Notes:**
- Based on the Lawson–Hanson active-set NNLS, with Bro–De Jong acceleration (reuse `A'*A` and `A'*B` across RHS).  
- Falls back to pseudoinverse when a passive subset is ill-conditioned.

---


## 🔬 References

1. **Tauler, R. (1995)** — Multivariate curve resolution applied to second order data. *Chemometrics and Intelligent Laboratory Systems*, 30(1), 133–146.  
   DOI: https://doi.org/10.1016/0169-7439(95)00047-X

2. **de Juan, A., & Tauler, R. (2021)** — Multivariate Curve Resolution: 50 years addressing the mixture analysis problem – A review. *Analytica Chimica Acta*, 1145, 59–78.  
   DOI: https://doi.org/10.1016/j.aca.2020.10.051

3. **Jaumot, J., de Juan, A., & Tauler, R. (2015)** — MCR-ALS GUI 2.0: New features and applications. *Chemometrics and Intelligent Laboratory Systems*, 140, 1–12.  
   DOI: https://doi.org/10.1016/j.chemolab.2014.10.003

4. **Lawton, W. H., & Sylvestre, E. A. (1971)** — Self modeling curve resolution. *Technometrics*, 13(3), 617–633.  
   DOI: https://doi.org/10.1080/00401706.1971.10488823

5. **Bro, R., & De Jong, S. (1997)** — A fast non-negativity-constrained least squares algorithm. *Journal of Chemometrics*, 11(5), 393–401.  
   DOI: https://doi.org/10.1002/%28SICI%291099-128X%28199709/10%2911%3A5%3C393%3A%3AAID-CEM483%3E3.0.CO%3B2-L

6. **Lawson, C. L., & Hanson, R. J. (1974)** — *Solving Least Squares Problems*. Prentice–Hall. (SIAM Classics reprint)  
   DOI: https://doi.org/10.1137/1.9781611971217

**Additional Reading**

- **Tauler, R., Smilde, A., & Kowalski, B. (1995)** — Selectivity, local rank, three-way data analysis and ambiguity in multivariate curve resolution. *Journal of Chemometrics*, 9(1), 31–58.  
  DOI: https://doi.org/10.1002/cem.1180090105

- **Golshan, A., Abdollahi, H., & Maeder, M. (2011)** — Resolution of rotational ambiguity for three-component systems. *Analytical Chemistry*, 83(3), 836–841.  
  DOI: https://doi.org/10.1021/ac102429q


---


## 📄 License

Released under the **MIT License**.

---

## 👤 Authors

- **Adrián Gómez-Sánchez**  
- **Date**: October 30, 2025  
- **Reviewed by**: Lovelace's Square

---


## 📧 Contributing

Contributions are welcome! Please:
1. Fork the repository  
2. Create a feature branch  
3. Submit a pull request

For issues or questions, please open an issue in this repository.

---

## 🏷️ Keywords

Multivariate Curve Resolution • MCR-ALS • Alternating Least Squares • Spectral unmixing • Chemometrics • Non-negative least squares • Bilinear decomposition • Mixture analysis • MATLAB
