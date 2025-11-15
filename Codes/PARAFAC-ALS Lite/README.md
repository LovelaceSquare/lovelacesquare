# PARAFAC-ALS Lite: Parallel Factor Analysis - Alternating Least Squares

**Lite** implementation of PARAFAC-ALS for decomposing three-way data into trilinear factor matrices.

## 📖 Overview

**PARAFAC-ALS Lite** is the fundamental, lightweight implementation of Parallel Factor Analysis using Alternating Least Squares. This Lite version provides the core algorithm with essential non-negativity constraints and factor normalization. It serves as the foundation for more advanced PARAFAC variants that may be added to this repository in the future.

PARAFAC (also known as CANDECOMP) is a powerful multi-way decomposition technique for analyzing three-way and higher-order data. It decomposes multi-way arrays (tensors) into a sum of rank-one tensors, providing unique solutions under mild conditions—a significant advantage over bilinear methods like PCA or MCR-ALS.

### The Trilinear Model

PARAFAC decomposes a three-way tensor according to the trilinear model:

**X ≈ ∑ᵣ aᵣ ⊗ bᵣ ⊗ cᵣ**

where:
- **X** (I × J × K): Three-way data tensor (e.g., samples × wavelengths × elution time)
- **A** (I × R): Factor matrix for mode 1 (e.g., concentration profiles)
- **B** (J × R): Factor matrix for mode 2 (e.g., spectral profiles)
- **C** (K × R): Factor matrix for mode 3 (e.g., elution profiles)
- **⊗** denotes outer product
- **R**: Number of components

In matricized form (mode-1 unfolding):

**X₁ ≈ A(C ⊙ B)ᵀ**

where **⊙** is the Khatri-Rao product (column-wise Kronecker product).

### Key Features of the Lite Version

- ✅ **Flexible initialization** with any combination of A, B, C (or none)
- ✅ **Non-negativity constraints** on all factor matrices (A, B, C)
- ✅ **Factor normalization** (B and C normalized; A absorbs scaling) to prevent scale ambiguity
- ✅ **Real-time visualization** of convergence progress
- ✅ **Lack of Fit (LOF)** monitoring at each iteration
- ✅ **Simple, lightweight implementation** with minimal dependencies
- ✅ **Foundation for future PARAFAC variants** in this repository

> **Note**: This is the **Lite** version. Future additions to this repository may include PARAFAC with additional constraints (unimodality, orthogonality), Tucker decomposition, weighted PARAFAC, and other advanced variants.

---

## Applications

PARAFAC is widely used in:

- **Fluorescence spectroscopy**: EEM (Excitation-Emission Matrix) analysis
- **Chromatography**: LC-DAD, GC-MS, HPLC with multi-detector systems
- **Process monitoring**: Batch process analysis, reaction monitoring
- **Environmental analysis**: Water quality monitoring, pollution source identification
- **Food quality**: Authenticity testing, quality control, adulteration detection
- **Pharmaceutical analysis**: Drug formulation, dissolution testing, stability studies
- **Neuroscience**: fMRI data analysis, EEG tensor decomposition

---

## 📊 Algorithm

### Alternating Least Squares (ALS) for PARAFAC

1. **Initialize** factor matrices **A**, **B**, **C** from provided initializations or randomly:
   - If factor provided: use it
   - If [] : initialize randomly (non-negative)
   - If all [] and A_init is scalar: use as number of components R
2. **Iterate** until convergence:
   - **Step 1**: Fix B and C, solve for A with non-negativity:
     `min ||X₁ - A(C ⊙ B)ᵀ||²  s.t.  A ≥ 0` (A not normalized)
   - **Step 2**: Fix A and C, solve for B with non-negativity:
     `min ||X₂ - B(C ⊙ A)ᵀ||²  s.t.  B ≥ 0`
   - **Step 3**: Normalize each column of B to unit norm, compensate in A
   - **Step 4**: Fix A and B, solve for C with non-negativity:
     `min ||X₃ - C(B ⊙ A)ᵀ||²  s.t.  C ≥ 0`
   - **Step 5**: Normalize each column of C to unit norm, compensate in A
   - **Step 6**: Calculate Lack of Fit: `LOF = 100 × ||X - X̂||_F / ||X||_F`
   - **Step 7**: Check convergence: if `|LOF(i-1) - LOF(i)| < tol`, stop
3. **Return** optimized A, B, C, and LOF history

### Tensor Unfolding and Khatri-Rao Product

**Tensor Unfolding (Matricization):**
- Mode-1 unfolding: **X₁** is I × (J·K), stacking mode-1 fibers as columns
- Mode-2 unfolding: **X₂** is J × (I·K), stacking mode-2 fibers as columns
- Mode-3 unfolding: **X₃** is K × (I·J), stacking mode-3 fibers as columns

**Khatri-Rao Product (C ⊙ B):**
Column-wise Kronecker product:
```
(C ⊙ B) = [c₁ ⊗ b₁,  c₂ ⊗ b₂,  ...,  cᵣ ⊗ bᵣ]
```

### Normalization Strategy

Only B and C factor columns are normalized by their Euclidean norms, with A absorbing all scaling:

```matlab
% Normalize factor column r in B and compensate in A
B(:,r) = B(:,r) / norm(B(:,r), 2);
A(:,r) = A(:,r) * norm(B(:,r), 2);

% Normalize factor column r in C and compensate in A
C(:,r) = C(:,r) / norm(C(:,r), 2);
A(:,r) = A(:,r) * norm(C(:,r), 2);
```

This ensures:
- Unique scaling across factor matrices
- B and C have unit intensity (normalized)
- Factor matrix A carries all the scaling information

---

## 🚀 Installation

### Prerequisites
- MATLAB R2016a or later
- No additional toolboxes required (uses custom `fnnls` implementation)

### Setup

1. Clone or download the repository
2. Add the **PARAFAC-ALS Lite** folder to your MATLAB path:

```matlab
addpath('path/to/Codes/PARAFAC-ALS Lite');
```

3. Verify installation:

```matlab
which PARAFAC_ALS_Lite
```

---

## 📝 Usage

### Example 1: Initialize with only A (modes B and C random)

```matlab
% Load your 3-way data tensor X (I × J × K)
% For example: X could be 50 samples × 40 wavelengths × 30 time points

% Initialize factor matrix A for mode-1 (e.g., 3 components)
nComponents = 3;
A_init = rand(size(X,1), nComponents);

% Run PARAFAC-ALS Lite (B_init and C_init are [], so they're random)
[A, B, C, lof] = PARAFAC_ALS_Lite(X, A_init, [], [], 100, 1e-6);

% A: Factor matrix for mode 1 (samples × components)
% B: Factor matrix for mode 2 (wavelengths × components)
% C: Factor matrix for mode 3 (time × components)
% lof: Lack of fit per iteration (%)
```

### Example 2: Initialize all three factors

```matlab
% Initialize all three factor matrices
nComponents = 3;
A_init = rand(size(X,1), nComponents);
B_init = rand(size(X,2), nComponents);
C_init = rand(size(X,3), nComponents);

% Run PARAFAC-ALS Lite with all initializations
[A, B, C, lof] = PARAFAC_ALS_Lite(X, A_init, B_init, C_init, 100, 1e-6);
```

### Example 3: All random initialization (specify number of components)

```matlab
% Specify number of components as scalar in A_init position
nComponents = 3;

% Run PARAFAC-ALS Lite with all factors randomly initialized
[A, B, C, lof] = PARAFAC_ALS_Lite(X, nComponents, [], [], 100, 1e-6);
```

### Example 4: Initialize only B and C (A random)

```matlab
% Initialize only modes 2 and 3
nComponents = 3;
B_init = rand(size(X,2), nComponents);
C_init = rand(size(X,3), nComponents);

% Run PARAFAC-ALS Lite (A_init is [], so A is random)
[A, B, C, lof] = PARAFAC_ALS_Lite(X, [], B_init, C_init, 100, 1e-6);
```

### Advanced Initialization

For better results, use informed initialization.

| Method | Description | When to Use |
|---|---|---|
| **SVD-based** | Unfold tensor, apply SVD | General purpose |
| **Random** | Random positive values | Simple problems, multiple runs |
| **DTLD** | Direct Trilinear Decomposition | When one mode is well-conditioned |
| **Pure variables** | Known regions where components are isolated | When pure variables exist |

---

### Running the Test Script

```matlab
% Run the provided test on synthetic 3-way data
test_PARAFAC_ALS_Lite
```

This will:
- Generate synthetic 3-way tensor data with 3 known components
- Apply PARAFAC-ALS Lite to recover factor matrices
- Display convergence plots in real-time
- Compare recovered vs. true factor profiles

---

## 📋 Function Reference

### `PARAFAC_ALS_Lite`

```matlab
[A, B, C, lof] = PARAFAC_ALS_Lite(X, A_init, B_init, C_init, maxIter, tol)
```

**Inputs:**
- `X` — Data tensor (I × J × K)
- `A_init` — Initial factor matrix for mode 1 (I × R), OR `[]` for random, OR scalar R if all factors are `[]`
- `B_init` — Initial factor matrix for mode 2 (J × R) OR `[]` for random
- `C_init` — Initial factor matrix for mode 3 (K × R) OR `[]` for random
- `maxIter` — Maximum iterations (default: 100)
- `tol` — Convergence tolerance for LOF change (default: 1e-6)

**Note:** Any combination of initializations can be provided (0, 1, 2, or 3 factors). All provided factors must have the same number of components R.

**Outputs:**
- `A` — Final factor matrix for mode 1 (I × R), absorbs all scaling (not normalized)
- `B` — Final factor matrix for mode 2 (J × R), columns normalized to unit norm
- `C` — Final factor matrix for mode 3 (K × R), columns normalized to unit norm
- `lof` — Lack of fit per iteration (%)

**Dependencies:**
- `fnnls.m` — Fast Non-Negative Least Squares solver (included)

### `fnnls`

```matlab
X = fnnls(A, B, tol, maxIter)
```

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

1. **Harshman, R. A. (1970)** — Foundations of the PARAFAC procedure: Models and conditions for an 'explanatory' multimodal factor analysis. *UCLA Working Papers in Phonetics*, 16, 1–84.

2. **Carroll, J. D., & Chang, J. J. (1970)** — Analysis of individual differences in multidimensional scaling via an N-way generalization of 'Eckart–Young' decomposition. *Psychometrika*, 35(3), 283–319.  
   DOI: https://doi.org/10.1007/BF02310791

3. **Bro, R. (1997)** — PARAFAC. Tutorial and applications. *Chemometrics and Intelligent Laboratory Systems*, 38(2), 149–171.  
   DOI: https://doi.org/10.1016/S0169-7439(97)00032-4

4. **Andersson, C. A., & Bro, R. (2000)** — The N-way Toolbox for MATLAB. *Chemometrics and Intelligent Laboratory Systems*, 52(1), 1–4.  
   DOI: https://doi.org/10.1016/S0169-7439(00)00071-X

5. **Bro, R., & De Jong, S. (1997)** — A fast non-negativity-constrained least squares algorithm. *Journal of Chemometrics*, 11(5), 393–401.  
   DOI:  https://doi.org/10.1002/(SICI)1099-128X(199709/10)11:5<393::AID-CEM483>3.0.CO;2-L

6. **Lawson, C. L., & Hanson, R. J. (1974)** — *Solving Least Squares Problems*. Prentice–Hall. (SIAM Classics reprint)  
   DOI: https://doi.org/10.1137/1.9781611971217

**Additional Reading**

- **Smilde, A., Bro, R., & Geladi, P. (2004)** — *Multi-way Analysis: Applications in the Chemical Sciences*. John Wiley & Sons.  
  DOI: https://doi.org/10.1002/0470012110

- **Kolda, T. G., & Bader, B. W. (2009)** — Tensor decompositions and applications. *SIAM Review*, 51(3), 455–500.  
  DOI: https://doi.org/10.1137/07070111X

- **Tomasi, G., & Bro, R. (2006)** — A comparison of algorithms for fitting the PARAFAC model. *Computational Statistics & Data Analysis*, 50(7), 1700–1734.  
  DOI: https://doi.org/10.1016/j.csda.2004.11.013


---

## 📄 License

Released under the **MIT License**.

---

## 👤 Authors

- **Adrián Gómez-Sánchez**
- **Date**: November 1, 2025
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

PARAFAC • Parallel Factor Analysis • CANDECOMP • Tensor decomposition • Trilinear decomposition • Alternating Least Squares • Multi-way analysis • Chemometrics • Non-negative least squares • Three-way data • MATLAB
