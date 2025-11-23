function [A, B, C, lof] = PARAFAC_ALS_Lite(X, N, A_init, B_init, C_init, maxIter, tol, nonnegMode)
% =========================================================================
% PARAFAC_ALS_Lite  —  Parallel Factor Analysis by Alternating Least Squares (Constrained/Unconstrained)
% =========================================================================
%
% PURPOSE:
%   Perform Parallel Factor Analysis (PARAFAC) using the Alternating Least
%   Squares (ALS) algorithm with optional non-negativity constraints on the
%   factor matrices. This "Lite" version implements the classical PARAFAC-ALS
%   core algorithm with compact code and simplified visualization.
%
% -------------------------------------------------------------------------
% REFERENCES:
%   • Harshman, R. A. (1970).
%     "Foundations of the PARAFAC procedure: Models and conditions for an
%      'explanatory' multimodal factor analysis." UCLA Working Papers in Phonetics.
%
%   • Bro, R. (1997).
%     "PARAFAC. Tutorial and applications." Chemometrics and Intelligent
%     Laboratory Systems, 38(2), 149–171.
%
% -------------------------------------------------------------------------
% AUTHORSHIP & VERSION:
%   Author:       Adrián Gómez-Sánchez
%   Created:      2025-11-01
%   Reviewed by:  Lovelace's Square
%   Version:      1.3 (Lite, with per-mode non-negativity and enhanced docs)
%   License:      MIT
%
% -------------------------------------------------------------------------
% DESCRIPTION:
%   PARAFAC_ALS_Lite decomposes a 3-way tensor X into trilinear factors:
%
%       X ≈ sum_{n=1}^N a_n ⊗ b_n ⊗ c_n
%
%   where:
%     • X (I × J × K): Data tensor (3-way array)
%     • A (I × N): Factor matrix for mode 1
%     • B (J × N): Factor matrix for mode 2
%     • C (K × N): Factor matrix for mode 3
%     • ⊗ denotes outer product
%
%   In matricized form (mode-1 unfolding):
%       X₁ ≈ A(C ⊙ B)ᵀ
%
%   where ⊙ is the Khatri-Rao product (column-wise Kronecker product).
%
%   The algorithm alternates between estimating A, B, and C using either:
%     • Non-negative least squares (FNNLS) in selected modes
%     • Unconstrained least squares in other modes
%
%   B and C are normalized each iteration to prevent scale indeterminacy,
%   with A absorbing the scaling. Lack of Fit (LOF) is computed per
%   iteration, and convergence plots are updated periodically.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
%   X         - (I × J × K) Data tensor [3-way array]
%   N         - Number of components (positive integer)
%   A_init    - (I × N) Initial estimate of factor matrix A (mode 1),
%               OR [] for random initialization
%   B_init    - (J × N) Initial estimate of factor matrix B (mode 2),
%               OR [] for random initialization
%   C_init    - (K × N) Initial estimate of factor matrix C (mode 3),
%               OR [] for random initialization
%   maxIter   - Maximum number of ALS iterations (positive integer)
%   tol       - (Optional) Convergence tolerance for LOF change (default 1e-6)
%   nonnegMode- (Optional) 1×3 logical vector [nA, nB, nC]:
%                 nA = true  -> A constrained non-negative (mode 1)
%                 nB = true  -> B constrained non-negative (mode 2)
%                 nC = true  -> C constrained non-negative (mode 3)
%               If scalar given, it is broadcast to all modes.
%               Default: [true true true] (original fully non-negative behaviour)
%
%   NOTE: Any combination of initializations can be provided (0, 1, 2, or 3 factors).
%         - Provided factor matrices must have N columns
%         - Pass [] for any factor to use random initialization
%
% -------------------------------------------------------------------------
% OUTPUT ARGUMENTS:
%   A   - (I × N) Final factor matrix for mode 1 (absorbs all scaling)
%   B   - (J × N) Final factor matrix for mode 2 (normalized columns)
%   C   - (K × N) Final factor matrix for mode 3 (normalized columns)
%   lof - (1 × iter) Lack of Fit (%) evolution
%
%     LOF = 100 × ||E||_F / ||X||_F
%     where ||·||_F denotes the Frobenius norm.
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW:
%   1. Initialize A, B, C from provided initializations or randomly:
%        - If A_init, B_init, or C_init provided: use those
%        - If [] : initialize randomly (with N components)
%        - If nonnegMode(m) = true for mode m: initial factor in that mode is made ≥ 0
%   2. Repeat until convergence or maxIter:
%        a) Solve for A:   min ||X₁ - A(C ⊙ B)ᵀ||
%             - via non-negative LS (FNNLS) if nonnegMode(1) = true
%             - via unconstrained LS otherwise
%        b) Solve for B:   min ||X₂ - B(C ⊙ A)ᵀ||
%        c) Normalize B columns to unit norm, compensate in A
%        d) Solve for C:   min ||X₃ - C(B ⊙ A)ᵀ||
%        e) Normalize C columns to unit norm, compensate in A
%        f) Compute LOF and test convergence (|ΔLOF| < tol)
%        g) Update plots periodically
%
% -------------------------------------------------------------------------
% KEY FEATURES:
%   • Flexible initialization: provide any combination of A, B, C (or none)
%   • Per-mode non-negativity via fast NNLS (FNNLS)
%   • B and C normalized to unit norm; A absorbs all scaling
%   • Convergence visualization for A, B, C, and LOF
%
% -------------------------------------------------------------------------
% DISCLAIMER:
%   Provided "as is" under the MIT License. The authors and Lovelace's
%   Square assume no liability. Validate results on your own datasets.
% =========================================================================


    %% ======================================================================
    %% STEP 1: VALIDATE INPUTS
    %% ======================================================================
    % Purpose: Ensure all inputs are valid before starting the algorithm
    % Strategy: Check dimensions, types, and logical consistency
    if nargin < 6
        error('PARAFAC_ALS_Lite requires at least 6 inputs: X, N, A_init, B_init, C_init, and maxIter.');
    end
    if nargin < 7 || isempty(tol)
        tol = 1e-6; % keep in sync with documentation
    end
    if nargin < 8 || isempty(nonnegMode)
        nonnegMode = [true true true];  % default: fully non-negative (original)
    end

    validateattributes(X,      {'double','single'}, {'3d','real','finite','nonempty'}, mfilename, 'X', 1);
    validateattributes(N,      {'numeric'},         {'scalar','real','finite','positive','integer'}, mfilename, 'N', 2);
    validateattributes(maxIter,{'numeric'},         {'scalar','real','finite','positive','integer'}, mfilename, 'maxIter', 6);
    validateattributes(tol,    {'numeric'},         {'scalar','real','finite','positive'}, mfilename, 'tol', 7);

    % handle nonnegMode: scalar or 1x3
    if isscalar(nonnegMode)
        nonnegMode = logical(nonnegMode) * [1 1 1];
    end
    validateattributes(nonnegMode, {'logical','numeric'}, {'vector','numel',3}, mfilename, 'nonnegMode', 8);
    nonnegMode = logical(nonnegMode(:)).'; % ensure 1x3 logical row
    nonnegA = nonnegMode(1);
    nonnegB = nonnegMode(2);
    nonnegC = nonnegMode(3);

    [I, J, K] = size(X);

    % Determine which factors are provided
    A_provided = ~isempty(A_init);
    B_provided = ~isempty(B_init);
    C_provided = ~isempty(C_init);

    % Validate provided initializations match N and tensor dimensions
    if A_provided
        validateattributes(A_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'A_init', 3);
        [IAinit, NA] = size(A_init);
        if IAinit ~= I
            error('A_init rows (%d) must match X mode-1 dimension (%d).', IAinit, I);
        end
        if NA ~= N
            error('A_init columns (%d) must match N (%d).', NA, N);
        end
    end

    if B_provided
        validateattributes(B_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'B_init', 4);
        [JBinit, NB] = size(B_init);
        if JBinit ~= J
            error('B_init rows (%d) must match X mode-2 dimension (%d).', JBinit, J);
        end
        if NB ~= N
            error('B_init columns (%d) must match N (%d).', NB, N);
        end
    end

    if C_provided
        validateattributes(C_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'C_init', 5);
        [KCinit, NC] = size(C_init);
        if KCinit ~= K
            error('C_init rows (%d) must match X mode-3 dimension (%d).', KCinit, K);
        end
        if NC ~= N
            error('C_init columns (%d) must match N (%d).', NC, N);
        end
    end

    %% ======================================================================
    %% STEP 2: INITIALIZE FACTOR MATRICES
    %% ======================================================================
    % Purpose: Set starting values for A, B, C before iterative optimization
    % Strategy: Use provided initializations OR generate random starting points

    % Pre-allocate LOF (Lack of Fit) tracking array
    lof = nan(1, maxIter);

    % Compute Frobenius norm of data tensor (for LOF calculation)
    % WHY: We need ||X||_F to normalize the residual error into a percentage
    normX = norm(X(:));

    % Edge case: empty tensor
    if normX == 0
        warning('||X||_F = 0; returning zeros for A, B, C and LOF=0.');
        A = zeros(I, N, 'like', X);
        B = zeros(J, N, 'like', X);
        C = zeros(K, N, 'like', X);
        lof = 0;
        return;
    end

    % --- Initialize factor matrices based on provided inputs ---
    % WHY: Good initialization can significantly improve convergence speed
    %      and help avoid poor local minima
    % --- Mode 1: A ---
    if A_provided
        A = A_init;
        if nonnegA
            A = max(A, 0);
        end
    else
        if nonnegA
            A = abs(randn(I, N));
        else
            A = randn(I, N);
        end
    end

    % --- Mode 2: B ---
    if B_provided
        B = B_init;
        if nonnegB
            B = max(B, 0);
        end
    else
        if nonnegB
            B = abs(randn(J, N));
        else
            B = randn(J, N);
        end
    end

    % --- Mode 3: C ---
    if C_provided
        C = C_init;
        if nonnegC
            C = max(C, 0);
        end
    else
        if nonnegC
            C = abs(randn(K, N));
        else
            C = randn(K, N);
        end
    end

    % --- Unfold (matricize) tensor along each mode ---
    % WHY: ALS solves matrix least-squares problems, so we need matrix forms
    %      Unfolding "flattens" the 3-way tensor into a 2D matrix for each mode
    X1 = unfold(X, 1);  % Mode-1 unfolding: I × (J*K)
    X2 = unfold(X, 2);  % Mode-2 unfolding: J × (I*K)
    X3 = unfold(X, 3);  % Mode-3 unfolding: K × (I*J)

    % --- Create figure for real-time visualization ---
    % WHY: Visual feedback helps monitor convergence and detect issues
    fig = figure('Name', 'PARAFAC-ALS Lite Convergence', ...
                 'NumberTitle', 'off', ...
                 'Color', 'w', ...
                 'Position', [100, 100, 1400, 600]);

    %% ======================================================================
    %% STEP 3: ALTERNATING LEAST SQUARES (ALS) ITERATIONS
    %% ======================================================================
    % Purpose: Iteratively optimize A, B, C to minimize ||X - X_reconstructed||
    % Strategy: Fix two factors, solve for the third; repeat until convergence
    %
    % KEY CONCEPT: ALS alternates between three steps:
    %   1. Fix B and C → solve for A
    %   2. Fix A and C → solve for B (then normalize)
    %   3. Fix A and B → solve for C (then normalize)
    %
    % WHY NORMALIZE? To prevent scale ambiguity (A could get smaller while B
    %                gets bigger, both representing the same decomposition)

    fprintf('PARAFAC-ALS Lite started...\n');
    fprintf('Tensor size: %d × %d × %d\n', I, J, K);
    fprintf('Number of components: %d\n', N);
    fprintf('Max iterations: %d, Tolerance: %.2e\n', maxIter, tol);
    fprintf('Non-negativity (A,B,C): [%d %d %d]\n\n', nonnegA, nonnegB, nonnegC);
    fprintf('Iter\tLOF (%%)\t\tChange\n');
    fprintf('----\t--------\t--------\n');

    for iter = 1:maxIter
        % ==================================================================
        % SUBSTEP 3a: Update factor matrix A (Mode 1)
        % ==================================================================
        % Goal: Solve for A while keeping B and C fixed
        % Math: X₁ ≈ A(C ⊙ B)ᵀ, where ⊙ is Khatri-Rao product
        % WHY KHATRI-RAO? It's the key operation that relates the matricized
        %                 tensor to the factor matrices in PARAFAC

        % Compute Khatri-Rao product: column-wise Kronecker product
        Z = khatrirao(C, B);  % (J*K) × N

        if nonnegA
            % Non-negative least squares: enforce A ≥ 0
            % WHY? Physical/chemical constraints (e.g., concentrations can't be negative)
            A = fnnls(Z, X1')';   % I × N
        else
            % Unconstrained least squares: A = X1 / Z'
            A = X1 / Z';
        end
        % NOTE: A is NOT normalized—it will absorb all scaling from B and C

        % ==================================================================
        % SUBSTEP 3b: Update factor matrix B (Mode 2)
        % ==================================================================
        % Goal: Solve for B while keeping A and C fixed
        % Math: X₂ ≈ B(C ⊙ A)ᵀ
        Z = khatrirao(C, A);  % (I*K) × N
        if nonnegB
            B = fnnls(Z, X2')';   % J × N
        else
            B = X2 / Z';
        end

        % ==================================================================
        % SUBSTEP 3c: Normalize B to prevent scale indeterminacy
        % ==================================================================
        % Strategy: Normalize each column of B to unit norm
        % WHY? Without normalization, A could shrink while B grows (or vice versa)
        %      giving the same reconstruction but unstable factor matrices
        % Compensation: Transfer B's scaling to A to maintain X ≈ A(C⊙B)ᵀ
        for n = 1:N
            b_norm = norm(B(:, n));
            if b_norm > eps
                B(:, n) = B(:, n) / b_norm;  % Normalize B
                A(:, n) = A(:, n) * b_norm;  % Compensate in A
            end
        end

        % ==================================================================
        % SUBSTEP 3d: Update factor matrix C (Mode 3)
        % ==================================================================
        % Goal: Solve for C while keeping A and B fixed
        % Math: X₃ ≈ C(B ⊙ A)ᵀ
        Z = khatrirao(B, A);  % (I*J) × N
        if nonnegC
            C = fnnls(Z, X3')';   % K × N
        else
            C = X3 / Z';
        end

        % ==================================================================
        % SUBSTEP 3e: Normalize C to prevent scale indeterminacy
        % ==================================================================
        % Strategy: Same as for B—normalize C and compensate in A
        for n = 1:N
            c_norm = norm(C(:, n));
            if c_norm > eps
                C(:, n) = C(:, n) / c_norm;  % Normalize C
                A(:, n) = A(:, n) * c_norm;  % Compensate in A
            end
        end

        % ==================================================================
        % SUBSTEP 3f: Compute Lack of Fit (LOF) - Quality Metric
        % ==================================================================
        % Purpose: Measure how well the model fits the data
        % Formula: LOF = 100 × ||X - X_reconstructed||_F / ||X||_F
        % Interpretation: Lower LOF = better fit (0% = perfect fit)
        X_reconstructed = reconstruct(A, B, C);
        E = X(:) - X_reconstructed(:);
        lof(iter) = 100 * norm(E) / normX;

        % ==================================================================
        % SUBSTEP 3g: Check for convergence
        % ==================================================================
        % Convergence criterion: LOF change between iterations < tolerance
        % WHY? If LOF isn't improving significantly, we've reached the optimum
        if iter > 1
            lof_change = abs(lof(iter-1) - lof(iter));
            fprintf('%4d\t%8.4f\t%8.6f\n', iter, lof(iter), lof_change);
            if lof_change < tol
                fprintf('\nConverged at iteration %d (LOF change < %.2e)\n', iter, tol);
                lof = lof(1:iter);  % Trim unused LOF entries
                break;
            end
        else
            fprintf('%4d\t%8.4f\t%8s\n', iter, lof(iter), '-');
        end

        % ==================================================================
        % SUBSTEP 3h: Update visualization (every 3 iterations for performance)
        % ==================================================================
        % WHY? Real-time plots help monitor convergence and spot issues early
        if mod(iter, 3) == 0 || iter == 1 || iter == maxIter
            updatePlots(fig, A, B, C, lof, iter, N);
            drawnow;
        end
    end

    % Final plot update
    updatePlots(fig, A, B, C, lof, min(iter, numel(lof)), N);

    fprintf('\nPARAFAC-ALS Lite completed.\n');
    fprintf('Final LOF: %.4f%%\n\n', lof(end));
end

%% ========================================================================
%% HELPER FUNCTION: Unfold tensor into matrix
%% ========================================================================
function Xn = unfold(X, mode)
    % Unfolds 3-way tensor X along specified mode
    % mode=1: I × (J*K)
    % mode=2: J × (I*K)
    % mode=3: K × (I*J)
    [I, J, K] = size(X);

    if mode == 1
        Xn = reshape(X, I, J*K);
    elseif mode == 2
        Xn = reshape(permute(X, [2 1 3]), J, I*K);
    elseif mode == 3
        Xn = reshape(permute(X, [3 1 2]), K, I*J);
    else
        error('Mode must be 1, 2, or 3 for 3-way tensor.');
    end
end

%% ========================================================================
%% HELPER FUNCTION: Khatri-Rao product (column-wise Kronecker product)
%% ========================================================================
function Z = khatrirao(A, B)
    % Computes Khatri-Rao product: Z = A ⊙ B
    % If A is m×r and B is n×r, then Z is (m*n)×r
    % Each column of Z is kron(A(:,i), B(:,i))

    [m, r1] = size(A);
    [n, r2] = size(B);

    if r1 ~= r2
        error('A and B must have the same number of columns.');
    end

    Z = zeros(m*n, r1, 'like', A);
    for i = 1:r1
        Z(:, i) = kron(A(:, i), B(:, i));
    end
end

%% ========================================================================
%% HELPER FUNCTION: Reconstruct tensor from factors
%% ========================================================================
function X_hat = reconstruct(A, B, C)
    % Reconstructs tensor from PARAFAC factors
    % X ≈ sum_{n=1}^N a_n ⊗ b_n ⊗ c_n

    [I, N] = size(A);
    [J, ~] = size(B);
    [K, ~] = size(C);

    X_hat = zeros(I, J, K, 'like', A);
    for n = 1:N
        % Outer product: a_n ⊗ b_n ⊗ c_n
        X_hat = X_hat + outerprod3(A(:,n), B(:,n), C(:,n));
    end
end

%% ========================================================================
%% HELPER FUNCTION: 3-way outer product
%% ========================================================================
function T = outerprod3(a, b, c)
    % Computes 3-way outer product: a ⊗ b ⊗ c
    % Result is a tensor of size length(a) × length(b) × length(c)

    I = length(a);
    J = length(b);
    K = length(c);

    T = zeros(I, J, K, 'like', a);
    for i = 1:I
        for j = 1:J
            for k = 1:K
                T(i,j,k) = a(i) * b(j) * c(k);
            end
        end
    end
end

%% ========================================================================
%% HELPER FUNCTION: Update Visualization
%% ========================================================================
function updatePlots(fig, A, B, C, lof, iter, N)
    if ~ishandle(fig), return; end
    figure(fig); clf(fig);

    % Factor matrix A (mode 1)
    subplot(2,4,[1,5]);
    plot(A, 'LineWidth', 2);
    xlabel('Mode-1 Index', 'FontSize', 11);
    ylabel('Loading', 'FontSize', 11);
    title('Factor A (Mode 1)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:N, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Factor matrix B (mode 2)
    subplot(2,4,[2,6]);
    plot(B, 'LineWidth', 2);
    xlabel('Mode-2 Index', 'FontSize', 11);
    ylabel('Loading', 'FontSize', 11);
    title('Factor B (Mode 2)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:N, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Factor matrix C (mode 3)
    subplot(2,4,[3,7]);
    plot(C, 'LineWidth', 2);
    xlabel('Mode-3 Index', 'FontSize', 11);
    ylabel('Loading', 'FontSize', 11);
    title('Factor C (Mode 3)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:N, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Lack of Fit evolution (log scale on y-axis)
    subplot(2,4,[4,8]);
    % Only plot non-NaN values
    valid_idx = ~isnan(lof);
    if any(valid_idx)
        semilogy(find(valid_idx), lof(valid_idx), 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    end
    xlabel('Iteration', 'FontSize', 11);
    ylabel('LOF (%) (log scale)', 'FontSize', 11);
    % Handle NaN in title
    if iter > 0 && iter <= length(lof) && ~isnan(lof(iter))
        title(sprintf('Lack of Fit (Iter %d, LOF=%.4f%%)', iter, lof(iter)), ...
            'FontSize', 12, 'FontWeight', 'bold');
    else
        title('Lack of Fit', 'FontSize', 12, 'FontWeight', 'bold');
    end
    grid on;
    xlim([1, max(10, length(lof))]);
end

%% ===================== HELPER FUNCTION ==================================
function X = fnnls(A, B, tol, maxIter)
% FNNLS  —  Fast Non-Negative Least Squares for multiple RHS
%   X = fnnls(A, B) solves, for each column b of B,
%       min_x ||A*x - b||_2^2  subject to  x >= 0
%   and returns X whose columns are the solutions.
%
%   X = fnnls(A, B, tol, maxIter) lets you set a stationarity tolerance
%   and a maximum number of active-set expansions per column.
%
% Inputs
%   A       (n x p) design matrix
%   B       (n x q) right-hand sides (one NNLS solve per column)
%   tol     (optional) stationarity/zero tolerance (default: 1e-12 * ||A||_F)
%   maxIter (optional) max active-set iterations per column (default: 5*p)
%
% Output
%   X       (p x q) non-negative solutions
%
% Notes
%   - Based on the classic Lawson–Hanson active-set NNLS method, with the
%     Bro–De Jong acceleration that reuses A'*A and A'*B across columns
%   - Falls back to pseudoinverse if a passive subset becomes ill-conditioned.
%
% REFERENCES (core methods & implementations)
%   • Lawson, C. L., & Hanson, R. J. (1974).
%     Solving Least Squares Problems. Prentice–Hall, Englewood Cliffs, NJ.
%     (Original NNLS active-set algorithm.)
%
%   • Bro, R., & De Jong, S. (1997).
%     "A fast non-negativity-constrained least squares algorithm."
%     Journal of Chemometrics, 11(5), 393–401.
%

    % ---- argument checks ----
    if nargin < 2
        error('fnnls requires A and B.');
    end
    [nA, p] = size(A);
    [nB, q] = size(B);
    if nA ~= nB
        error('Row mismatch: size(A,1)=%d must equal size(B,1)=%d.', nA, nB);
    end

    if nargin < 3 || isempty(tol)
        tol = 1e-12 * norm(A, 'fro');
    end
    if nargin < 4 || isempty(maxIter)
        maxIter = 5 * p;
    end

    % ---- precomputations reused for all columns ----
    G = A' * A;        % (p x p)
    H = A' * B;        % (p x q)

    X = zeros(p, q, 'like', B);
    % For numerical safety, force symmetry on G
    G = (G + G.') * 0.5;

    for j = 1:q
        % Active-set bookkeeping for column j
        passive = false(p, 1);  % true => in the model
        x = zeros(p, 1, 'like', B);
        w = H(:, j) - G * x;    % reduced gradient (Lagrange multipliers)

        it = 0;
        % Outer loop: add variables with positive gradient to passive set
        while any(~passive & w > tol) && it < maxIter
            it = it + 1;

            % Choose the most positive candidate to enter the passive set
            wMasked = w;
            wMasked(passive) = -Inf;
            [~, t] = max(wMasked);
            passive(t) = true;

            % Solve restricted LS over the passive set
            z = zeros(p, 1, 'like', B);
            z(passive) = safeSolveSubset(G, H(:, j), passive);

            % If any passive coeffs are non-positive, step back to boundary
            while any(z(passive) <= tol)
                negIdx = passive & (z <= tol);
                denom = (x(negIdx) - z(negIdx));
                denom(abs(denom) < eps) = eps;
                alpha = min(x(negIdx) ./ denom);
                alpha = max(0, min(1, alpha));

                x = x + alpha * (z - x);
                zeroIdx = passive & (abs(x) <= tol);
                x(zeroIdx) = 0;
                passive(zeroIdx) = false;

                z = zeros(p, 1, 'like', B);
                if any(passive)
                    z(passive) = safeSolveSubset(G, H(:, j), passive);
                end
            end

            % Accept new solution and recompute multipliers
            x = z;
            w = H(:, j) - G * x;
        end

        X(:, j) = max(x, 0); % clip tiny negatives
    end
end

function zP = safeSolveSubset(G, h, passiveMask)
    GP = G(passiveMask, passiveMask);
    hP = h(passiveMask);
    if isempty(GP)
        zP = zeros(0,1);
        return
    end
    % Try backslash; if ill-conditioned, fall back to pinv
    if rcond(GP) > 1e-12
        zP = GP \ hP;
    else
        zP = pinv(GP) * hP;
    end
end
