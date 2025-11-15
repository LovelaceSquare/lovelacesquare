function [A, B, C, lof] = PARAFAC_ALS_Lite(X, A_init, B_init, C_init, maxIter, tol)
% =========================================================================
% PARAFAC_ALS_Lite  —  Parallel Factor Analysis by Alternating Least Squares (Constrained)
% =========================================================================
%
% PURPOSE:
%   Perform Parallel Factor Analysis (PARAFAC) using the Alternating Least
%   Squares (ALS) algorithm with non-negativity constraints on all factor
%   matrices. This "Lite" version implements the classical PARAFAC-ALS core
%   algorithm with compact code and simplified visualization.
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
%   Version:      1.0 (Lite)
%   License:      MIT
%
% -------------------------------------------------------------------------
% DESCRIPTION:
%   PARAFAC_ALS_Lite decomposes a 3-way tensor X into trilinear factors:
%
%       X ≈ sum_{r=1}^R a_r ⊗ b_r ⊗ c_r
%
%   where:
%     • X (I × J × K): Data tensor (3-way array)
%     • A (I × R): Factor matrix for mode 1
%     • B (J × R): Factor matrix for mode 2
%     • C (K × R): Factor matrix for mode 3
%     • ⊗ denotes outer product
%
%   In matricized form (mode-1 unfolding):
%       X₁ ≈ A(C ⊙ B)ᵀ
%
%   where ⊙ is the Khatri-Rao product (column-wise Kronecker product).
%
%   The algorithm alternates between estimating A, B, and C using non-negative
%   least squares (FNNLS), normalizing each factor to prevent scale
%   indeterminacy. Lack of Fit (LOF) is computed per iteration, and
%   convergence plots are updated periodically.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
%   X        - (I × J × K) Data tensor [3-way array]
%   A_init   - (I × R) Initial estimate of factor matrix A (mode 1),
%              OR [] for random initialization,
%              OR scalar R (number of components) if all B_init, C_init are also []
%   B_init   - (J × R) Initial estimate of factor matrix B (mode 2) OR [] for random
%   C_init   - (K × R) Initial estimate of factor matrix C (mode 3) OR [] for random
%   maxIter  - Maximum number of ALS iterations (positive integer)
%   tol      - (Optional) Convergence tolerance for LOF change (default 1e-6)
%
%   NOTE: Any combination of initializations can be provided (0, 1, 2, or 3 factors).
%         - If at least one factor matrix is provided, R is determined from it
%         - If all are [], A_init must be a scalar specifying R
%         - All provided factor matrices must have the same number of components R
%
% -------------------------------------------------------------------------
% OUTPUT ARGUMENTS:
%   A   - (I × R) Final factor matrix for mode 1 (non-negative, absorbs scaling)
%   B   - (J × R) Final factor matrix for mode 2 (non-negative, normalized)
%   C   - (K × R) Final factor matrix for mode 3 (non-negative, normalized)
%   lof - (1 × iter) Lack of Fit (%) evolution
%
%     LOF = 100 × ||E||_F / ||X||_F
%     where ||·||_F denotes the Frobenius norm.
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW:
%   1. Initialize A, B, C from provided initializations or randomly:
%        - If A_init, B_init, or C_init provided: use those
%        - If [] : initialize randomly (non-negative)
%        - If all [] and A_init is scalar: A_init specifies R (number of components)
%   2. Repeat until convergence or maxIter:
%        a) Solve for A:   min ||X₁ - A(C ⊙ B)ᵀ||²,  A ≥ 0  (A not normalized)
%        b) Solve for B:   min ||X₂ - B(C ⊙ A)ᵀ||²,  B ≥ 0
%        c) Normalize B columns to unit norm, compensate in A
%        d) Solve for C:   min ||X₃ - C(B ⊙ A)ᵀ||²,  C ≥ 0
%        e) Normalize C columns to unit norm, compensate in A
%        f) Compute LOF and test convergence (|ΔLOF| < tol)
%        g) Update plots periodically
%
% -------------------------------------------------------------------------
% KEY FEATURES:
%   • Flexible initialization: provide any combination of A, B, C (or none)
%   • Non-negativity enforced via fast NNLS (FNNLS)
%   • B and C normalized to unit norm; A absorbs all scaling
%   • Convergence visualization for A, B, C, and LOF
%
% -------------------------------------------------------------------------
% DISCLAIMER:
%   Provided "as is" under the MIT License. The authors and Lovelace's
%   Square assume no liability. Validate results on your own datasets.
% =========================================================================


    %% ------------------------- 1) Validate inputs -------------------------
    if nargin < 5
        error('PARAFAC_ALS_Lite requires at least 5 inputs: X, A_init, B_init, C_init, and maxIter.');
    end
    if nargin < 6 || isempty(tol)
        tol = 1e-6; % keep in sync with documentation
    end

    validateattributes(X,      {'double','single'}, {'3d','real','finite','nonempty'}, mfilename, 'X', 1);
    validateattributes(maxIter,{'numeric'},         {'scalar','real','finite','positive','integer'}, mfilename, 'maxIter', 5);
    validateattributes(tol,    {'numeric'},         {'scalar','real','finite','positive'}, mfilename, 'tol', 6);

    [I, J, K] = size(X);

    % Determine which factors are provided
    A_provided = ~isempty(A_init) && ~isscalar(A_init);
    B_provided = ~isempty(B_init);
    C_provided = ~isempty(C_init);
    R_only = isscalar(A_init) && ~B_provided && ~C_provided;

    % Determine R and validate dimensions
    R = [];

    if R_only
        % A_init is a scalar specifying R
        R = A_init;
        validateattributes(R, {'numeric'}, {'scalar','real','finite','positive','integer'}, mfilename, 'R (A_init as scalar)', 2);
    else
        % Extract R from provided initializations and validate dimensions
        if A_provided
            validateattributes(A_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'A_init', 2);
            [IAinit, RA] = size(A_init);
            if IAinit ~= I
                error('A_init rows (%d) must match X mode-1 dimension (%d).', IAinit, I);
            end
            R = RA;
        end

        if B_provided
            validateattributes(B_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'B_init', 3);
            [JBinit, RB] = size(B_init);
            if JBinit ~= J
                error('B_init rows (%d) must match X mode-2 dimension (%d).', JBinit, J);
            end
            if isempty(R)
                R = RB;
            elseif R ~= RB
                error('B_init components (%d) must match other provided initializations (%d).', RB, R);
            end
        end

        if C_provided
            validateattributes(C_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'C_init', 4);
            [KCinit, RC] = size(C_init);
            if KCinit ~= K
                error('C_init rows (%d) must match X mode-3 dimension (%d).', KCinit, K);
            end
            if isempty(R)
                R = RC;
            elseif R ~= RC
                error('C_init components (%d) must match other provided initializations (%d).', RC, R);
            end
        end

        % If no initializations provided at all
        if isempty(R)
            error('At least one factor initialization must be provided, or A_init must be a scalar specifying R.');
        end
    end

    %% ------------------------- 2) Initialize ------------------------------
    lof = nan(1, maxIter);
    normX = norm(X(:));               % Frobenius norm of X

    if normX == 0
        warning('||X||_F = 0; returning zeros for A, B, C and LOF=0.');
        A = zeros(I, R, 'like', X);
        B = zeros(J, R, 'like', X);
        C = zeros(K, R, 'like', X);
        lof = 0;
        return;
    end

    % Initialize factors based on what was provided
    if A_provided
        A = max(A_init, 0);           % enforce non-negativity
    else
        A = abs(randn(I, R));         % random initialization
    end

    if B_provided
        B = max(B_init, 0);           % enforce non-negativity
    else
        B = abs(randn(J, R));         % random initialization
    end

    if C_provided
        C = max(C_init, 0);           % enforce non-negativity
    else
        C = abs(randn(K, R));         % random initialization
    end

    % Unfold tensor for each mode
    X1 = unfold(X, 1);  % I × (J*K)
    X2 = unfold(X, 2);  % J × (I*K)
    X3 = unfold(X, 3);  % K × (I*J)

    % Create figure for real-time visualization
    fig = figure('Name', 'PARAFAC-ALS Lite Convergence', ...
                 'NumberTitle', 'off', ...
                 'Color', 'w', ...
                 'Position', [100, 100, 1400, 600]);

    %% ------------------------- 3) ALS Iterations --------------------------
    fprintf('PARAFAC-ALS Lite started...\n');
    fprintf('Tensor size: %d × %d × %d\n', I, J, K);
    fprintf('Number of components: %d\n', R);
    fprintf('Max iterations: %d, Tolerance: %.2e\n\n', maxIter, tol);
    fprintf('Iter\tLOF (%%)\t\tChange\n');
    fprintf('----\t--------\t--------\n');

    for iter = 1:maxIter
        % --- 3a. Update A: Non-negative least squares ---
        %     X₁ ≈ A(C ⊙ B)ᵀ  =>  A = X₁ × pinv((C ⊙ B)ᵀ)
        %     Solve: A = fnnls((C ⊙ B), X₁')'
        Z = khatrirao(C, B);  % (J*K) × R
        A = fnnls(Z, X1')';   % I × R
        % Note: A is NOT normalized; it absorbs all scaling

        % --- 3b. Update B: Non-negative least squares ---
        %     X₂ ≈ B(C ⊙ A)ᵀ
        Z = khatrirao(C, A);  % (I*K) × R
        B = fnnls(Z, X2')';   % J × R

        % --- 3c. Normalize B (column-wise) and compensate in A ---
        for r = 1:R
            b_norm = norm(B(:, r));
            if b_norm > eps
                B(:, r) = B(:, r) / b_norm;
                A(:, r) = A(:, r) * b_norm;
            end
        end

        % --- 3d. Update C: Non-negative least squares ---
        %     X₃ ≈ C(B ⊙ A)ᵀ
        Z = khatrirao(B, A);  % (I*J) × R
        C = fnnls(Z, X3')';   % K × R

        % --- 3e. Normalize C (column-wise) and compensate in A ---
        for r = 1:R
            c_norm = norm(C(:, r));
            if c_norm > eps
                C(:, r) = C(:, r) / c_norm;
                A(:, r) = A(:, r) * c_norm;
            end
        end

        % --- 3f. Compute Lack of Fit ---
        X_reconstructed = reconstruct(A, B, C);
        E = X(:) - X_reconstructed(:);
        lof(iter) = 100 * norm(E) / normX;

        % --- 3g. Check convergence ---
        if iter > 1
            lof_change = abs(lof(iter-1) - lof(iter));
            fprintf('%4d\t%8.4f\t%8.6f\n', iter, lof(iter), lof_change);
            if lof_change < tol
                fprintf('\nConverged at iteration %d (LOF change < %.2e)\n', iter, tol);
                lof = lof(1:iter);
                break;
            end
        else
            fprintf('%4d\t%8.4f\t%8s\n', iter, lof(iter), '-');
        end

        % --- 3h. Visualization update ---
        if mod(iter, 3) == 0 || iter == 1 || iter == maxIter
            updatePlots(fig, A, B, C, lof, iter, R);
            drawnow;
        end
    end

    % Final plot update
    updatePlots(fig, A, B, C, lof, iter, R);

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

    Z = zeros(m*n, r1);
    for i = 1:r1
        Z(:, i) = kron(A(:, i), B(:, i));
    end
end

%% ========================================================================
%% HELPER FUNCTION: Reconstruct tensor from factors
%% ========================================================================
function X_hat = reconstruct(A, B, C)
    % Reconstructs tensor from PARAFAC factors
    % X ≈ sum_{r=1}^R a_r ⊗ b_r ⊗ c_r

    [I, R] = size(A);
    [J, ~] = size(B);
    [K, ~] = size(C);

    X_hat = zeros(I, J, K);
    for r = 1:R
        % Outer product: a_r ⊗ b_r ⊗ c_r
        X_hat = X_hat + outerprod3(A(:,r), B(:,r), C(:,r));
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

    T = zeros(I, J, K);
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
function updatePlots(fig, A, B, C, lof, iter, R)
    if ~ishandle(fig), return; end
    figure(fig); clf(fig);

    % Factor matrix A (mode 1)
    subplot(2,4,[1,5]);
    plot(A, 'LineWidth', 2);
    xlabel('Mode-1 Index', 'FontSize', 11);
    ylabel('Loading', 'FontSize', 11);
    title('Factor A (Mode 1)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:R, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Factor matrix B (mode 2)
    subplot(2,4,[2,6]);
    plot(B, 'LineWidth', 2);
    xlabel('Mode-2 Index', 'FontSize', 11);
    ylabel('Loading', 'FontSize', 11);
    title('Factor B (Mode 2)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:R, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Factor matrix C (mode 3)
    subplot(2,4,[3,7]);
    plot(C, 'LineWidth', 2);
    xlabel('Mode-3 Index', 'FontSize', 11);
    ylabel('Loading', 'FontSize', 11);
    title('Factor C (Mode 3)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:R, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Lack of Fit evolution
    subplot(2,4,[4,8]);
    plot(1:length(lof), lof, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Iteration', 'FontSize', 11);
    ylabel('LOF (%)', 'FontSize', 11);
    title(sprintf('Lack of Fit (Iter %d, LOF=%.4f%%)', iter, lof(end)), ...
        'FontSize', 12, 'FontWeight', 'bold');
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
