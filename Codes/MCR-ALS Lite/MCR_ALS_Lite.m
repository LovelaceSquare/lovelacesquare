function [C, S, lof] = MCR_ALS_Lite(D, N, C_init, S_init, maxIter, tol, nonnegMode)
% =========================================================================
% MCR_ALS_Lite  —  Multivariate Curve Resolution by Alternating Least Squares (Constrained)
% =========================================================================
%
% PURPOSE:
%   Perform Multivariate Curve Resolution (MCR) using the Alternating Least
%   Squares (ALS) algorithm with non-negativity constraints on both the
%   concentration (C) and spectral (S) profiles. This “Lite” version
%   implements the classical MCR-ALS core algorithm with compact code and
%   simplified visualization.
%
% -------------------------------------------------------------------------
% REFERENCES:
%   • Lawton, W. H., & Sylvestre, E. A. (1971).
%     "Self modeling curve resolution." Technometrics, 13(3), 617–633.
%
%   • de Juan, A., & Tauler, R. (2021).
%     "Multivariate Curve Resolution: 50 years addressing the mixture
%      analysis problem – A review."
%      Analytica Chimica Acta, 1145, 59–78. Elsevier.
%
% -------------------------------------------------------------------------
% AUTHORSHIP & VERSION:
%   Author:       Adrián Gómez-Sánchez
%   Created:      2025-10-30
%   Reviewed by:  Lovelace's Square
%   Version:      1.2 (Lite, with optional non-negativity and improved notation)
%   License:      MIT
%
% -------------------------------------------------------------------------
% DESCRIPTION:
%   MCR_ALS_Lite decomposes a data matrix D into bilinear factors:
%
%       D = C * S + E
%
%   where:
%     • D (I × J): Data matrix (samples × variables)
%     • C (I × N): Concentration profiles of N components
%     • S (N × J): Spectral profiles of N components
%     • E (I × J): Residuals (noise, model error)
%
%   The algorithm alternates between estimating S and C using either:
%     • Non-negative least squares (FNNLS) in selected modes
%     • Unconstrained least squares in other modes
%   Spectral rows (S) are normalized to prevent scale indeterminacy.
%   Lack of Fit (LOF) is computed per iteration, and convergence plots
%   are updated periodically.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
%   D        - (I × J) Data matrix [samples × variables]
%   N        - Number of components (positive integer)
%   C_init   - (I × N) Initial estimate of concentration profiles OR [] for random
%   S_init   - (N × J) Initial estimate of spectral profiles OR [] for random
%   maxIter  - Maximum number of ALS iterations (positive integer)
%   tol      - (Optional) Convergence tolerance for LOF change (default 1e-6)
%   nonnegMode- (Optional) 1×2 logical vector [nC, nS]:
%                 nC = true  -> C constrained non-negative (concentration)
%                 nS = true  -> S constrained non-negative (spectral)
%               If scalar given, it is broadcast to both modes.
%               Default: [true true] (original fully non-negative behaviour)
%
%   NOTE: Provide C_init, S_init, both, or neither (all will be random if both []).
%         Provided matrices must have N components.
%
% -------------------------------------------------------------------------
% OUTPUT ARGUMENTS:
%   C   - (I × N) Final concentration profiles
%   S   - (N × J) Final spectral profiles (normalized)
%   lof - (1 × iter) Lack of Fit (%) evolution
%
%     LOF = 100 × ||E||_F / ||D||_F
%     where ||·||_F denotes the Frobenius norm.
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW:
%   1. Initialize C and/or S from provided initializations or randomly:
%        - If nonnegMode(1) = true: C is initialized ≥ 0
%        - If nonnegMode(2) = true: S is initialized ≥ 0
%   2. Repeat until convergence or maxIter:
%
%      If C_init provided:
%        a) Solve for S:   min ||D - C*S||²
%             - via non-negative LS (FNNLS) if nonnegMode(2) = true
%             - via unconstrained LS otherwise
%        b) Normalize each S(i,:) to unit Euclidean norm
%        c) Solve for C:   min ||D - C*S||²
%             - via non-negative LS (FNNLS) if nonnegMode(1) = true
%             - via unconstrained LS otherwise
%
%      If S_init provided (or both [] ):
%        a) Solve for C:   min ||D - C*S||²
%             - via non-negative LS (FNNLS) if nonnegMode(1) = true
%             - via unconstrained LS otherwise
%        b) Solve for S:   min ||D - C*S||²
%             - via non-negative LS (FNNLS) if nonnegMode(2) = true
%             - via unconstrained LS otherwise
%        c) Normalize each S(i,:) to unit Euclidean norm
%
%      d) Compute LOF and test convergence (|ΔLOF| < tol)
%      e) Update plots periodically
%
% -------------------------------------------------------------------------
% KEY FEATURES:
%   • Flexible initialization: C_init, S_init, both, or neither
%   • Optional non-negativity via fast NNLS (FNNLS)
%   • Spectral normalization avoids scale ambiguity
%   • Convergence visualization for C, S, and LOF
%
% -------------------------------------------------------------------------
% DISCLAIMER:
%   Provided “as is” under the MIT License. The authors and Lovelace’s
%   Square assume no liability. Validate results on your own datasets.
% =========================================================================


    %% ======================================================================
    %% STEP 1: VALIDATE INPUTS
    %% ======================================================================
    % Purpose: Ensure all inputs are valid before starting the algorithm
    % Strategy: Check dimensions, types, and logical consistency
    if nargin < 5
        error('MCR_ALS_Lite requires at least 5 inputs: D, N, C_init, S_init, and maxIter.');
    end
    if nargin < 6 || isempty(tol)
        tol = 1e-6; % keep in sync with documentation
    end
    if nargin < 7 || isempty(nonnegMode)
        nonnegMode = [true true];  % default: fully non-negative (original)
    end

    validateattributes(D,      {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'D', 1);
    validateattributes(N,      {'numeric'},         {'scalar','real','finite','positive','integer'}, mfilename, 'N', 2);
    validateattributes(maxIter,{'numeric'},         {'scalar','real','finite','positive','integer'}, mfilename, 'maxIter', 5);
    validateattributes(tol,    {'numeric'},         {'scalar','real','finite','positive'}, mfilename, 'tol', 6);

    % handle nonnegMode: scalar or 1x2
    if isscalar(nonnegMode)
        nonnegMode = logical(nonnegMode) * [1 1];
    end
    validateattributes(nonnegMode, {'logical','numeric'}, {'vector','numel',2}, mfilename, 'nonnegMode', 7);
    nonnegMode = logical(nonnegMode(:)).'; % ensure 1x2 logical row
    nonnegC = nonnegMode(1);
    nonnegS = nonnegMode(2);

    % Check which initializations are provided
    C_provided = ~isempty(C_init);
    S_provided = ~isempty(S_init);

    [I, J] = size(D);

    % Validate provided initializations match N and data dimensions
    if C_provided
        validateattributes(C_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'C_init', 3);
        [ICinit, NC] = size(C_init);
        if ICinit ~= I
            error('C_init rows (%d) must match D rows (%d).', ICinit, I);
        end
        if NC ~= N
            error('C_init columns (%d) must match N (%d).', NC, N);
        end
    end

    if S_provided
        validateattributes(S_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'S_init', 4);
        [NS, JSinit] = size(S_init);
        if JSinit ~= J
            error('S_init columns (%d) must match D columns (%d).', JSinit, J);
        end
        if NS ~= N
            error('S_init rows (%d) must match N (%d).', NS, N);
        end
    end

    %% ======================================================================
    %% STEP 2: INITIALIZE CONCENTRATION (C) AND SPECTRAL (S) PROFILES
    %% ======================================================================
    % Purpose: Set starting values for C and S before iterative optimization
    % Strategy: Use provided initializations OR generate random starting points
    %
    % KEY CONCEPT: MCR decomposes D = C×S + E
    %   C (samples × components): How much of each component in each sample
    %   S (components × variables): Pure spectral profile of each component

    % Pre-allocate LOF (Lack of Fit) tracking array
    lof = nan(1, maxIter);

    % Compute Frobenius norm of data matrix (for LOF calculation)
    % WHY: We need ||D||_F to normalize the residual error into a percentage
    normD = norm(D, 'fro');

    % Edge case: empty data matrix
    if normD == 0
        warning('||D||_F = 0; returning zeros for C, S and LOF=0.');
        C = zeros(I, N, 'like', D);
        S = zeros(N, J, 'like', D);
        lof = 0;
        return;
    end

    % --- Initialize profiles based on provided inputs ---
    % WHY: Good initialization is critical—poor start can lead to wrong solutions
    %      MCR suffers from rotational ambiguity (many C,S pairs give same D)
    if C_provided
        C = C_init;
        if nonnegC
            C = max(C, 0);            % enforce non-negativity on C_init if required
        end
        if S_provided
            S = S_init;
            if nonnegS
                S = max(S, 0);        % enforce non-negativity on S_init if required
            end
        else
            S = [];                    % S will be computed in first iteration
        end
    elseif S_provided
        S = S_init;
        if nonnegS
            S = max(S, 0);            % enforce non-negativity on S_init if required
        end
        C = [];                        % C will be computed in first iteration
    else
        % Both empty: random initialization
        C = [];
        if nonnegS
            S = abs(randn(N, J));  % Random non-negative S
        else
            S = randn(N, J);       % Random S (can be negative)
        end
    end

    % Create figure for real-time visualization
    fig = figure('Name', 'MCR-ALS Lite Convergence', ...
                 'NumberTitle', 'off', ...
                 'Color', 'w', ...
                 'Position', [100, 100, 1400, 600]);

    %% ======================================================================
    %% STEP 3: ALTERNATING LEAST SQUARES (ALS) ITERATIONS
    %% ======================================================================
    % Purpose: Iteratively refine C and S to minimize ||D - C×S||
    % Strategy: Fix one matrix, solve for the other; repeat until convergence
    %
    % KEY CONCEPT: ALS alternates between two steps:
    %   1. Fix C → solve for S (spectral profiles)
    %   2. Fix S → solve for C (concentration profiles)
    %
    % ORDER MATTERS: Which one to update first depends on initialization:
    %   - If C_init provided: solve S first (we trust C more initially)
    %   - If S_init provided (or random): solve C first
    %
    % WHY NORMALIZE S? To prevent scale ambiguity (C could get bigger while S
    %                  gets smaller, both giving the same product D = C×S)

    fprintf('MCR-ALS Lite started...\n');
    fprintf('Data size: %d samples × %d variables\n', I, J);
    fprintf('Number of components: %d\n', N);
    fprintf('Max iterations: %d, Tolerance: %.2e\n', maxIter, tol);
    fprintf('Non-negativity (C,S): [%d %d]\n\n', nonnegC, nonnegS);
    fprintf('Iter\tLOF (%%)\t\tChange\n');
    fprintf('----\t--------\t--------\n');

    for iter = 1:maxIter
        if C_provided
            % ===============================================================
            % BRANCH A: C was initialized → solve for S first
            % ===============================================================
            % Rationale: Start by estimating the pure spectra using the
            %            provided concentration information

            % ==============================================================
            % SUBSTEP 3a: Update S (Spectral profiles)
            % ==============================================================
            % Goal: Find S that minimizes ||D - C×S||²
            % Math: Non-negative least squares if nonnegS = true
            if nonnegS
                S = fnnls(C, D);  % S is (N × nVars)
            else
                S = C \ D;        % Unconstrained least squares
            end

            % ==============================================================
            % SUBSTEP 3b: Normalize spectral profiles to prevent scaling issues
            % ==============================================================
            % Strategy: Normalize each row of S (each component's spectrum)
            % WHY? Without normalization, C could grow while S shrinks
            % Compensation: Transfer S's scaling to C to maintain D ≈ C×S
            for n = 1:N
                s_norm = norm(S(n, :));
                if s_norm > eps
                    S(n, :)  = S(n, :) / s_norm;
                    C(:, n)  = C(:, n) * s_norm;
                end
            end

            % ==============================================================
            % SUBSTEP 3c: Update C (Concentration profiles)
            % ==============================================================
            % Goal: Find C that minimizes ||D - C×S||²
            % Math trick: Transpose to use fnnls efficiently
            if nonnegC
                C = fnnls(S', D')';  % C is (nSamples × N)
                C(C < 0) = 0;  % Clip tiny negative values from numerical errors
            else
                C = D / S;           % Unconstrained least squares
            end

        else % S_provided or both empty
            % ===============================================================
            % BRANCH B: S was initialized (or random) → solve for C first
            % ===============================================================
            % Rationale: Start by estimating concentrations using the
            %            provided (or random) spectral information

            % ==============================================================
            % SUBSTEP 3a: Update C (Concentration profiles)
            % ==============================================================
            % Goal: Find C that minimizes ||D - C×S||²
            if nonnegC
                C = fnnls(S', D')';  % C is (nSamples × N)
                C(C < 0) = 0;  % Clip tiny negative values from numerical errors
            else
                C = D / S;           % Unconstrained least squares
            end

            % ==============================================================
            % SUBSTEP 3b: Update S (Spectral profiles)
            % ==============================================================
            % Goal: Find S that minimizes ||D - C×S||²
            if nonnegS
                S = fnnls(C, D);  % S is (N × nVars)
            else
                S = C \ D;        % Unconstrained least squares
            end

            % ==============================================================
            % SUBSTEP 3c: Normalize spectral profiles to prevent scaling issues
            % ==============================================================
            % Strategy: Same as Branch A—normalize S and compensate in C
            for n = 1:N
                s_norm = norm(S(n, :));
                if s_norm > eps
                    S(n, :)  = S(n, :) / s_norm;
                    C(:, n)  = C(:, n) * s_norm;
                end
            end
        end

        % ==================================================================
        % SUBSTEP 3d: Compute Lack of Fit (LOF) - Quality Metric
        % ==================================================================
        % Purpose: Measure how well the bilinear model C×S fits the data D
        % Formula: LOF = 100 × ||D - C×S||_F / ||D||_F
        % Interpretation: Lower LOF = better fit (0% = perfect fit)
        E = D - C * S;  % Residual matrix
        lof(iter) = 100 * norm(E, 'fro') / normD;

        % ==================================================================
        % SUBSTEP 3e: Check for convergence
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
        % SUBSTEP 3f: Update visualization (every 3 iterations for performance)
        % ==================================================================
        % WHY? Real-time plots help monitor convergence and spot issues early
        if mod(iter, 3) == 0 || iter == 1 || iter == maxIter
            updatePlots(fig, C, S, lof, iter, N);
            drawnow;
        end
    end

    % Final plot update
    updatePlots(fig, C, S, lof, iter, N);

    fprintf('\nMCR-ALS Lite completed.\n');
    fprintf('Final LOF: %.4f%%\n\n', lof(end));
end

%% ========================================================================
%% HELPER FUNCTION: Update Visualization
%% ========================================================================
function updatePlots(fig, C, S, lof, iter, N)
    if ~ishandle(fig), return; end
    figure(fig); clf(fig);

    % Concentration profiles
    subplot(2,3,[1,4]);
    plot(C, 'LineWidth', 2);
    xlabel('Sample Index', 'FontSize', 11);
    ylabel('Concentration', 'FontSize', 11);
    title('Concentration Profiles (C)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:N, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Spectral profiles
    subplot(2,3,[2,5]);
    plot(S', 'LineWidth', 2);
    xlabel('Variable Index', 'FontSize', 11);
    ylabel('Intensity (normalized)', 'FontSize', 11);
    title('Spectral Profiles (S)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:N, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Lack of Fit evolution (log scale on y-axis)
    subplot(2,3,[3,6]);
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
