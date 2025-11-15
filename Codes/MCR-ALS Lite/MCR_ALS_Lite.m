function [C, S, lof] = MCR_ALS_Lite(D, C_init, S_init, maxIter, tol)
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
%   Version:      1.0 (Lite)
%   License:      MIT
%
% -------------------------------------------------------------------------
% DESCRIPTION:
%   MCR_ALS_Lite decomposes a data matrix D into bilinear factors:
%
%       D = C * S + E
%
%   where:
%     • D (n × m): Data matrix (samples × variables)
%     • C (n × k): Concentration profiles of k components
%     • S (k × m): Spectral profiles of k components
%     • E (n × m): Residuals (noise, model error)
%
%   The algorithm alternates between estimating S and C using non-negative
%   least squares (FNNLS), normalizing each S row to prevent scale
%   indeterminacy. Lack of Fit (LOF) is computed per iteration, and
%   convergence plots are updated periodically.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS:
%   D        - (n × m) Data matrix [samples × variables]
%   C_init   - (n × k) Initial estimate of concentration profiles OR [] if using S_init
%   S_init   - (k × m) Initial estimate of spectral profiles OR [] if using C_init
%   maxIter  - Maximum number of ALS iterations (positive integer)
%   tol      - (Optional) Convergence tolerance for LOF change (default 1e-6)
%
%   NOTE: Provide either C_init OR S_init (set the other to []).
%         Exactly one must be non-empty.
%
% -------------------------------------------------------------------------
% OUTPUT ARGUMENTS:
%   C   - (n × k) Final concentration profiles (non-negative)
%   S   - (k × m) Final spectral profiles (non-negative, normalized)
%   lof - (1 × iter) Lack of Fit (%) evolution
%
%     LOF = 100 × ||E||_F / ||D||_F
%     where ||·||_F denotes the Frobenius norm.
%
% -------------------------------------------------------------------------
% ALGORITHM OVERVIEW:
%   1. Initialize with either C = C_init OR S = S_init
%   2. Repeat until convergence or maxIter:
%
%      If C_init provided:
%        a) Solve for S:   min ||D - C*S||²,  S ≥ 0
%        b) Normalize each S(i,:) to unit Euclidean norm
%        c) Solve for C:   min ||D - C*S||²,  C ≥ 0
%
%      If S_init provided:
%        a) Solve for C:   min ||D - C*S||²,  C ≥ 0
%        b) Solve for S:   min ||D - C*S||²,  S ≥ 0
%        c) Normalize each S(i,:) to unit Euclidean norm
%
%      d) Compute LOF and test convergence (|ΔLOF| < tol)
%      e) Update plots periodically
%
% -------------------------------------------------------------------------
% KEY FEATURES:
%   • Flexible initialization: C_init OR S_init
%   • Non-negativity enforced via fast NNLS (FNNLS)
%   • Spectral normalization avoids scale ambiguity
%   • Convergence visualization for C, S, and LOF
%
% -------------------------------------------------------------------------
% DISCLAIMER:
%   Provided “as is” under the MIT License. The authors and Lovelace’s
%   Square assume no liability. Validate results on your own datasets.
% =========================================================================


    %% ------------------------- 1) Validate inputs -------------------------
    if nargin < 4
        error('MCR_ALS_Lite requires at least 4 inputs: D, C_init, S_init, and maxIter.');
    end
    if nargin < 5 || isempty(tol)
        tol = 1e-6; % keep in sync with documentation
    end

    validateattributes(D,      {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'D', 1);
    validateattributes(maxIter,{'numeric'},         {'scalar','real','finite','positive','integer'}, mfilename, 'maxIter', 4);
    validateattributes(tol,    {'numeric'},         {'scalar','real','finite','positive'}, mfilename, 'tol', 5);

    % Check that exactly one initialization is provided
    C_provided = ~isempty(C_init);
    S_provided = ~isempty(S_init);

    if ~C_provided && ~S_provided
        error('Either C_init or S_init must be provided (both are empty).');
    end
    if C_provided && S_provided
        error('Provide only ONE initialization: either C_init or S_init (not both).');
    end

    [nSamples, nVars] = size(D);

    % Validate and extract dimensions based on which init is provided
    if C_provided
        validateattributes(C_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'C_init', 2);
        [nCinit, nComp] = size(C_init);
        if nCinit ~= nSamples
            error('C_init rows (%d) must match D rows (%d).', nCinit, nSamples);
        end
    else % S_provided
        validateattributes(S_init, {'double','single'}, {'2d','real','finite','nonempty'}, mfilename, 'S_init', 3);
        [nComp, nSinit] = size(S_init);
        if nSinit ~= nVars
            error('S_init columns (%d) must match D columns (%d).', nSinit, nVars);
        end
    end

    %% ------------------------- 2) Initialize ------------------------------
    lof = nan(1, maxIter);
    normD = norm(D, 'fro');           % Frobenius norm of D

    if normD == 0
        warning('||D||_F = 0; returning zeros for C, S and LOF=0.');
        if C_provided
            C = C_init;
            S = zeros(nComp, nVars, 'like', D);
        else
            S = S_init;
            C = zeros(nSamples, nComp, 'like', D);
        end
        lof = 0;
        return;
    end

    % Initialize based on which input was provided
    if C_provided
        C = max(C_init, 0);           % enforce non-negativity on C_init
        S = [];                        % S will be computed in first iteration
    else % S_provided
        S = max(S_init, 0);           % enforce non-negativity on S_init
        C = [];                        % C will be computed in first iteration
    end

    % Create figure for real-time visualization
    fig = figure('Name', 'MCR-ALS Lite Convergence', ...
                 'NumberTitle', 'off', ...
                 'Color', 'w', ...
                 'Position', [100, 100, 1400, 600]);

    %% ------------------------- 3) ALS Iterations --------------------------
    fprintf('MCR-ALS Lite started...\n');
    fprintf('Data size: %d samples x %d variables\n', nSamples, nVars);
    fprintf('Number of components: %d\n', nComp);
    fprintf('Max iterations: %d, Tolerance: %.2e\n\n', maxIter, tol);
    fprintf('Iter\tLOF (%%)\t\tChange\n');
    fprintf('----\t--------\t--------\n');

    for iter = 1:maxIter
        if C_provided
            % --- Standard order: C provided, solve for S first ---

            % --- 3a. Update S: Non-negative least squares ---
            %     A=C (nSamples×nComp), B=D (nSamples×nVars) -> S (nComp×nVars)
            S = fnnls(C, D);

            % --- 3b. Normalize spectral profiles (row-wise) ---
            for g = 1:nComp
                s_norm = norm(S(g, :));
                if s_norm > eps
                    S(g, :)  = S(g, :) / s_norm;
                    C(:, g)  = C(:, g) * s_norm;
                end
            end

            % --- 3c. Update C: Non-negative least squares ---
            %     Use multi-RHS trick: C = fnnls(S', D')'
            C = fnnls(S', D')';
            C(C < 0) = 0; % gentle clip for numerical noise

        else % S_provided
            % --- Reverse order: S provided, solve for C first ---

            % --- 3a. Update C: Non-negative least squares ---
            %     Use multi-RHS trick: C = fnnls(S', D')'
            C = fnnls(S', D')';
            C(C < 0) = 0; % gentle clip for numerical noise

            % --- 3b. Update S: Non-negative least squares ---
            %     A=C (nSamples×nComp), B=D (nSamples×nVars) -> S (nComp×nVars)
            S = fnnls(C, D);

            % --- 3c. Normalize spectral profiles (row-wise) ---
            for g = 1:nComp
                s_norm = norm(S(g, :));
                if s_norm > eps
                    S(g, :)  = S(g, :) / s_norm;
                    C(:, g)  = C(:, g) * s_norm;
                end
            end
        end

        % --- 3d. Compute Lack of Fit ---
        E = D - C * S;
        lof(iter) = 100 * norm(E, 'fro') / normD;

        % --- 3e. Check convergence ---
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

        % --- 3f. Visualization update ---
        if mod(iter, 3) == 0 || iter == 1 || iter == maxIter
            updatePlots(fig, C, S, lof, iter, nComp);
            drawnow;
        end
    end

    % Final plot update
    updatePlots(fig, C, S, lof, iter, nComp);

    fprintf('\nMCR-ALS Lite completed.\n');
    fprintf('Final LOF: %.4f%%\n\n', lof(end));
end

%% ========================================================================
%% HELPER FUNCTION: Update Visualization
%% ========================================================================
function updatePlots(fig, C, S, lof, iter, nComp)
    if ~ishandle(fig), return; end
    figure(fig); clf(fig);

    % Concentration profiles
    subplot(2,3,[1,4]);
    plot(C, 'LineWidth', 2);
    xlabel('Sample Index', 'FontSize', 11);
    ylabel('Concentration', 'FontSize', 11);
    title('Concentration Profiles (C)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:nComp, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Spectral profiles
    subplot(2,3,[2,5]);
    plot(S', 'LineWidth', 2);
    xlabel('Variable Index', 'FontSize', 11);
    ylabel('Intensity (normalized)', 'FontSize', 11);
    title('Spectral Profiles (S)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(arrayfun(@(x) sprintf('Comp %d', x), 1:nComp, 'UniformOutput', false), ...
        'Location', 'best');
    grid on;

    % Lack of Fit evolution
    subplot(2,3,[3,6]);
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
