function [Dr, T, P, r2, lofc] = OALS(D, iter, nPC, P)
% OALS  Perform Orthogonalized Alternating Least Squares (OALS) on the input data.
%
% REFERENCE:
%   Gómez-Sánchez, Adrián, et al. "Solving the missing value problem in 
%   PCA by Orthogonalized-Alternating Least Squares (O-ALS)." Chemometrics 
%   and Intelligent Laboratory Systems (2024): 105153.
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-01-02
% License:        MIT
% Reviewed by:    Lovelace's Square
% Version:        v 2.0
%
% This function implements Orthogonalized Alternating Least Squares (OALS) to
% decompose a data matrix D into scores (T) and loadings (P). The user
% specifies:
%   - D:    [R x C] data matrix (may contain NaNs).
%   - iter: Maximum number of OALS iterations (positive integer).
%   - nPC:  Desired number of components (positive integer).
%   - P:    (Optional) Initial guess for the loadings [nPC x C]. If empty,
%           OALS will generate 5 random guesses and pick the best solution
%           based on lowest final LOF (lack of fit).
%
% INPUTS:
%   D    - (R x C) Data matrix to be decomposed.
%   iter - Maximum number of iterations (positive integer).
%   nPC  - Number of components (positive integer).
%   P    - (Optional) Initial loadings [nPC x C]. If empty or omitted,
%          random initialization is used.
%
% OUTPUTS:
%   Dr   - (R x C) Reconstructed data matrix after OALS.
%   T    - (R x nPC) Scores matrix.
%   P    - (nPC x C) Final loadings.
%   r2   - (1 x #iters) Explained variance per iteration.
%   lofc - (1 x #iters) Lack of fit per iteration.
%
% DEPENDENCIES:
%   This function requires the following supporting scripts:
%     - ScoresLS.m  : Computes scores matrix from data and loadings.
%     - LoadingsLS.m: Computes loadings matrix from data and scores.
%     - lofNaN.m    : Calculates Lack of Fit (LOF) and explained variance.
%
% EXAMPLES:
%   % 1) Random initialization with nPC=3, pick best of 5 runs:
%   [Dr, T, P, r2, lofc] = OALS(D, 100, 3);
%
%   % 2) Provide your own loadings guess (nPC=2 must match size(P,1)):
%   P_init = rand(2, size(D,2));
%   [Dr, T, P, r2, lofc] = OALS(D, 100, 2, P_init);
%
% DISCLAIMER:
%   The authors and Lovelace's Square are not responsible for any damage,
%   inaccuracies, or data loss arising from the use of this function. 
%   Test thoroughly before using in production.

    %% --- 1. Validate inputs -------------------------------------------
    if nargin < 3
        error(['OALS requires at least 3 inputs: D, iter, and nPC.\n',...
               'Usage: OALS(D, iter, nPC, [P])']);
    end
    if nargin < 4
        P = [];  % If not provided, default to empty
    end

    if iter <= 0 || mod(iter, 1) ~= 0
        error('Number of iterations (iter) must be a positive integer.');
    end
    if nPC <= 0 || mod(nPC, 1) ~= 0
        error('Number of components (nPC) must be a positive integer.');
    end

    %% --- 2. Main logic -----------------------------------------------
    % If P is non-empty, ensure it matches the specified nPC
    if ~isempty(P)
        if size(P,1) ~= nPC
            error(['Provided loadings P have size(P,1)=%d, but nPC=%d.\n',...
                   'They must match for a valid OALS run.'], size(P,1), nPC);
        end
        % Single run with user-provided P
        [Dr, T, P, r2, lofc] = runOALS(D, iter, P);
    else
        % Multiple random runs (5) => pick best based on LOF
        nRuns = 5;
        bestLOF = Inf;
        bestSolution = struct('Dr',[],'T',[],'P',[],'r2',[],'lofc',[]);

        for runIdx = 1:nRuns
            % Build random initial loadings
            P_try = zeros(nPC, size(D,2));
            % First row is the col-wise mean of D (omitting NaNs)
            P_try(1,:) = mean(D, 1, 'omitnan');
            % Remaining rows random
            for rowIdx = 2:nPC
                P_try(rowIdx,:) = rand(1, size(D,2));
            end

            % Run OALS with this initial guess
            [Dr_try, T_try, P_try, r2_try, lofc_try] = runOALS(D, iter, P_try);

            % Check final LOF to see if best so far
            finalLOF = lofc_try(end);
            if finalLOF < bestLOF
                bestLOF = finalLOF;
                bestSolution.Dr   = Dr_try;
                bestSolution.T    = T_try;
                bestSolution.P    = P_try;
                bestSolution.r2   = r2_try;
                bestSolution.lofc = lofc_try;
            end
        end

        % Use best solution
        Dr   = bestSolution.Dr;
        T    = bestSolution.T;
        P    = bestSolution.P;
        r2   = bestSolution.r2;
        lofc = bestSolution.lofc;
    end
end % End of OALS function


%% ------------------------------------------------------------------------
% Subfunction: runOALS
%% ------------------------------------------------------------------------
function [Dr, T, P, r2, lofc] = runOALS(D, iter, P)
% RUNOALS  Single-run routine for OALS with a given initial loadings P.

    %% 1. Initialization
    it   = 0;
    conv = 99;   
    lofc(1) = 99;

    %% 2. Orthogonalize initial P
    a = P';
    q = P';
    for i = 2:size(P,1)
        p = 0;
        for j = 1:(i-1)
            denom = (q(:,j)' * q(:,j));
            if denom > 0
                p = p + (a(:,i)' * q(:,j)) / denom * q(:,j);
            end
        end
        q(:,i) = a(:,i) - p;
    end
    P = q';

    % Normalize each row in P
    for i = 1:size(P,1)
        rowNorm = norm(P(i,:));
        if rowNorm > 0
            P(i,:) = P(i,:) / rowNorm;
        end
    end

    %% 3. Iterations
    while abs(conv) > 1e-12 && it < iter
        it = it + 1;

        % --- a) Update Scores
        T = ScoresLS(D, P);

        % --- b) Orthogonalize T
        a = T;
        q = T;
        for i = 2:size(T,2)
            p = 0;
            for j = 1:(i-1)
                denom = (q(:,j)' * q(:,j));
                if denom > 0
                    p = p + (a(:,i)' * q(:,j)) / denom * q(:,j);
                end
            end
            q(:,i) = a(:,i) - p;
        end
        T = q;

        % --- c) Update Loadings
        P = LoadingsLS(D, T);

        % --- d) Orthogonalize P
        a = P';
        q = P';
        for i = 2:size(P,1)
            p = 0;
            for j = 1:(i-1)
                denom = (q(:,j)'*q(:,j));
                if denom > 0
                    p = p + (a(:,i)'*q(:,j)) / denom * q(:,j);
                end
            end
            q(:,i) = a(:,i) - p;
        end
        P = q';

        % --- e) Normalize loadings & rescale T accordingly
        for in = 1:size(P,1)
            lnorm = norm(P(in,:));
            if lnorm > 0
                T(:,in) = T(:,in) * lnorm;
                P(in,:) = P(in,:) / lnorm;
            end
        end

        % --- f) Compute r2 & LOF
        [r2(it), lofc(it)] = lofNaN(D, T, P);

        % --- g) Check convergence
        if it > 2
            conv = (lofc(it-1) - lofc(it)) / lofc(it);
        end
    end

    %% 4. Final reconstruction
    Dr = T * P;
end
