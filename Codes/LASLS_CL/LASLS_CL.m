function [baseline, weights] = LASLS_CL( ...
    y,                 ... % (n x 1) data vector
    intervals,         ... % (m x 2) or cell of intervals; each row [startIdx, endIdx]
    pVals,             ... % (m x 1) each interval's p_j
    lambdasAsym,       ... % (m x 1) each interval's "asym" lambda_j
    lambdaWhit,        ... % default smoothing penalty outside intervals
    mu,                ... % global first-derivative penalty
    maxIter,           ... % IRLS max iterations
    tol)               ... % IRLS stopping tolerance (relative baseline change)
% LAsLS. Perform Local Asymmetric Least Squares (LAsLS) baseline correction with per-interval parameters.
%
% REFERENCES:
%   Eilers, Paul H.C., and Hans F.M. Boelens. 
%   "Baseline correction with asymmetric least squares smoothing." 
%   Leiden University Medical Centre Report 1.1 (2005): 5.
%
% Authors: Adrián Gómez-Sánchez, Berta Torres-Cobos, Rodrigo Rocha de Oliveira
% Date Created: 2024-12-16
% License: MIT
% Repository: https://github.com/LovelaceSquare/lovelacesquare
% Reviewed by Lovelace's Square: Yes
% Version: 1.1
%
% The Local Asymmetric Least Squares (LAsLS) algorithm is a modification of the
% standard Asymmetric Least Squares (AsLS) smoothing method, adapted for
% one-dimensional data vectors. It permits distinct asymmetry parameters
% (`pVals`) and second-derivative smoothing penalties (`lambdasAsym`) on user-defined
% intervals. Outside these intervals, a uniform smoothing penalty (`lambdaWhit`)
% and a global first-derivative penalty (`mu`) are enforced. Through an Iteratively
% Reweighted Least Squares (IRLS) procedure, weights are updated based on residuals
% at each iteration to produce a smooth, locally tailored baseline.
%
% Specifically, for each data point, the baseline is estimated by solving:
%       (W + D' * diag(lambdaVec) * D + mu * L' * L) * baseline = W * y
% where:
%   - y is the input data vector,
%   - W is a diagonal weight matrix that is updated iteratively,
%   - D is the second-order finite difference operator,
%   - L is the first-order finite difference operator,
%   - lambdaVec contains local smoothing parameters (set to the maximum of overlapping lambdas),
%   - mu is the global first-derivative penalty,
%   - pVals and lambdasAsym define the local asymmetry and smoothing inside each specified interval.
%
% The function returns:
%   - baseline: an (n x 1) vector containing the estimated baseline.
%   - weights : an (n x 1) vector with the final IRLS weights.
%
% INPUTS:
%   y           - (n x 1) numeric vector representing the data (e.g., intensity values).
%   intervals   - either an (m x 2) numeric array or a cell array of [startIdx, endIdx] indices defining intervals.
%   pVals       - (m x 1) numeric vector; each element is the asymmetry parameter p_j for the corresponding interval.
%   lambdasAsym - (m x 1) numeric vector; each element is the local smoothing penalty lambda_asym_j for the corresponding interval.
%   lambdaWhit  - scalar, the smoothing penalty used outside any specified interval.
%   mu          - scalar, the global first-derivative penalty.
%   maxIter     - scalar, the maximum number of iterations allowed for the IRLS algorithm.
%   tol         - scalar, the stopping tolerance based on the relative change in the baseline estimate.
%
% OUTPUTS:
%   baseline    - (n x 1) numeric vector containing the estimated baseline.
%   weights     - (n x 1) numeric vector containing the final IRLS weights.
%
% EXAMPLE:
%   y = rand(100,1)*10;
%   intervals = [10 20; 40 50];
%   pVals = [0.01; 0.001];
%   lambdasAsym = [1e4; 2e5];
%   lambdaWhit = 100;
%   mu = 1e3;
%   maxIter = 50;
%   tol = 1e-6;
%   [baseline, weights] = LASLS_CL(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);
%
% DISCLAIMER:
%   The author and Lovelace's Square are not responsible for any issues, inaccuracies, 
%   or data loss arising from the use of this function.

    %% 1) Input checks & setup
    n = length(y);
    if size(y,2) > 1
        error('Input y must be a single column or row vector.');
    end

    m = length(pVals);
    if length(lambdasAsym) ~= m
        error('pVals and lambdasAsym must match in length.');
    end

    % Convert intervals to a consistent (m x 2) numeric array if it's a cell
    if iscell(intervals)
        % E.g. intervals{k} = [startIdx, endIdx]
        intMat = zeros(m,2);
        for k = 1:m
            intMat(k,1) = intervals{k}(1);
            intMat(k,2) = intervals{k}(2);
        end
    else
        intMat = intervals;
    end

    %% 2) Build arrays for local parameters for each index
    % pVec: local asymmetry parameters; value -1 indicates outside any interval.
    pVec = -ones(n,1);
    % lambdaVec: local second-derivative smoothing penalty for the "center" of the operator.
    lambdaVec = lambdaWhit * ones(n-2,1);  % default value outside intervals

    % For each interval j, update pVec and lambdaVec accordingly.
    for j = 1:m
        iStart = intMat(j,1);
        iEnd   = intMat(j,2);
        this_p = pVals(j);
        this_lam = lambdasAsym(j);

        % Sanity check for interval boundaries
        if iStart < 1 || iEnd > n || iStart > iEnd
            warning('Skipping invalid interval [%d, %d].', iStart, iEnd);
            continue;
        end

        % Update pVec for indices within the interval (overwrite if overlapping)
        pVec(iStart:iEnd) = this_p;

        % Update lambdaVec for rows corresponding to indices within the interval.
        % Note: row i in lambdaVec corresponds to data index i+1.
        rowStart = max(1, iStart-1);
        rowEnd   = min(n-2, iEnd-1);
        for r = rowStart:rowEnd
            lambdaVec(r) = max(lambdaVec(r), this_lam);
        end
    end

    %% 3) Construct difference operators
    % Second derivative operator, size = (n-2) x n.
    e = ones(n,1);
    D = spdiags([e, -2*e, e], 0:2, n-2, n);

    % First derivative operator, size = (n-1) x n.
    L = spdiags([-ones(n,1), ones(n,1)], [0,1], n-1, n);

    %% 4) IRLS initialization
    weights = ones(n,1);
    baseline = zeros(n,1);

    %% 5) Iterative Reweighted Least Squares (IRLS)
    for iter = 1:maxIter
        % Diagonal weight matrix
        W = spdiags(weights, 0, n, n);

        % Construct normal equation matrix:
        %   A = W + D' * diag(lambdaVec) * D + mu * (L' * L)
        A = W + D' * spdiags(lambdaVec, 0, n-2, n-2) * D + mu * (L' * L);
        rhs = W * y;

        % Solve for new baseline estimate
        newBaseline = A \ rhs;

        % Check relative change for convergence
        relChange = norm(newBaseline - baseline, 2) / (norm(baseline, 2) + eps);
        baseline = newBaseline;
        if relChange < tol
            break;
        end

        % Update weights only for indices within specified intervals
        % For points outside intervals (pVec < 0), weight remains 1.
        for i = 1:n
            if pVec(i) > 0
                if y(i) > baseline(i)
                    weights(i) = pVec(i);
                else
                    weights(i) = 1 - pVec(i);
                end
            else
                weights(i) = 1;
            end
        end
    end
end
