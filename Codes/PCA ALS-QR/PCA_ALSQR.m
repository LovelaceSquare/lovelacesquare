function [Dr, T, P, r2, lofc] = PCA_ALSQR(D, iter, nPC, P_init, useSVD)
% PCA_ALSQR  Perform Principal Component Analysis using Alternating Least
% Squares with final QR (and optional SVD) orthogonalization.
%
% Based on:
%   Gómez-Sánchez, Adrián, et al. "Solving the missing value problem in 
%   PCA by Orthogonalized-Alternating Least Squares (O-ALS)." Chemometrics 
%   and Intelligent Laboratory Systems (2024): 105153.
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-03-01
% License:        MIT
% Reviewed by:    Lovelace's Square
% Version:        v 1.0
%
% DESCRIPTION:
%   PCA_ALSQR performs principal component analysis on incomplete datasets
%   by iterative Alternating Least Squares (ALS) followed by a single
%   orthogonalization step using QR decomposition. This variant achieves 
%   the same PCA subspace as O-ALS, with reduced computational cost since 
%   orthogonality is enforced only once after ALS convergence.
%
% INPUTS:
%   D       - (R x C) Data matrix (may contain NaNs)
%   iter    - Maximum number of ALS iterations (positive integer)
%   nPC     - Number of principal components (positive integer)
%   P_init  - (Optional) Initial guess for loadings [nPC x C]
%   useSVD  - (Optional) Logical flag (default=true). If true, performs
%              additional SVD refinement after QR orthogonalization.
%
% OUTPUTS:
%   Dr   - (R x C) Reconstructed data matrix
%   T    - (R x nPC) Orthogonal scores matrix
%   P    - (nPC x C) Orthogonal loadings matrix
%   r2   - (1 x #iters) Explained variance per iteration
%   lofc - (1 x #iters) Lack of fit per iteration
%
% DEPENDENCIES:
%   Requires:
%     - ScoresLS.m
%     - LoadingsLS.m
%     - lofNaN.m
%
% EXAMPLE:
%   [Dr, T, P, r2, lofc] = PCA_ALSQR(D, 100, 3);
%   [Dr, T, P, r2, lofc] = PCA_ALSQR(D, 100, 3, [], false);
%
% DISCLAIMER:
%   The authors and Lovelace's Square assume no liability for the use of
%   this code. Validate results on your own datasets.

    %% --- 1. Validate inputs -------------------------------------------
    if nargin < 3
        error(['PCA_ALSQR requires at least 3 inputs: D, iter, and nPC.\n', ...
               'Usage: PCA_ALSQR(D, iter, nPC, [P_init], [useSVD])']);
    end
    if nargin < 4, P_init = []; end
    if nargin < 5, useSVD = true; end

    if iter <= 0 || mod(iter,1)~=0
        error('Number of iterations (iter) must be a positive integer.');
    end
    if nPC <= 0 || mod(nPC,1)~=0
        error('Number of components (nPC) must be a positive integer.');
    end

    %% --- 2. Initialization --------------------------------------------
    [R,C] = size(D);
    if isempty(P_init)
        P = zeros(nPC, C);
        P(1,:) = mean(D, 'omitnan');     % First component: dominant trend
        for i = 2:nPC
            P(i,:) = rand(1, C);
        end
    else
        P = P_init;
        if size(P,1) ~= nPC
            error('Provided P_init rows (%d) ≠ nPC (%d)', size(P,1), nPC);
        end
    end

    lofc(1) = 99;
    it = 0; conv = 99;

    %% --- 3. Alternating Least Squares Iterations ----------------------
    while abs(conv) > 1e-12 && it < iter
        it = it + 1;

        % (a) Update Scores
        T = ScoresLS(D, P);

        % (b) Update Loadings
        P = LoadingsLS(D, T);

        % (c) Compute explained variance and lack of fit
        [r2(it), lofc(it)] = lofNaN(D, T, P);

        % (d) Convergence criterion
        if it > 2
            conv = (lofc(it-1) - lofc(it)) / lofc(it);
        end
    end

%% --- 4. Final Orthogonalization ----------------------------------
[Q, Rq] = qr(T, 0);
M = Rq * P;

if useSVD
    [U, S, V] = svd(M, 'econ');

    % --- Reorder components by descending singular value ------------
    [svals, idx] = sort(diag(S), 'descend');
    S = diag(svals);
    U = U(:, idx);
    V = V(:, idx);

    % --- Build orthogonal scores and loadings -----------------------
    T = Q * U * sqrt(S);
    P = sqrt(S) * V';
else
    T = Q;
    P = M;
end


    %% --- 5. Normalization and Reconstruction --------------------------
    for i = 1:size(P,1)
        lnorm = norm(P(i,:));
        if lnorm > 0
            T(:,i) = T(:,i) * lnorm;
            P(i,:) = P(i,:) / lnorm;
        end
    end
    Dr = T * P;
end

