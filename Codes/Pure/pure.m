function [pureProfiles, pureIndices] = pure(dataMatrix, numPureComponents, noisePercent, showPlots)
% PURE.  Legacy “weighted-SIMPLISMA” search for the purest variables.
%
% REFERENCES:
%   Windig, W.; Guilment, J. “Interactive Self-Modeling Mixture Analysis.”
%   Analytical Chemistry 63 (1991): 1425–1432.  — Original SIMPLISMA paper:
%   purity based on σ/(μ + noise) and independence via a determinant/
%   correlation-about-the-origin criterion to avoid collinearity.
%
%
% NOTE ON THIS IMPLEMENTATION:
%   This routine reproduces the legacy *pure.m* behavior widely circulated in the
%   mid-1990s: (i) classic SIMPLISMA purity p_j = σ_j/(μ_j + n), (ii) an extra
%   static weight w_j = (σ_j^2 + μ_j^2)/(σ_j^2 + (μ_j + n)^2) used in some early
%   scripts, and (iii) independence weighting via the determinant of a correlation-
%   about-the-origin submatrix formed from scaled columns (to penalize collinearity).
%
% Author    : Legacy SIMPLISMA pure.m implementations circulated from multiple sources; 
	      the most widely used open version appears in the Barcelona MCR-ALS toolboxes (Jaumot, Tauler, de Juan). 
              Proprietary implementations exist in PLS_Toolbox (Eigenvector).

Adrián Gómez-Sánchez polished and documented this legacy-style implementation.

Date of Creation: August 4, 2025

Reviewed by: Lovelace’s SquareAdrián Gómez-Sánchez  (legacy algorithm wrapped & documented)
% Date      : 2025-08-04
% License   : MIT
% Reviewed  : Lovelace’s Square – Yes
% Version   : 1.0
%
% DESCRIPTION
%   Given a data matrix **D** of size *nRows × nCols*, **PURE** extracts
%   **K** variables whose intensities are dominated by single chemical
%   components according to the historical weighting scheme:
%
%   • Unweighted purity                p   = σ / (µ + n)
%   • Legacy extra weight              w   = (σ² + µ²) / (σ² + (µ+n)²)
%   • Weighted purity (step 1)         P₁  = w · p
%   • For k ≥ 2 the candidate j that maximises
%
%        P_k(j) = P₁(j) · det ( R_k(j) ),                             (Eq. L1)
%
%     is selected, where R_k(j) is the k × k correlation sub-matrix built
%     from variables already chosen and candidate j **after** scaling each
%     column by l_j = √(σ² + (µ+n)²).
%
%   The routine returns the indices of the K pure variables (*pureIndices*)
%   and their raw profiles (*pureProfiles*), each row normalised to unit
%   maximum (as in the original script).
%
% MATHEMATICAL FORMULATION
%   • l_j             = √(σ_j² + (µ_j + n)²)
%   • Scaled data     dl(:,j) = D(:,j) / l_j
%   • Correlation     C = (dlᵀ dl) / nRows    (correlation-about-the-origin)
%   • Weight det      det( C(subset,subset) )
%
% INPUTS
%   dataMatrix         – Numeric matrix (nRows × nCols).
%   numPureComponents  – Positive integer K ≤ nCols.
%   noisePercent       – Noise level *f* (0–100 %), expressed relative to
%                        max (mean(dataMatrix)).
%   showPlots (opt.)   – true (default) to reproduce the diagnostic plots
%                        of the legacy script; false for silent mode.
%
% OUTPUTS
%   pureProfiles  – Numeric matrix.  Each **row** is the profile of one
%                   selected variable (nr × nRows if rows = spectra, or
%                   nr × nCols if rows = variables).
%   pureIndices   – 1 × K vector, column indices of the chosen variables.
%
% EXAMPLE
%   [sp, idx] = pure(X, 3, 1);          % 3 pure vars, 1 % noise threshold
%   figure, plot(sp'); title('Legacy PURE profiles');
%
% DISCLAIMER
%   This implementation is provided “as is.”  Neither the author nor
%   Lovelace’s Square accept liability for any loss or damage arising from
%   its use.
% -------------------------------------------------------------------------

    if nargin < 4, showPlots = true; end
    [nRows, nCols] = size(dataMatrix);

    %% Step 1 – purity spectrum (legacy definition) ----------------------
    n          = (noisePercent/100) * max(mean(dataMatrix,1));
    sigma      = std(dataMatrix,0,1);
    mu         = mean(dataMatrix,1);
    purity0    = sigma ./ (mu + n);                        % unweighted p_j

    % First pure variable
    [~, pureIndices(1)] = max(purity0);
    if showPlots
        fprintf('First pure variable : %d\n', pureIndices(1));
    end

    %% Step 2 – build legacy correlation matrix --------------------------
    l          = sqrt( sigma.^2 + (mu + n).^2 );
    dl         = dataMatrix ./ l;                          % column scaling
    C          = (dl' * dl) / nRows;                       % legacy “corr”

    %% Step 3 – initialise weights & spectra -----------------------------
    K          = numPureComponents;
    weights    = zeros(K, nCols);
    weights(1,:)= (sigma.^2 + mu.^2) ./ l.^2;              % extra w_j
    P          = zeros(K, nCols);
    P(1,:)     = weights(1,:) .* purity0;
    Sig        = zeros(K, nCols);
    Sig(1,:)   = weights(1,:) .* sigma;

    if showPlots
        figure('Name','PURE legacy – step 1');
        subplot(3,1,1), plot(mu),   title('Mean spectrum (µ_j)');
        subplot(3,1,2), plot(sigma),title('Std-dev spectrum (σ_j)');
        subplot(3,1,3), plot(P(1,:)),title('Weighted purity spectrum P₁');
    end

    %% Step 4 – iterative search for further pure variables --------------
    for k = 2:K
        for j = 1:nCols
            subM        = C([j, pureIndices(1:k-1)], [j, pureIndices(1:k-1)]);
            weights(k,j)= det(subM);
            P(k,j)      = P(1,j) .* weights(k,j);
            Sig(k,j)    = Sig(1,j) .* weights(k,j);
        end
        [~, pureIndices(k)] = max(P(k,:));
        if showPlots
            fprintf('Next pure variable  : %d\n', pureIndices(k));
            figure('Name',sprintf('PURE legacy – step %u',k));
            subplot(2,1,1), plot(P(k,:)),  title('Weighted purity spectrum');
            subplot(2,1,2), plot(Sig(k,:)),title('Weighted σ-spectrum');
        end
    end

    %% Step 5 – collect & normalise pure profiles ------------------------
    pureProfiles = zeros(K, nRows);
    for k = 1:K
        pureProfiles(k,:) = dataMatrix(:, pureIndices(k));
    end
    pureProfiles = rowNormMax(pureProfiles);               % legacy normv2
end
% ========================================================================

function A = rowNormMax(A)
% Scale each row to unit maximum (legacy normv2 behaviour)
    mx = max(abs(A),[],2);
    mx(mx==0) = 1;
    A = A ./ mx;
end
