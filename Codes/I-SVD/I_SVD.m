function [Dimp, T, P, r2, lofc] = I_SVD(D, nComp, maxIter, tol)
% I_SVD  Perform Iterative SVD-based PCA Imputation on the input data.
%
% REFERENCE:
%   Hastie, Trevor, et al. "Imputing missing data for gene expression arrays." (1999).
%   Technical Report, Division of Biostatistics, Stanford University
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-01-02
% License:        MIT
% Reviewed by:    Lovelace's Square
% Version:        v1.0
%
% This function implements an iterative Singular Value Decomposition (SVD)
% approach to handle missing data in a data matrix 'D' (size R x C). It replaces
% NaNs with initial estimates, computes a rank-nComp approximation of the fully
% imputed matrix via SVD, then updates the missing entries, repeating until
% convergence is achieved based on a specified tolerance or until 'maxIter' is reached.
%
% INPUTS:
%   D       - (R x C) Data matrix with possible missing values (NaNs).
%   nComp   - (1 x 1) Number of components to retain in the SVD-based approximation (positive integer).
%   maxIter - (1 x 1) Maximum number of global iterations for imputation (positive integer).
%   tol     - (1 x 1, optional) Tolerance for convergence (default: 1e-12).
%
% OUTPUTS:
%   Dimp - (R x C) Final imputed data matrix after iterative SVD.
%   T    - (R x nComp) Scores matrix representing the projection of data onto principal components.
%   P    - (nComp x C) Loadings matrix representing the principal components.
%   r2   - (1 x maxIter) Explained variance for each iteration.
%   lofc - (1 x maxIter) Lack of fit at each iteration.
%
% DEPENDENCIES:
%   None.
%
% EXAMPLES:
%   % 1) Basic usage with 3 components and 50 maximum iterations:
%   [Dimp, T, P, r2, lofc] = I_SVD(D, 3, 50);
%
%   % 2) Specify a different tolerance for convergence:
%   [Dimp, T, P, r2, lofc] = I_SVD(D, 5, 100, 1e-10);
%
% DISCLAIMER:
%   The authors and Lovelace's Square are not responsible for any damage,
%   inaccuracies, or data loss arising from the use of this function. 
%   Test thoroughly before using in production.

    %% --- 1. Validate Inputs -------------------------------------------
    if nargin < 3
        error(['I_SVD requires at least 3 inputs: D, nComp, and maxIter.\n',...
               'Usage: I_SVD(D, nComp, maxIter, [tol])']);
    end
    if nargin < 4
        tol = 1e-12;  % Default tolerance
    end
    
    if ~isnumeric(D) || ~ismatrix(D)
        error('Input D must be a numeric matrix.');
    end
    [R, C] = size(D);
    
    if ~isscalar(nComp) || nComp <= 0 || floor(nComp) ~= nComp
        error('Number of components (nComp) must be a positive integer.');
    end
    if nComp > min(R, C)
        error('Number of components (nComp) cannot exceed the smaller dimension of D.');
    end
    
    if ~isscalar(maxIter) || maxIter <= 0 || floor(maxIter) ~= maxIter
        error('Maximum number of iterations (maxIter) must be a positive integer.');
    end
    
    if ~isscalar(tol) || tol <= 0
        error('Tolerance (tol) must be a positive scalar.');
    end

    %% --- 2. Initialization -------------------------------------------
    missingMask = isnan(D);
    Dimp = D;
    
    % Initial imputation by column means
    for c = 1:C
        colValid = ~missingMask(:,c);
        if any(colValid)
            colMean = mean(Dimp(colValid, c));
            Dimp(~colValid, c) = colMean;
        else
            Dimp(~colValid, c) = 0;  % If entire column is NaN, set to zero
        end
    end
    
    r2 = zeros(1, maxIter);
    lofc = zeros(1, maxIter);
    
    prevLOF = Inf;
    
    %% --- 3. Iterative Imputation Loop -------------------------------
    for it = 1:maxIter
        % --- a) SVD on the fully imputed matrix
        [U, S, V] = svd(Dimp, 'econ');
        
        % Truncate to nComp components
        U_n = U(:, 1:nComp);
        S_n = S(1:nComp, 1:nComp);
        V_n = V(:, 1:nComp);
        
        % --- b) Reconstruct the rank-nComp approximation
        Dapprox = U_n * S_n * V_n';
        
        % --- c) Update missing values only
        Dimp(missingMask) = Dapprox(missingMask);
        
        % --- d) Compute r2 & LOF
        % Assuming lofNaN calculates R^2 and lack of fit based on original D and imputed Dimp
        [r2(it), lofc(it)] = lofNaN(D, Dapprox);
        
        % --- e) Check for convergence
        if it > 1
            convVal = abs(lofc(it-1) - lofc(it)) / max(eps, lofc(it));
            if convVal < tol
                % Trim the results if convergence is achieved early
                r2 = r2(1:it);
                lofc = lofc(1:it);
                break;
            end
        end
        
        prevLOF = lofc(it);
    end
    
    %% --- 4. Final Outputs ---------------------------------------------
    T = U_n * S_n;
    P = V_n';
    
end % End of I_SVD function

%% ------------------------------------------------------------------------
% Subfunction: lofNaN
% -------------------------------------------------------------------------
function [r2, lofc] = lofNaN(D_original, D_approx)
% LOFNaN  Calculate R^2 and Lack of Fit (LOF) between original and approximated data.
%
    % Calculate residuals only where original data is not NaN
    validMask = ~isnan(D_original);
    residuals = D_original(validMask) - D_approx(validMask);
    ss_res = sum(residuals.^2);
    ss_tot = sum((D_original(validMask) - mean(D_original(validMask))).^2);
    
    r2 = 1 - ss_res / ss_tot;
    lofc = ss_res;
end
