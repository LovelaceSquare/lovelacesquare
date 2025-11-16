function [Dr,T,P] = NIPALS(D, it, n)
% Perform a NIPALS on the input matrix with missing data.
%
% Original publication: Wold, H. Nonlinear estimation by iterative least 
% square procedures. 1968.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 02 Jan 2025
% License: MIT
% Reviewed by Lovelace's Square: Yes
%
% Dependencies: ScoresLS, LoadingsLS, lofNaN
%
% Detailed function description:
% This function computes 'n' principal components from the data matrix 'D'
% (size R x C) using a NIPALS-like approach that handles missing data. We do:
%    1) For each component j = 1..n:
%        a) Initialize the loading profile P_j with the column means of D.
%        b) Iterate up to 'it' times:
%            - Calculate T_j using ScoresLS (skips NaNs).
%            - Calculate P_j using LoadingsLS (skips NaNs).
%            - Normalize P_j and rescale T_j accordingly.
%            - Check convergence based on change in P_j.
%        c) Deflate D by subtracting T_j * P_j.
%    2) Store T_j in T(:,j) and P_j in P(j,:).
%    3) Build the reconstructed matrix Dr = T * P at the end.
%
% Args:
%    D (matrix): R x C data matrix, possibly with NaNs.
%    it (int):   Maximum iterations per component.
%    n (int):    Number of components to extract.
%
% Returns:
%    Dr (matrix): Reconstructed data (R x C) using the 'n' extracted components.
%    T (matrix):  Scores (R x n).
%    P (matrix):  Loadings (n x C).
%
% Example:
%   [Dr, T, P] = NIPALS(D, 100, 3);
%
% Lovelace's Square is not responsible for any issues or errors
% that may arise from the use of this function. Use it at your own risk.

% Validate inputs
if nargin < 3
    error('All input arguments (D, it, n) must be specified.');
end
if it <= 0 || mod(it,1)~=0
    error('it (iterations) must be a positive integer.');
end
if n <= 0 || mod(n,1)~=0
    error('n (number of components) must be a positive integer.');
end

[R,C] = size(D);

% Pre-allocate results
T = zeros(R, n);
P = zeros(n, C);

% Loop over each component
for j = 1:n
    
    % 1) Initialize loadings for component j with mean of columns (ignoring NaNs)
    Px = mean(D, 1, 'omitnan');  % 1 x C
    if norm(Px) < eps
        warning('Initial P for component %d is near zero. Skipping.', j);
        continue
    end
    Px = Px / norm(Px);  % Normalize
    
    % Keep an old copy for convergence check
    oldPx = Px;
    
    % 2) Iterative refinement of T_j and P_j
    convThreshold = 1e-12;  % You can adjust as needed
    for iIter = 1:it
        
        % a) Calculate T_j given P_j
        Tx = ScoresLS(D, Px);
        
        % b) Calculate P_j given T_j
        Px = LoadingsLS(D, Tx);
        if norm(Px) < eps
            warning('Norm of loading for component %d is zero. Breaking early.', j);
            break;
        end
        
        % c) Normalize P_j, rescale T_j
        pNorm = norm(Px);
        Px = Px / pNorm;
        Tx = Tx * pNorm;
        
        % d) Convergence check on P_j
        if norm(Px - oldPx) < convThreshold
            break;
        end
        oldPx = Px;
    end
    
    % 3) Deflate the data
    %    D <- D - T_j * P_j
    for ccol = 1:C
        validRows = ~isnan(D(:, ccol));
        D(validRows, ccol) = D(validRows, ccol) - Tx(validRows) * Px(ccol);
    end
    
    % 4) Store final T_j, P_j
    T(:, j) = Tx;
    P(j, :) = Px;
    
end

% 5) Final reconstruction
Dr = T * P;

end
