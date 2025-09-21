function [correctedData, baseline] = AsLS(data, lambda, p)
% AsLS. Perform Asymmetric Least Squares (AsLS) baseline correction on 
% spectral data.
%
% REFERENCES:
%   Eilers, Paul H.C., and Hans F.M. Boelens. 
%   "Baseline correction with asymmetric least squares smoothing." 
%   Leiden University Medical Centre Report 1.1 (2005): 5.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The AsLS function applies asymmetric least squares baseline correction to
% each row of a 2D data matrix, where each row is treated as an individual 
% spectrum or signal. The method iteratively updates a set of weights to 
% penalize large positive residuals more than negative ones, governed by 
% the asymmetry parameter (p), while controlling the baseline smoothness 
% through a penalty term scaled by lambda.
%
% Specifically, the baseline is obtained by solving a weighted least squares 
% problem:
%     min || w .* (data - baseline) ||^2 + lambda * ||D^2 baseline||^2
% where:
%   - w is an iteratively updated weight vector based on residuals,
%   - D is the second-order finite difference operator,
%   - lambda is the smoothing parameter,
%   - p is the asymmetry parameter.
%
% The function returns:
%   - correctedData: the original data with the estimated baseline subtracted,
%   - baseline: the estimated baseline for each row, matching the size of data.
%
% INPUTS:
%   data (array)   : A 2D numeric matrix (nRows x nCols). Each row represents 
%                    one spectrum or signal.
%   lambda (double): A positive smoothing parameter that controls how smooth 
%                    the baseline will be. Larger values => smoother baseline.
%   p (double)     : The asymmetry parameter in the range (0, 1). Smaller values 
%                    (e.g., 0.001) place greater emphasis on negative residuals.
%
% OUTPUTS:
%   correctedData (array): A 2D numeric matrix (nRows x nCols), the baseline-
%                          corrected signals (data - baseline).
%   baseline (array)     : A 2D numeric matrix (nRows x nCols), each row 
%                          containing the estimated baseline for the 
%                          corresponding row in 'data'.
%
% EXAMPLE:
%   % Suppose 'X' is your spectral data, each row is a spectrum:
%   lambdaVal = 1e6;   % Adjust based on data smoothness
%   pVal      = 0.001; % Typical value for asymmetry
%   [corrected, base] = AsLS(X, lambdaVal, pVal);
%   % 'corrected' now contains the baseline-corrected spectra, 
%   % 'base' contains the estimated baselines.
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues, 
%   inaccuracies, or data loss arising from the use of this function.

    % Validate input arguments
    if nargin < 3
        error('Usage: AsLS(data, lambda, p). All inputs are required.');
    end
    if ~isnumeric(data) || ndims(data) ~= 2 || isempty(data)
        error('Input "data" must be a nonempty 2D numeric matrix.');
    end
    if lambda <= 0
        error('The smoothing parameter "lambda" must be positive.');
    end
    if p <= 0 || p >= 1
        error('The asymmetry parameter "p" must be between 0 and 1.');
    end

    % Pre-allocate for speed
    [nRows, nCols] = size(data);
    baseline = zeros(nRows, nCols);
    correctedData = zeros(nRows, nCols);

    % Apply AsLS baseline correction row by row
    for r = 1:nRows
        y = data(r, :)';
        z = asls_baseline(y, lambda, p);
        baseline(r, :) = z';
    end

    % Subtract baseline from data
    correctedData = data - baseline;
end

%--------------------------------------------------------------------------
function z = asls_baseline(y, lambda, p)
% ASLS_BASELINE computes the Asymmetric Least Squares baseline for a single 
% row (spectrum) y using a fixed number of iterations.
%
%   y       : Column vector representing the spectrum.
%   lambda  : Positive smoothing parameter for baseline smoothness.
%   p       : Asymmetry parameter in (0, 1). 
%
% Returns:
%   z       : The estimated baseline for the spectrum y, same size as y.

    y = y(:);  % Ensure column vector
    N = length(y);

    % Second derivative operator for smoothing penalty
    D = diff(speye(N), 2);
    D2 = D' * D;

    % Initialize weights
    w = ones(N, 1);

    % Typically ~10 iterations suffice for AsLS
    numIterations = 10;
    for i = 1:numIterations
        W = spdiags(w, 0, N, N);
        % Solve for z: (W + lambda * D2) * z = W * y
        z = (W + lambda * D2) \ (W * y);

        % Update weights
        %   p weight if data above baseline => penalize positive residuals
        %   (1-p) weight if data below baseline
        w = p * (y > z) + (1 - p) * (y < z);
    end
end
