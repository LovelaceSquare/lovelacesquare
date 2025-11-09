function [smoothedMatrix] = WhittakerSmoother(inputMatrix, lambda, d)
% WhittakerSmoother. Apply the Whittaker smoother for signal smoothing.
%
% REFERENCES:
% Whittaker, Edmund T. "On a new method of graduation." 
% Proceedings of the Edinburgh Mathematical Society 41 (1922): 63-75.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%
% The Whittaker smoother is a signal processing method that reduces noise 
% while preserving the overall trends in the input data. It uses a penalized 
% least squares approach, controlled by a smoothing parameter (lambda) and 
% a difference operator of order (d). 
%
% Specifically, the smoothed solution is obtained by solving the linear system:
% (I + lambda * D' * D) * z = y, where:
%   - I is the identity matrix,
%   - D is a difference matrix of order 'd',
%   - y is the input signal.
%
% The output is a smoothed version of the input data, which captures broad 
% trends and suppresses noise or rapid fluctuations.
%
%
% INPUTS:
%    inputMatrix (array): A 2D array (nRows x nCols) containing the data 
%       to be smoothed. Each row (y) is treated independently as a 
%       separate signal.
%    lambda (double): A positive smoothing parameter that controls the 
%       degree of smoothness. Larger values of lambda yield smoother signals.
%    d (int): The order of the finite difference operator used in the 
%       penalty term. Commonly set to 2 for second-order differences to 
%       enforce curvature-based smoothing.
%
% OUTPUT:
%    smoothedMatrix (array): A 2D array (nRows x nCols) of the same size 
%       as inputMatrix, containing the smoothed signals. Each row corresponds 
%       to the smoothed version of the original input row.
%
% EXAMPLE:
%    Suppose you have a data matrix 'X' with each row representing a noisy signal:
%       lambdaValue = 1e3; % Adjust based on smoothing needs
%       diffOrder = 2; 
%       smoothedData = WhittakerSmoother(X, lambdaValue, diffOrder);
%    The resulting 'smoothedData' will contain smoothed versions of the 
%    input signals.
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues, 
%   inaccuracies, or data loss arising from the use of this function.

    % Validate input arguments
    if nargin < 3
        error('All input arguments (inputMatrix, lambda, d) must be specified.');
    end
    if lambda <= 0
        error('Smoothing parameter lambda must be positive.');
    end
    if d <= 0 || mod(d, 1) ~= 0
        error('Order of differences d must be a positive integer.');
    end

    % Get the size of the input matrix
    [nRows, nCols] = size(inputMatrix);
    
    % Initialize the output smoothed matrix
    smoothedMatrix = zeros(nRows, nCols);

    % Create the identity matrix and difference matrix for penalization
    E = speye(nCols);
    D = diff(E, d);
    penaltyMatrix = lambda * (D' * D);
    
    % Precompute the coefficient matrix (Z)
    Z = E + penaltyMatrix;

    % Apply the Whittaker smoother row-wise
    for i = 1:nRows
        y = inputMatrix(i, :).'; % Current row as a column vector
        smoothedRow = Z \ y;     % Solve the linear system
        smoothedMatrix(i, :) = smoothedRow.'; % Store result as a row
    end
end
