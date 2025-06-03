function [smoothedMatrix, imputedMatrix] = whittakerSmootherImpute(inputMatrix, lambda, d, maxIter, tol)
% WHITTAKERSMOOTHERIMPUTE  Apply the Whittaker smoother for signal smoothing,
%                          with iterative imputation of NaN values.
%
% REFERENCES:
%    Whittaker, Edmund T. "On a new method of graduation."
%    Proceedings of the Edinburgh Mathematical Society 41 (1922): 63-75.
%
% Author:       Adrián Gómez-Sánchez (modified by [Your Name])
% Date Created: 2024-12-14
% License:      MIT
% Reviewed by:  Lovelace's Square
% Version:      v 1.0
%
% DESCRIPTION:
%    The function applies Whittaker smoothing in an iterative fashion to
%    handle missing data (NaNs). It proceeds as follows:
%      1) Identifies missing values and imputes them (e.g., via linear interpolation).
%      2) Applies Whittaker smoothing to each row of the data.
%      3) Replaces missing values in the original data with the newly smoothed values.
%      4) Repeats Steps 2 and 3 until convergence or until maxIter iterations.
%
%    The Whittaker smoother is a penalized least squares method controlled by
%    parameters:
%        lambda (smoothness) and d (order of differences).
%    Larger lambda yields smoother signals, while d = 2 (commonly) enforces curvature-based smoothness.
%
% USAGE:
%    [smoothedMatrix, imputedMatrix] = whittakerSmootherImpute(inputMatrix, lambda, d, maxIter, tol);
%
% INPUTS:
%    inputMatrix (double, 2D):  The data matrix (nRows x nCols). Each row is
%                               treated as a separate 1D signal.
%    lambda (double):           Smoothing parameter, larger = smoother.
%    d (int):                   Order of finite difference (e.g., 1 or 2).
%    maxIter (int):             Maximum number of iterations for the
%                               iterative imputation process.
%    tol (double):              Convergence tolerance. If the change in
%                               imputed values between iterations is less
%                               than tol, the process stops.
%
% OUTPUTS:
%    smoothedMatrix (double, 2D): The final smoothed data matrix (nRows x nCols),
%                                 after the last iteration.
%    imputedMatrix (double, 2D):  The final data matrix with missing values
%                                 replaced by imputed (smoothed) values.
%                                 (In most cases, this is effectively the
%                                  same as smoothedMatrix, but kept separate
%                                  for clarity.)
%
% EXAMPLE:
%    data = [1 2 3 NaN 5; 2 3 NaN 5 6];
%    lambdaValue = 10;
%    diffOrder = 2;
%    maxIterations = 5;
%    tolerance = 1e-4;
%    [smoothData, finalImpute] = whittakerSmootherImpute(data, lambdaValue, diffOrder, ...
%                                                        maxIterations, tolerance);
%
% DISCLAIMER:
%    The author and Lovelace's Square are not responsible for inaccuracies
%    in the results derived from this function.

    % Validate input arguments
    if nargin < 5
        error('All input arguments (inputMatrix, lambda, d, maxIter, tol) must be specified.');
    end
    if lambda <= 0
        error('Smoothing parameter lambda must be positive.');
    end
    if d <= 0 || mod(d,1) ~= 0
        error('Order of differences d must be a positive integer.');
    end
    if maxIter <= 0 || mod(maxIter,1) ~= 0
        error('maxIter must be a positive integer.');
    end
    if tol <= 0
        error('tol must be a positive number.');
    end

    % Matrix dimensions
    [nRows, nCols] = size(inputMatrix);

    %----------------------------------------------------------------------
    % STEP 1: Initialize any missing values
    %         (Here we do a simple linear interpolation. Other methods
    %          could be used, e.g., spline, nearest, etc.)
    %----------------------------------------------------------------------

    % Copy the input to an imputation matrix
    imputedMatrix = inputMatrix;

    for iRow = 1:nRows
        rowData = imputedMatrix(iRow, :);
        nanMask = isnan(rowData);
        if all(nanMask)
            % If entire row is NaN, skip or fill with zeros (or some default).
            % For demonstration, fill with zero:
            imputedMatrix(iRow, :) = zeros(1, nCols);
        else
            % Perform linear interpolation for NaNs in the row
            validX = find(~nanMask);
            validY = rowData(~nanMask);
            allX   = 1:nCols;
            % Interpolate only for missing positions
            rowData(nanMask) = interp1(validX, validY, find(nanMask), 'linear', 'extrap');
            % Update imputed values in the matrix
            imputedMatrix(iRow, :) = rowData;
        end
    end

    %----------------------------------------------------------------------
    % STEP 2: Iterative Whittaker smoothing & re-imputation
    %----------------------------------------------------------------------

    % Precompute difference-based penalty for dimension nCols
    E = speye(nCols);
    D = diff(E, d);             % (nCols-d) x nCols
    penaltyMatrix = lambda * (D' * D);
    Z = E + penaltyMatrix;      % nCols x nCols

    % We will store the difference between iterations to check convergence
    prevImputed = imputedMatrix;

    for iter = 1:maxIter

        % Apply Whittaker smoothing row-wise using the current (imputed) data
        smoothedMatrix = zeros(nRows, nCols);
        for iRow = 1:nRows
            yRow = imputedMatrix(iRow, :).';

            % Solve the penalized system: (I + lambda * D'D) z = y
            zRow = Z \ yRow;
            smoothedMatrix(iRow, :) = zRow.';
        end

        % Re-impute missing values using the newly smoothed data
        %   i.e., if the original was NaN, replace with smoothed
        newImputedMatrix = imputedMatrix;  % start with the old imputed
        missingMask = isnan(inputMatrix);  % where original data is NaN
        newImputedMatrix(missingMask) = smoothedMatrix(missingMask);

        % Check convergence: max absolute difference across missing positions
        delta = max(abs(newImputedMatrix(missingMask) - prevImputed(missingMask)), [], 'all');

        % Update for next iteration
        imputedMatrix = newImputedMatrix;
        prevImputed   = imputedMatrix;

        % If below tolerance, stop
        if delta < tol
            fprintf('Converged at iteration %d (delta = %g).\n', iter, delta);
            break;
        end

    end

    % Final smoothing pass (optional, in case you want final "pure" smoothing):
    % If you want the final output to be smoothed (not just imputed) values:
    smoothedMatrix = zeros(nRows, nCols);
    for iRow = 1:nRows
        yRow = imputedMatrix(iRow, :).';
        zRow = Z \ yRow;
        smoothedMatrix(iRow, :) = zRow.';
    end

end
