function [filteredMatrix] = pcaFilter(inputMatrix, numComponents)
% pcaFilter Performs PCA filtering on the input matrix.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%
% This function performs PCA filtering on the input matrix. It reduces the
% dimensionality of the data by projecting it onto the specified number of
% principal components and then reconstructs the matrix.
%
% INPUTS:
%    inputMatrix (array): Input matrix (nrows, ncols).
%    numComponents (int): Number of principal components to keep.
%
% OUTPUT:
%    filteredMatrix (array): Matrix reconstructed using the specified number
%    of principal components.
%    pcaParams (struct): Structure containing PCA parameters, including
%    eigenvectors and eigenvalues.
%
% Use example: [filteredMatrix, pcaParams] = pcaFilter(inputMatrix, 5);
%
% Disclaimer:
% Authors and Lovelace's Square are not responsible for any issues, inaccuracies, or data loss arising
% from the use of this function.

    % Check and initialize variables
    if nargin < 2
        error('Number of principal components must be specified.');
    end

    % Validate numComponents
    if numComponents <= 0 || numComponents > min(size(inputMatrix))
        error(['Number of principal components must be a positive integer ' ...
            'and less than or equal to the minimum dimension of the input matrix.']);
    end

    % Perform SVD (Singular Value Decomposition)
    [U, S, V] = svd(inputMatrix, 'econ');

    % Keep only the specified number of components
    U_reduced = U(:, 1:numComponents);
    S_reduced = S(1:numComponents, 1:numComponents);
    V_reduced = V(:, 1:numComponents);

    % Reconstruct the matrix using the reduced components
    filteredMatrix = U_reduced * S_reduced * V_reduced';
end
