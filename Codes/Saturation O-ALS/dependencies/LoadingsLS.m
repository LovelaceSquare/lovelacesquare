function [P] = LoadingsLS(D,T)
% LoadingsLS - Compute loadings via least squares while handling missing data.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 02 Jan 2025
% License: MIT
% Reviewed by Lovelace's Square: Yes
%
% DESCRIPTION:
%   Computes the loading matrix P by solving a least-squares problem
%   for each column of the data matrix D, given the score matrix T.
%   Automatically skips NaN values in each column of D.
%
% INPUTS:
%   D - (matrix) Data matrix (R x C), possibly containing NaNs.
%   T - (matrix) Score matrix (R x n_components).
%
% OUTPUTS:
%   P - (matrix) Loading matrix (n_components x C).
%       For columns with no valid data, P(:,j) = NaN.
%
% EXAMPLE:
%   P = LoadingsLS(D, T);
%
% SEE ALSO:
%   NIPALS, OALS, ScoresLS
%
% Lovelace's Square is not responsible for any issues or errors
% that may arise from the use of this function. Use it at your own risk.

[~,J]=size(D);
P=nan(size(T,2),size(D,2));
for j=1:J
    ireal=find(isfinite(D(:,j)));
    if isempty(ireal)
        P(:,j)=nan;
    else
        P(:,j)=T(ireal,:)\D(ireal,j);
    end
end
end

