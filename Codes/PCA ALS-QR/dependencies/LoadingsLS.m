function [P] = LoadingsLS(D,T)
% LOADINGLS  Compute loadings matrix from data and scores (least squares)
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-01-02
% License:        MIT
% Reviewed by:    Lovelace's Square
%
% DESCRIPTION:
%   Computes variable loadings given current scores matrix T, ignoring NaNs.
%   Each variable j is solved independently using observed samples.
%
% INPUTS:
%   D - (I x J) data matrix
%   T - (I x K) scores
%
% OUTPUT:
%   P - (K x J) loadings matrix

[~,J] = size(D);
P = nan(size(T,2), size(D,2));
for j = 1:J
    ireal = find(isfinite(D(:,j)));
    if isempty(ireal)
        P(:,j) = nan;
    else
        P(:,j) = T(ireal,:) \ D(ireal,j);
    end
end
end
