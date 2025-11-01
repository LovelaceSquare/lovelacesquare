function [C] = ScoresLS(D,S)
% SCORESLS  Compute scores matrix from data and loadings (least squares)
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-01-02
% License:        MIT
% Reviewed by:    Lovelace's Square
%
% DESCRIPTION:
%   Computes sample scores given current loadings matrix S, ignoring NaNs.
%   Each sample i is solved independently by least squares over its
%   observed variables.
%
% INPUTS:
%   D - (I x J) data matrix
%   S - (K x J) loadings
%
% OUTPUT:
%   C - (I x K) scores matrix

[I,~] = size(D);
C = nan(size(D,1), size(S,1));
for i = 1:I
    jreal = find(isfinite(D(i,:)));
    if isempty(jreal)
        C(i,:) = nan;
    else
        C(i,:) = D(i,jreal) / S(:,jreal);
    end
end
end
