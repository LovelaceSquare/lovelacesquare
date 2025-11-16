function [C] = ScoresLS(D,S)
% ScoresLS - Compute scores via least squares while handling missing data.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 02 Jan 2025
% License: MIT
% Reviewed by Lovelace's Square: Yes
%
% DESCRIPTION:
%   Computes the score matrix C by solving a least-squares problem
%   for each row of the data matrix D, given the loading matrix S.
%   Automatically skips NaN values in each row of D.
%
% INPUTS:
%   D - (matrix) Data matrix (R x C), possibly containing NaNs.
%   S - (matrix) Loading matrix (n_components x C).
%
% OUTPUTS:
%   C - (matrix) Score matrix (R x n_components).
%       For rows with no valid data, C(i,:) = NaN.
%
% EXAMPLE:
%   C = ScoresLS(D, S);
%
% SEE ALSO:
%   NIPALS, OALS, LoadingsLS
%
% Lovelace's Square is not responsible for any issues or errors
% that may arise from the use of this function. Use it at your own risk.

[I,~]=size(D);
C=nan(size(D,1),size(S,1));
for i=1:I
    jreal=find(isfinite(D(i,:)));
    if isempty(jreal)
        C(i,:)=nan;
    else
        C(i,:)=D(i,jreal)/S(:,jreal);
    end
end
end

