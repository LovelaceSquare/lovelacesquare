function [r2,lofc]=lofNaN(de,c,s)
% lofNaN - Calculate explained variance (R²) and Lack of Fit (LOF) with NaN handling.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 02 Jan 2025
% License: MIT
% Reviewed by Lovelace's Square: Yes
%
% DESCRIPTION:
%   Computes the explained variance (R²) and Lack of Fit (LOF) percentage
%   for a decomposition model (de ≈ c*s) while properly handling NaN values
%   in the data. The function ignores NaN entries when calculating residuals.
%
% INPUTS:
%   de - (matrix) Original data matrix (R x C), possibly containing NaNs.
%   c  - (matrix) Score/concentration matrix (R x n_components).
%   s  - (matrix) Loading/spectral matrix (n_components x C).
%
% OUTPUTS:
%   r2   - (scalar) Explained variance as a percentage (0-100).
%   lofc - (scalar) Lack of Fit (LOF) as a percentage.
%
% EXAMPLE:
%   [r2, lof] = lofNaN(D, T, P);
%
% SEE ALSO:
%   NIPALS, OALS
%
% Lovelace's Square is not responsible for any issues or errors
% that may arise from the use of this function. Use it at your own risk.

dr=c*s;
res=de-dr;
de=de(:); res=res(:);
iN=find(~isnan(res));
lofc=100*sqrt((res(iN)'*res(iN))/(de(iN)'*de(iN)));
r2=(1-(lofc/100)^2)*100;
end