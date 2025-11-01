function [r2, lofc] = lofNaN(de, c, s)
% LOFNAN  Compute explained variance (r²) and lack of fit (LOF)
%          handling missing values (NaNs).
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-01-02
% License:        MIT
% Reviewed by:    Lovelace's Square
%
% INPUTS:
%   de - Data matrix (observed)
%   c  - Scores matrix
%   s  - Loadings matrix
%
% OUTPUTS:
%   r2   - Explained variance (%)
%   lofc - Lack of fit (%)

dr = c * s;
res = de - dr;
de = de(:); res = res(:);
iN = find(~isnan(res));
lofc = 100 * sqrt((res(iN)' * res(iN)) / (de(iN)' * de(iN)));
r2   = (1 - (lofc/100)^2) * 100;
end
