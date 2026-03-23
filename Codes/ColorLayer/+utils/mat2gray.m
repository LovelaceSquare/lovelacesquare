function out = mat2gray(M)
% mat2gray - Normalize matrix to [0, 1] range (vanilla MATLAB implementation)
%
% SYNTAX:
%   out = mat2gray(M)
%
% INPUTS:
%   M - Numeric array
%
% OUTPUTS:
%   out - Normalized array in [0, 1]
%
% DESCRIPTION:
%   Vanilla MATLAB implementation of mat2gray to avoid Image Processing
%   Toolbox dependency.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

minVal = min(M(:));
maxVal = max(M(:));

if maxVal == minVal
    % Constant matrix
    out = zeros(size(M));
else
    out = (M - minVal) / (maxVal - minVal);
end

end
