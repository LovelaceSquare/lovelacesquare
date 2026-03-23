function out = imresize(M, scaleFactor, method)
% imresize - Resize 2D matrix (vanilla MATLAB implementation)
%
% SYNTAX:
%   out = imresize(M, scaleFactor)
%   out = imresize(M, scaleFactor, method)
%   out = imresize(M, [newH, newW])
%   out = imresize(M, [newH, newW], method)
%
% INPUTS:
%   M           - 2D or 3D numeric array
%   scaleFactor - Scalar scaling factor OR [newH, newW] target size
%   method      - Interpolation method (default: 'bicubic')
%                 Options: 'nearest', 'bilinear', 'bicubic'
%
% OUTPUTS:
%   out - Resized array
%
% DESCRIPTION:
%   Vanilla MATLAB implementation of imresize to avoid Image Processing
%   Toolbox dependency. Uses interp2 for interpolation.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

if nargin < 3
    method = 'cubic';  % interp2 uses 'cubic' instead of 'bicubic'
end

% Convert method names
switch method
    case 'bicubic'
        method = 'cubic';
    case 'bilinear'
        method = 'linear';
end

[H, W, C] = size(M);

% Determine new size
if isscalar(scaleFactor)
    newH = round(H * scaleFactor);
    newW = round(W * scaleFactor);
else
    newH = scaleFactor(1);
    newW = scaleFactor(2);
end

% Create coordinate grids
[X, Y] = meshgrid(1:W, 1:H);
[Xq, Yq] = meshgrid(linspace(1, W, newW), linspace(1, H, newH));

% Resize each channel
if C == 1
    out = interp2(X, Y, M, Xq, Yq, method);
else
    out = zeros(newH, newW, C);
    for c = 1:C
        out(:,:,c) = interp2(X, Y, M(:,:,c), Xq, Yq, method);
    end
end

% Preserve data type
if ~isa(M, 'double')
    out = cast(out, 'like', M);
end

end
