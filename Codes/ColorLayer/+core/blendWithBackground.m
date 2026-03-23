function blendedImage = blendWithBackground(compositeImage, Original, alpha)
% blendWithBackground - Blend composite with background image
%
% SYNTAX:
%   blendedImage = blendWithBackground(compositeImage, Original, alpha)
%
% INPUTS:
%   compositeImage - RGB image [0,1] double (H x W x 3)
%   Original       - Background image (H x W) or (H x W x 3)
%   alpha          - Blending factor [0,1]
%                    0 = pure background, 1 = pure composite
%
% OUTPUTS:
%   blendedImage - Blended RGB image [0,1] double (H x W x 3)
%
% DESCRIPTION:
%   Blends composite image with background using:
%   blendedImage = alpha * compositeImage + (1-alpha) * background
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Validate inputs
if nargin < 3
    alpha = 0.5;
end

if isempty(Original)
    blendedImage = compositeImage;
    return;
end

%% Get target size
[H, W, ~] = size(compositeImage);

%% Prepare background
bg = double(Original);

% Resize if needed
if ~isequal(size(bg, 1), H) || ~isequal(size(bg, 2), W)
    bg = utils.imresize(bg, [H, W], 'bicubic');
end

% Convert to RGB if grayscale
if size(bg, 3) == 1
    bg = repmat(bg, [1, 1, 3]);
end

% Normalize to [0, 1]
bg = utils.mat2gray(bg);

%% Blend
blendedImage = alpha * compositeImage + (1 - alpha) * bg;

% Ensure [0, 1] range
blendedImage = max(0, min(1, blendedImage));

end
