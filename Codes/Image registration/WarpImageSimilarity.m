function imgW = WarpImageSimilarity(img, pOpt)
% WarpImageSimilarity. Apply similarity transformation to grayscale images 
% using backward warping with bilinear interpolation.
%
% REFERENCES:
% Implementation follows standard computer vision backward warping approach
% with bilinear interpolation for sub-pixel accuracy if needed.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2025-08-02
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The WarpImageSimilarity function applies a similarity transformation 
% (translation, rotation, scaling) to a grayscale input image using backward 
% warping. The transformation parameters are converted internally to the
% required inverse transformation matrix.
%
% For each output pixel, the function computes the corresponding source 
% coordinates using the inverse transformation, then uses bilinear 
% interpolation to determine the pixel value. Out-of-bounds regions are 
% filled with NaN values.
%
% The output size is automatically determined by the transformation to
% contain the entire warped image without cropping any data.
%
% INPUTS:
% img (array)  : A 2D numeric matrix (HxW) containing the grayscale
%                image to be transformed. Values should be in [0,1].
% pOpt (array) : A 1x4 double vector [tx ty thetaDeg scale] specifying
%                the transformation parameters.
%
% OUTPUTS:
% imgW (array) : A 2D numeric matrix containing the transformed
%                grayscale image with adaptive dimensions and NaN
%                for undefined regions.
%
% EXAMPLE:
% % Transform image using registration parameters:
% pOpt = [50 -20 15 1.2];  % [tx ty thetaDeg scale]
% warpedImg = WarpImageSimilarity(inputImg, pOpt);
%
% % Hyperspectral band transformation:
% for band = 1:nBands
%     aligned(:,:,band) = WarpImageSimilarity(hyperCube(:,:,band), pOpt);
% end
%
% DISCLAIMER:
% Authors and Lovelace's Square are not responsible for any issues,
% inaccuracies, or data loss arising from the use of this function.

% ----- transform matrices -----
T = paramsToMatrix(pOpt);
Tinv = inv(T);

% ----- sizes and output grid (fixed to input size) -----
if ndims(img) == 2
    [H_in, W_in] = size(img);
    [X, Y] = meshgrid(1:W_in, 1:H_in);  % fixed-size output grid

    % map output -> input
    xs = X .* Tinv(1,1) + Y .* Tinv(2,1) + Tinv(3,1);
    ys = X .* Tinv(1,2) + Y .* Tinv(2,2) + Tinv(3,2);

    % sample
    imgW = bilinearSample(img, xs, ys);

elseif ndims(img) == 3
    [H_in, W_in, nBands] = size(img);
    [X, Y] = meshgrid(1:W_in, 1:H_in);  % fixed-size output grid

    % map output -> input (compute ONCE and reuse)
    xs = X .* Tinv(1,1) + Y .* Tinv(2,1) + Tinv(3,1);
    ys = X .* Tinv(1,2) + Y .* Tinv(2,2) + Tinv(3,2);

    % allocate and warp each band using the SAME grid
    imgW = zeros(H_in, W_in, nBands, 'like', img);
    for band = 1:nBands
        if mod(band, 20) == 0 || band == nBands
            fprintf('  Processing band %d/%d\n', band, nBands);
        end
        imgW(:,:,band) = bilinearSample(img(:,:,band), xs, ys);
    end
else
    error('Input must be 2D (grayscale) or 3D (hyperspectral) image');
end
end

% ---------------------------- helpers ---------------------------------
function T = paramsToMatrix(p)
% 3x3 row-vector transform from [tx ty thetaDeg scale].
tx = p(1); ty = p(2); th = p(3); sc = p(4);
c = cosd(th); s = sind(th);
T = [ sc*c  sc*s  0; ...
     -sc*s  sc*c  0; ...
        tx    ty  1 ];
end

function Z = bilinearSample(I, xs, ys)
% 2D bilinear interpolation with out-of-bounds = NaN.
[H, W] = size(I);
x0 = floor(xs); y0 = floor(ys);
x1 = x0 + 1; y1 = y0 + 1;
dx = xs - x0; dy = ys - y0;

m00 = x0>=1 & x0<=W & y0>=1 & y0<=H;
m10 = x1>=1 & x1<=W & y0>=1 & y0<=H;
m01 = x0>=1 & x0<=W & y1>=1 & y1<=H;
m11 = x1>=1 & x1<=W & y1>=1 & y1<=H;

Z = NaN(size(xs), 'like', I);

V00 = zeros(size(xs), 'like', I); V00(m00) = I(sub2ind([H W], y0(m00), x0(m00)));
V10 = zeros(size(xs), 'like', I); V10(m10) = I(sub2ind([H W], y0(m10), x1(m10)));
V01 = zeros(size(xs), 'like', I); V01(m01) = I(sub2ind([H W], y1(m01), x0(m01)));
V11 = zeros(size(xs), 'like', I); V11(m11) = I(sub2ind([H W], y1(m11), x1(m11)));

W00 = (1 - dx) .* (1 - dy);
W10 = dx .* (1 - dy);
W01 = (1 - dx) .* dy;
W11 = dx .* dy;

val = W00.*V00 + W10.*V10 + W01.*V01 + W11.*V11;
anyMask = m00 | m10 | m01 | m11;
Z(anyMask) = val(anyMask);
end
