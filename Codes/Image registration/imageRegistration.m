function [registeredImg, pOpt] = imageRegistration( ...
    movingImg, fixedImg, p0, lb, ub, doShow, optScale, useInteractive)
% imageRegistration. Perform bounded 2-D similarity registration between 
% grayscale images using toolbox-free optimization.
%
% REFERENCES:
% Sara Piqueras et al. "Handling Different Spatial Resolutions in Image 
% Fusion by Multivariate Curve Resolution–Alternating Least Squares for 
% Incomplete Image Multisets." Analytical Chemistry 2018, 90(11), 6757–6765. 
% DOI: 10.1021/acs.analchem.8b00630.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2025-08-02
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The imageRegistration function estimates optimal similarity transformation
% parameters [tx ty thetaDeg scale] that align a moving grayscale image to 
% a fixed reference grayscale image. The method uses bounded optimization 
% with logit transformation to minimize NaN-masked mean squared error (MSE) 
% between the warped moving image and fixed image.
%
% The transformation follows row-vector, post-multiply convention:
% [x' y' 1] = [x y 1] * T, where T(1:2,1:2) = s * [cos sin; -sin cos].
% Non-overlapping regions are handled by filling with NaN values during 
% warping operations.
%
% The function supports both minimal 2-argument calls (with interactive 
% initialization and loose bounds) and full parameter specification for 
% complete control over the optimization process.
%
% The output registered image uses adaptive sizing to contain the entire
% transformed image.
%
% INPUTS:
% movingImg (array)     : A 2D numeric matrix (HxW) containing the grayscale
%                         image to be transformed to align with the fixed image.
% fixedImg (array)      : A 2D numeric matrix (HxW) containing the grayscale
%                         reference image that movingImg will be aligned to.
% p0 (array, optional)  : A 1x4 double vector [tx ty thetaDeg scale] providing
%                         initial parameter estimates. Default: [0 0 0 1].
% lb (array, optional)  : A 1x4 double vector [tx ty thetaDeg scale] specifying
%                         lower bounds. Default: very loose bounds.
% ub (array, optional)  : A 1x4 double vector [tx ty thetaDeg scale] specifying
%                         upper bounds. Default: very loose bounds.
% doShow (logical, optional)      : Display false-color overlay (default true).
% optScale (logical, optional)    : Optimize scale parameter (default true).
% useInteractive (logical, optional): Enable interactive point-based
%                         initialization (default true for 2-arg calls).
%
% OUTPUTS:
% registeredImg (array) : A 2D numeric matrix containing the moving image 
%                         warped to align with the fixed image. Output uses
%                         adaptive sizing with NaN for out-of-bounds regions.
% pOpt (array)          : A 1x4 double vector [tx ty thetaDeg scale] containing
%                         the optimal transformation parameters found by
%                         optimization.
%
% EXAMPLE:
% % Minimal interactive call with loose bounds:
% [regImg, pOpt] = imageRegistration(movingGray, fixedGray);
%
% % Full control call:
% p0 = [0 0 0 1];  lb = [-100 -100 -45 0.5];  ub = [100 100 45 2.0];
% [regImg, pOpt] = imageRegistration(movingGray, fixedGray, p0, lb, ub, ...
%                                   true, true, false);
%
% % Silent registration with custom bounds:
% [regImg, pOpt] = imageRegistration(movingGray, fixedGray, [0 0 0 1], ...
%                                   [-50 -50 -30 0.8], [50 50 30 1.2], false);
%
% DEPENDENCIES:
% Requires WarpImageSimilarity.m
%
% DISCLAIMER:
% Authors and Lovelace's Square are not responsible for any issues,
% inaccuracies, or data loss arising from the use of this function.

    % -------------------- Handle different calling modes ---------------
    if nargin == 2
        % Minimal call: use interactive mode with loose bounds
        p0 = [0, 0, 0, 1];                    % No transformation initial guess
        lb = [-5000, -5000, -360, 0.1];      % Very loose lower bounds  
        ub = [ 5000,  5000,  360, 10.0];     % Very loose upper bounds
        doShow = true;                        % Show results
        optScale = true;                      % Optimize scale
        useInteractive = true;                % Interactive initialization
        
        fprintf('\n=== Interactive Registration Mode ===\n');
        fprintf('Using loose bounds and interactive initialization.\n');
        fprintf('Bounds: tx∈[%.0f,%.0f], ty∈[%.0f,%.0f], θ∈[%.0f°,%.0f°], scale∈[%.1f,%.1f]\n', ...
                lb(1), ub(1), lb(2), ub(2), lb(3), ub(3), lb(4), ub(4));
    else
        % Full specification mode - handle individual defaults
        if nargin < 3 || isempty(p0)
            p0 = [0, 0, 0, 1];
        end
        if nargin < 4 || isempty(lb)
            lb = [-1000, -1000, -360, 0.1];
        end
        if nargin < 5 || isempty(ub)
            ub = [ 1000,  1000,  360, 10.0];
        end
        if nargin < 6 || isempty(doShow)
            doShow = true;
        end
        if nargin < 7 || isempty(optScale)
            optScale = true;
        end
        if nargin < 8 || isempty(useInteractive)
            useInteractive = false;  % Default to non-interactive for full calls
        end
    end

    validateInputs(movingImg, fixedImg, p0, lb, ub);

    % Convert to double in [0,1]
    movingGray = toDouble01(movingImg);
    fixedGray = toDouble01(fixedImg);

    outSize = [size(fixedGray,1), size(fixedGray,2)];

    % -------------------- Interactive initialization -------------------
    if useInteractive
        fprintf('\n=== Interactive Parameter Initialization ===\n');
        p0 = interactivePointSelection(fixedGray, movingGray, p0);
        fprintf('Updated initial parameters from user input:\n');
        disp(p0);
    end

    % -------------------- Variable packing (bounds) -------------------
    if optScale
        x0     = p0(:).';
        lo     = lb(:).';
        hi     = ub(:).';
        toX    = @(y) invLogitToBox(y, lo, hi);
        toY    = @(x) boxToLogit(x, lo, hi);
        objY   = @(y) objMaskedMse(toX(y), movingGray, fixedGray, outSize);
        y0     = toY(x0);
    else
        sFix   = p0(4);
        x0s    = p0(1:3);
        lo     = lb(1:3);
        hi     = ub(1:3);
        toX    = @(y) invLogitToBox(y, lo, hi);
        toY    = @(x) boxToLogit(x, lo, hi);
        objY   = @(y) objMaskedMse([toX(y), sFix], ...
                                    movingGray, fixedGray, outSize);
        y0     = toY(x0s);
    end

    % -------------------- Unconstrained search (simplex) -------------
    opts = optimset('Display','iter', ...
                    'MaxFunEvals', 5000, ...
                    'MaxIter',     5000, ...
                    'TolX',        1e-14, ...
                    'TolFun',      1e-14);

    yOpt = fminsearch(objY, y0, opts);

    if optScale
        pOpt = toX(yOpt);
    else
        pOpt = [toX(yOpt), sFix];
    end

    % -------------------- Produce outputs ----------------------------
    % Output registered image with adaptive sizing
    registeredImg = WarpImageSimilarity(movingGray, pOpt);

    if doShow
        % For visualization, create version that matches fixed image size
        T = paramsToMatrix(pOpt);
        Tinv = inv(T);
        [X, Y] = meshgrid(1:outSize(2), 1:outSize(1));
        xs = X .* Tinv(1,1) + Y .* Tinv(2,1) + Tinv(3,1);
        ys = X .* Tinv(1,2) + Y .* Tinv(2,2) + Tinv(3,2);
        regViz = bilinearSampleFixed(movingGray, xs, ys);
        
        showFalsecolorOverlay(fixedGray, regViz);
        ttl = sprintf(['tx=%.3f  ty=%.3f  \\theta=%.3f^\\circ  ' ...
                       's=%.6f%s'], ...
                      pOpt(1), pOpt(2), pOpt(3), pOpt(4), ...
                      ternary(~optScale,' (fixed)',''));
        title(ttl);
    end
end

% ====================== Local helper functions =========================

function mseVal = objMaskedMse(p, movG, fixG, outSz)
% NaN-masked MSE between warped moving and fixed images.
% Handles different input image sizes by using overlap regions.
    T = paramsToMatrix(p);
    Tinv = inv(T);
    
    % Create fixed-size warped image for comparison with fixed image
    [X, Y] = meshgrid(1:outSz(2), 1:outSz(1));
    xs = X .* Tinv(1,1) + Y .* Tinv(2,1) + Tinv(3,1);
    ys = X .* Tinv(1,2) + Y .* Tinv(2,2) + Tinv(3,2);
    
    Iw = bilinearSampleFixed(movG, xs, ys);
    
    % Create mask for valid pixels in both images
    mask = ~isnan(Iw) & ~isnan(fixG) & isfinite(Iw) & isfinite(fixG);
    
    % Compute MSE only on overlap region
    d = Iw(mask) - fixG(mask);
    if isempty(d)
        mseVal = 1e6;  % discourage empty overlap
    else
        mseVal = mean(d.^2);
    end
end

function Z = bilinearSampleFixed(I, xs, ys)
% 2D bilinear interpolation with out-of-bounds = NaN for optimization.
[H, W] = size(I);
x0 = floor(xs); y0 = floor(ys);
x1 = x0 + 1; y1 = y0 + 1;
dx = xs - x0; dy = ys - y0;

% Check bounds for all four corners
m00 = x0>=1 & x0<=W & y0>=1 & y0<=H;
m10 = x1>=1 & x1<=W & y0>=1 & y0<=H;
m01 = x0>=1 & x0<=W & y1>=1 & y1<=H;
m11 = x1>=1 & x1<=W & y1>=1 & y1<=H;

% Initialize with NaN
Z = NaN(size(xs));

% Sample values at corners (0 where out of bounds)
V00 = zeros(size(xs)); V00(m00) = I(sub2ind([H W], y0(m00), x0(m00)));
V10 = zeros(size(xs)); V10(m10) = I(sub2ind([H W], y0(m10), x1(m10)));
V01 = zeros(size(xs)); V01(m01) = I(sub2ind([H W], y1(m01), x0(m01)));
V11 = zeros(size(xs)); V11(m11) = I(sub2ind([H W], y1(m11), x1(m11)));

% Bilinear weights
W00 = (1 - dx) .* (1 - dy);
W10 = dx .* (1 - dy);
W01 = (1 - dx) .* dy;
W11 = dx .* dy;

% Interpolate only where at least one corner is valid
val = W00.*V00 + W10.*V10 + W01.*V01 + W11.*V11;
anyMask = m00 | m10 | m01 | m11;
Z(anyMask) = val(anyMask);
end

function T = paramsToMatrix(p)
% 3x3 row-vector transform from [tx ty thetaDeg scale].
    tx = p(1);  ty = p(2);  th = p(3);  sc = p(4);
    c  = cosd(th);  s = sind(th);
    T  = [ sc*c   sc*s   0; ...
          -sc*s   sc*c   0; ...
            tx      ty   1 ];
end

function y = boxToLogit(x, lb, ub)
% Map x in (lb,ub) -> y in (-inf,inf) elementwise.
    epsv = 1e-9;
    z = (x - lb) ./ max(ub - lb, epsv);
    z = min(max(z, epsv), 1 - epsv);
    y = log(z ./ (1 - z));
end

function x = invLogitToBox(y, lb, ub)
% Map y in (-inf,inf) -> x in (lb,ub) elementwise.
    z = 1 ./ (1 + exp(-y));
    x = lb + (ub - lb) .* z;
end

function imgD = toDouble01(img)
% Convert numeric image to double in [0,1].
    imgD = double(img);
    if ~isfloat(img)
        mx = double(intmax(class(img)));
        if mx > 1
            imgD = imgD ./ mx;
        end
    end
    imgD = min(max(imgD, 0), 1);
end

function validateInputs(movingImg, fixedImg, p0, lb, ub)
% Strict checks for shapes, classes, and bounds.
    validateattributes(movingImg, {'numeric','logical'}, ...
        {'nonempty','2d'}, mfilename, 'movingImg', 1);
    validateattributes(fixedImg, {'numeric','logical'}, ...
        {'nonempty','2d'}, mfilename, 'fixedImg', 2);

    % Note: Different image sizes are now allowed for registration
    % The algorithm will find the optimal alignment using overlap regions
    
    validateattributes(p0, {'double'}, ...
        {'vector','numel',4,'finite','real'}, mfilename, 'p0', 3);
    validateattributes(lb, {'double'}, ...
        {'vector','numel',4,'finite','real'}, mfilename, 'lb', 4);
    validateattributes(ub, {'double'}, ...
        {'vector','numel',4,'finite','real'}, mfilename, 'ub', 5);

    if any(lb > ub)
        error('Each lb(i) must be <= ub(i).');
    end
end

function showFalsecolorOverlay(fixedGray, regGray)
% Simple falsecolor visualization handling different image sizes.
    fixG = fixedGray;
    regG = regGray;
    fixG(~isfinite(fixG)) = 0;
    regG(~isfinite(regG)) = 0;

    % Handle different sizes by cropping/padding to match fixed image
    [H_fix, W_fix] = size(fixG);
    [H_reg, W_reg] = size(regG);
    
    if H_fix == H_reg && W_fix == W_reg
        % Same size - direct overlay
        overlay = zeros([H_fix, W_fix, 3]);
        overlay(:,:,1) = fixG;        % red   = fixed
        overlay(:,:,2) = regG;        % green = registered
    else
        % Different sizes - crop/pad registered to match fixed
        regG_sized = zeros(H_fix, W_fix);
        h_end = min(H_fix, H_reg);
        w_end = min(W_fix, W_reg);
        
        if h_end > 0 && w_end > 0
            regG_sized(1:h_end, 1:w_end) = regG(1:h_end, 1:w_end);
        end
        
        overlay = zeros([H_fix, W_fix, 3]);
        overlay(:,:,1) = fixG;          % red   = fixed
        overlay(:,:,2) = regG_sized;    % green = registered (resized)
    end

    figure('Name','Registration Overlay','NumberTitle','off');
    imagesc(overlay);
    axis image off;
    
    % Add informative title
    if H_fix == H_reg && W_fix == W_reg
        title(sprintf('Registration Result: %dx%d\n(Red=Fixed, Green=Registered)', H_fix, W_fix));
    else
        title(sprintf('Registration Result: Fixed %dx%d, Registered %dx%d\n(Red=Fixed, Green=Registered)', ...
              H_fix, W_fix, H_reg, W_reg));
    end
end

function out = ternary(cond, a, b)
% Inline conditional string helper.
    if cond, out = a; else, out = b; end
end

function p0_updated = interactivePointSelection(fixedImg, movingImg, p0_default)
% Interactive point-based parameter initialization interface.
    fprintf('Interactive initialization: Select number of point pairs to use:\n');
    fprintf('  1 - Translation only (preserves rotation/scale from p0)\n');
    fprintf('  2 - Full similarity transform (translation, rotation, scale)\n');
    fprintf('  3+ - Least-squares similarity transform\n');
    
    % Get number of points from user
    numPairs = input('Enter number of point pairs (1-5): ');
    numPairs = max(1, min(5, round(numPairs))); % Clamp to valid range
    
    fprintf('\n--- Point Selection Process ---\n');
    fprintf('Step 1: Click %d point(s) in the FIXED (reference) image\n', numPairs);
    
    % Display fixed image for point selection
    fig1 = figure('Name','Select Points in FIXED Image','NumberTitle','off', ...
                  'Position',[100 300 600 500]);
    imagesc(fixedImg); axis image off; colormap gray;
    title(sprintf('Click %d point(s) in FIXED image, then press Enter', numPairs));
    [xf, yf] = ginput(numPairs);
    close(fig1);
    
    fprintf('Step 2: Click corresponding %d point(s) in the MOVING image\n', numPairs);
    
    % Display moving image for point selection
    fig2 = figure('Name','Select Corresponding Points in MOVING Image','NumberTitle','off', ...
                  'Position',[100 300 600 500]);
    imagesc(movingImg); axis image off; colormap gray;
    title(sprintf('Click %d corresponding point(s) in MOVING image, then press Enter', numPairs));
    [xm, ym] = ginput(numPairs);
    close(fig2);
    
    % Organize point data
    fixPts = [xf, yf];
    movPts = [xm, ym];
    
    % Compute initial parameters based on number of points
    try
        if numPairs == 1
            % Single point: translation only
            theta0 = p0_default(3); 
            s0 = p0_default(4);
            c = cosd(theta0); 
            s = sind(theta0);
            A = s0 * [c s; -s c];
            p0_updated = p0_default;
            p0_updated(1:2) = fixPts(1,:) - movPts(1,:) * A;
            fprintf('Using translation-only initialization (1 point pair)\n');
            fprintf('Preserved rotation: %.2f°, scale: %.3f\n', theta0, s0);
            
        elseif numPairs == 2
            % Two points: closed-form similarity
            p0_updated = estimateSimilarityTwoPoints(movPts(1,:), movPts(2,:), ...
                                                    fixPts(1,:), fixPts(2,:));
            fprintf('Using closed-form similarity initialization (2 point pairs)\n');
            
        else
            % Three or more points: least-squares similarity
            p0_updated = estimateSimilarityLeastSquares(movPts, fixPts);
            fprintf('Using least-squares similarity initialization (%d point pairs)\n', numPairs);
        end
        
        fprintf('Estimated parameters: [tx=%.2f, ty=%.2f, theta=%.2f°, scale=%.3f]\n', ...
                p0_updated(1), p0_updated(2), p0_updated(3), p0_updated(4));
                
    catch ME
        warning('Point-based initialization failed: %s');
        fprintf('Falling back to provided p0 parameters\n');
        p0_updated = p0_default;
    end
end

function p = estimateSimilarityTwoPoints(m1, m2, f1, f2)
% Closed-form similarity estimation from two point pairs.
    u = m2 - m1; 
    v = f2 - f1; 
    nu2 = dot(u, u);
    
    if nu2 < 1e-14
        error('Point pairs too close together for reliable estimation');
    end
    
    sc = (u(1)*v(1) + u(2)*v(2)) / nu2;  % s*cos(theta)
    ss = (u(1)*v(2) - u(2)*v(1)) / nu2;  % s*sin(theta)
    sVal = hypot(sc, ss);
    th = atan2d(ss, sc);
    
    A = [sc ss; -ss sc];
    t = f1 - m1 * A;
    p = [t, th, sVal];
end

function p = estimateSimilarityLeastSquares(movPts, fixPts)
% Least-squares similarity estimation from multiple point pairs.
    if size(movPts,1) ~= size(fixPts,1) || size(movPts,1) < 3
        error('Need at least 3 point pairs with matching counts');
    end
    
    % Center point sets
    muM = mean(movPts, 1); 
    muF = mean(fixPts, 1);
    Mc = movPts - muM;    
    Fc = fixPts - muF;
    
    % Cross-covariance matrix
    C = (Mc' * Fc) / size(movPts, 1);
    
    % SVD for optimal rotation
    [U, S, V] = svd(C);
    R = V * U';
    
    % Ensure proper rotation (det(R) = 1)
    if det(R) < 0
        V(:,2) = -V(:,2); 
        R = V * U';
    end
    
    % Estimate scale
    varM = trace((Mc' * Mc) / size(movPts, 1));
    sVal = trace(S) / max(varM, 1e-14);
    
    % Estimate translation
    t = muF - muM * (sVal * R);
    
    % Convert rotation matrix to angle
    th = atan2d(R(1,2), R(1,1));
    
    p = [t, th, sVal];
end