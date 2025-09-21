% demo_ImageRegistration. Demonstration of hyperspectral image registration 
% using sum-based grayscale alignment.
%
% REFERENCES:
% Demonstrates the imageRegistration function developed based on methods 
% from Sara Piqueras et al. "Handling Different Spatial Resolutions in 
% Image Fusion by Multivariate Curve Resolution–Alternating Least Squares 
% for Incomplete Image Multisets." Analytical Chemistry 2018, 90(11), 
% 6757–6765.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2025-08-02
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The demo_ImageRegistration script demonstrates the complete hyperspectral
% image registration workflow. It creates a realistic synthetic hyperspectral
% cube with 100 spectral bands using Gaussian-based pure spectra, applies a
% known geometric distortion, then recovers the alignment using sum-based
% registration and applies the correction to the entire hyperspectral dataset.
%
% The workflow follows the standard approach: HSI → sum to grayscale → 
% register → apply transformation to all bands. If images have clear
% contour, it is recommended binarize them before registration.
%
% INPUTS:
% None. All parameters are set internally for demonstration purposes.
%
% OUTPUTS:
% The function generates figure windows showing:
% - Sum images (fixed and moving) used for registration
% Console output includes:
% - Registration parameters and accuracy assessment
% - Processing progress for the 100-band hyperspectral cube
%
% EXAMPLE:
% % Run the hyperspectral registration demonstration:
% demo_ImageRegistration;
%
% DISCLAIMER:
% Authors and Lovelace's Square are not responsible for any issues,
% inaccuracies, or data loss arising from the use of this function.
clear, clc; close all; rng(0);

%% Generate synthetic hyperspectral data
H = 256; W = 256; nBands = 100;
[x, y] = meshgrid(1:W, 1:H);

fprintf('Creating %d-band hyperspectral cube...\n', nBands);

% Simple spatial patterns
circle = (x - W/4).^2 + (y - H/2).^2 <= (H/6)^2;
rectangle = x > 3*W/5 & x < 4*W/5 & y > H/4 & y < 3*H/4;
triangle = inpolygon(x, y, [W/2, 0.9*W, 0.1*W], [H/6, 5*H/6, 5*H/6]);

% Spectral signatures
wavelengths = linspace(400, 900, nBands);
veg_spectrum = 0.3 * exp(-((wavelengths - 550)/50).^2) + 0.8 * exp(-((wavelengths - 720)/80).^2);
soil_spectrum = 0.6 * exp(-((wavelengths - 650)/120).^2);
water_spectrum = 0.9 * exp(-((wavelengths - 450)/200).^2);

% Create hyperspectral cube
hyperspectral_original = zeros(H, W, nBands);
for band = 1:nBands
    hyperspectral_original(:,:,band) = ...
        double(circle) * veg_spectrum(band) + ...
        double(rectangle) * soil_spectrum(band) + ...
        double(triangle) * water_spectrum(band);
end

% Add noise
noise_level = 0.1;
hyperspectral_original = hyperspectral_original + noise_level * randn(size(hyperspectral_original));
hyperspectral_original = max(0, hyperspectral_original);

%% Apply distortion to create moving image
p_true = [60, -40, 20, 1.2];  % [tx ty thetaDeg scale]
fprintf('Applying distortion: tx=%.1f, ty=%.1f, theta=%.1f°, scale=%.3f\n', p_true);

% Apply distortion to create moving hyperspectral image (adaptive sizing)
hyperspectral_moving = WarpImageSimilarity(hyperspectral_original, p_true);
fprintf('Original size: %dx%dx%d\n', size(hyperspectral_original));
fprintf('Moving size: %dx%dx%d\n', size(hyperspectral_moving));
hyperspectral_moving=hyperspectral_moving(20:end-30,10:end-40,:);

%% Registration workflow
fprintf('\n=== Registration Workflow ===\n');

% Step 1: Create sum images
fprintf('Step 1: Creating sum images...\n');
fixed_sum = sum(hyperspectral_original, 3);
moving_sum = sum(hyperspectral_moving, 3);

% Normalize
fixed_sum = fixed_sum / max(fixed_sum(:));
moving_sum = moving_sum / max(moving_sum(:));

% Step 2: Register using sum images
fprintf('Step 2: Registering sum images...\n');
[~, pOpt] = imageRegistration(moving_sum, fixed_sum);

% Step 3: Apply transformation to hyperspectral cube
fprintf('Step 3: Applying transformation to hyperspectral cube...\n');
hyperspectral_aligned = WarpImageSimilarity(hyperspectral_moving, pOpt);
fprintf('Aligned size: %dx%dx%d\n', size(hyperspectral_aligned));

%% Results
p_inv_expected = invertSimilarity(p_true);
error_vals = abs(pOpt - p_inv_expected);
pct_errors = abs(error_vals ./ p_inv_expected) * 100;

fprintf('\n=== Results ===\n');
fprintf('Ground truth transform: [%.1f, %.1f, %.1f°, %.3f]\n', p_true);
fprintf('Expected recovery: [%.3f, %.3f, %.3f°, %.6f]\n', p_inv_expected);
fprintf('Recovered parameters: [%.3f, %.3f, %.3f°, %.6f]\n', pOpt);
fprintf('Percentage errors: [%.2f%%, %.2f%%, %.2f%%, %.2f%%]\n', pct_errors);

if all(pct_errors < 2.0)
    fprintf('✓ Registration successful!\n');
else
    fprintf('⚠ Registration accuracy could be improved\n');
end

%% Visualization
figure('Position', [100 100 1200 400]);

subplot(1,3,1);
imagesc(fixed_sum); axis image off; colormap gray;
title('Fixed (Original)');

subplot(1,3,2);
imagesc(moving_sum); axis image off; colormap gray;
title('Moving (Distorted)');

subplot(1,3,3);
aligned_sum = sum(hyperspectral_aligned, 3);
if ~isempty(aligned_sum)
    aligned_sum = aligned_sum / max(aligned_sum(:));
    
    % Create overlay
    [H_fix, W_fix] = size(fixed_sum);
    [H_align, W_align] = size(aligned_sum);
    
    if H_fix == H_align && W_fix == W_align
        overlay = zeros([H_fix, W_fix, 3]);
        overlay(:,:,1) = fixed_sum;    % Red
        overlay(:,:,2) = aligned_sum;  % Green
    else
        % Crop to common size
        min_H = min(H_fix, H_align);
        min_W = min(W_fix, W_align);
        overlay = zeros([min_H, min_W, 3]);
        overlay(:,:,1) = fixed_sum(1:min_H, 1:min_W);
        overlay(:,:,2) = aligned_sum(1:min_H, 1:min_W);
    end
    
    imagesc(overlay); axis image off;
    title('Registration Result');
end

fprintf('\nDemo completed!\n');


%% Helper function
function pinv = invertSimilarity(p)
    tx = p(1); ty = p(2); th = p(3); sc = p(4);
    c = cosd(th); s = sind(th);
    Ainv = (1/sc) * [c -s; s c];
    tInv = -[tx ty] * Ainv;
    pinv = [tInv, -th, 1/sc];
end