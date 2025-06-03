% test.m
%
% This script loads the default 'cameraman.tif' image, applies binning to
% reduce its resolution using both 'sum' and 'mean' methods, and visualizes
% the results for a clear comparison.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2025-04-04
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0

%% Load Default Image
% Read the default image and convert it to double for processing.
originalImage = imread('cameraman.tif');
originalImage = double(originalImage);

%% Display Original Image
figure;
subplot(1,3,1);
imagesc(originalImage);
colormap(gray);
axis image off;
title('Original Image');

%% Apply Binning
% Define a bin vector to group pixels in 2x2 blocks.
binVector = [7, 7];

% Apply binning using 'sum' mode, which sums the pixel values in each block.
binnedSum = binning(originalImage, binVector, 'sum');

% Apply binning using 'mean' mode, which averages the pixel values in each block.
binnedMean = binning(originalImage, binVector, 'mean');

%% Display Binned Images
subplot(1,3,2);
imagesc(binnedSum);
colormap(gray);
axis image off;
title('Binned Image (Sum Mode)');

subplot(1,3,3);
imagesc(binnedMean);
colormap(gray);
axis image off;
title('Binned Image (Mean Mode)');

%% Completion Message
disp('Binning image test completed. Please review the displayed images.');
