classdef ParticleDetector
% ParticleDetector. Automatic object detection for background cropping.
%
%   Detects particles, cells, or well-defined objects using intensity-based
%   thresholding and morphological cleanup. Requires the Image Processing
%   Toolbox for full functionality.
%
%   Methods (Static):
%     detect       - Run automatic detection, returns binary mask
%     checkToolbox - Verify Image Processing Toolbox is available
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-29
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.1

    methods (Static)

        function hasIPT = checkToolbox()
        % checkToolbox  Verify that the Image Processing Toolbox is available.
            hasIPT = ~isempty(ver('images'));
        end

        function [mask, numObjects] = detect(intensityImg, params)
        % detect  Run automatic particle/cell detection.
        %
        %   [mask, numObjects] = ParticleDetector.detect(intensityImg, params)
        %
        %   Input
        %     intensityImg  [rows x cols] 2-D intensity image
        %     params        struct with optional fields:
        %       .method       'otsu' | 'adaptive'         (default 'otsu')
        %       .sensitivity  adaptive threshold 0-1      (default 0.5)
        %       .minArea      min object area in pixels    (default 50)
        %       .watershed    logical, separate touching   (default false)
        %       .smoothSigma  Gaussian smoothing sigma     (default 1.0)
        %
        %   Output
        %     mask          [rows x cols] logical (true = retained)
        %     numObjects    number of detected objects

            if ~ParticleDetector.checkToolbox()
                error('ParticleDetector:noIPT', ...
                    'Image Processing Toolbox is required for auto-detection.');
            end

            % Defaults
            if ~isfield(params, 'method'),      params.method      = 'otsu'; end
            if ~isfield(params, 'sensitivity'),  params.sensitivity = 0.5;    end
            if ~isfield(params, 'minArea'),      params.minArea     = 50;     end
            if ~isfield(params, 'watershed'),    params.watershed   = false;  end
            if ~isfield(params, 'smoothSigma'),  params.smoothSigma = 1.0;    end

            % Normalise to [0, 1]
            img = mat2gray(double(intensityImg));

            % Gaussian smoothing
            if params.smoothSigma > 0
                img = imgaussfilt(img, params.smoothSigma);
            end

            % --- Thresholding ---
            switch lower(params.method)
                case 'otsu'
                    level = graythresh(img);
                    bw = imbinarize(img, level);

                case 'adaptive'
                    sens = max(0, min(1, params.sensitivity));
                    T = adaptthresh(img, sens, ...
                        'ForegroundPolarity', 'bright');
                    bw = imbinarize(img, T);

                otherwise
                    level = graythresh(img);
                    bw = imbinarize(img, level);
            end

            % --- Morphological cleanup ---
            bw = imopen(bw,  strel('disk', 2));
            bw = imclose(bw, strel('disk', 3));
            bw = imfill(bw, 'holes');
            bw = bwareaopen(bw, max(1, round(params.minArea)));
            bw = imclearborder(bw);

            % --- Watershed (optional, for touching particles) ---
            if params.watershed
                D = bwdist(~bw);
                D = imgaussfilt(D, 2);
                localMax  = imextendedmax(D, 2);
                D_mod     = imimposemin(-D, localMax);
                L         = watershed(D_mod);
                bw(L == 0) = false;
                bw = bwareaopen(bw, max(1, round(params.minArea)));
            end

            % --- Count objects ---
            CC = bwconncomp(bw);
            numObjects = CC.NumObjects;
            mask = bw;
        end

    end
end
