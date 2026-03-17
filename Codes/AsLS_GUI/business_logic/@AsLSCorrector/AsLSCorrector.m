classdef AsLSCorrector
% AsLSCorrector. Core AsLS operations used by the GUI backend.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Implements the Asymmetric Least Squares Smoothing algorithm (Eilers &
% Boelens, 2005) for baseline estimation. Solves the penalized system
% (I + lambda * D' * D) * z = y with iteratively reweighted asymmetric
% weights. Supports preview on mean spectrum and batch correction.
%
% METHODS:
%   previewBaseline - Compute baseline for the mean spectrum (real-time)
%   correct         - Apply baseline correction to all spectra
%
% EXAMPLE:
%   corrector = AsLSCorrector();
%   [meanSpec, baseline] = corrector.previewBaseline(data, 1e6, 1e-3, struct.empty, false, 1e5, 10);
%   [corrected, baselines] = corrector.correct(data, 1e6, 1e-3, struct.empty, false, 1e5, 10);
%
% Disclaimer:
%   Author and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this software.

    methods (Access = public)
        function [meanSpectrum, baseline] = previewBaseline(obj, data, lambda, p, intervals, smoothEnabled, smoothLambda, maxIter)
            arguments
                obj
                data (:,:) double {mustBeFinite}
                lambda (1,1) double {mustBePositive}
                p (1,1) double {mustBeGreaterThanOrEqual(p, 0), mustBeLessThanOrEqual(p, 1)}
                intervals struct
                smoothEnabled (1,1) logical
                smoothLambda (1,1) double {mustBePositive}
                maxIter (1,1) double {mustBeInteger, mustBePositive}
            end

            nCols = size(data, 2);
            y = mean(data, 1)';

            globalBaseline = aslsBaseline(obj, y, lambda, p, maxIter);

            if isempty(intervals)
                merged = globalBaseline';
            else
                localBaselines = nan(numel(intervals), nCols);
                for i = 1:numel(intervals)
                    idx = intervalToIndices(obj, intervals(i).start, intervals(i).end, nCols);
                    if isempty(idx)
                        continue;
                    end
                    zLocal = aslsBaseline(obj, y(idx), intervals(i).lambda, intervals(i).p, maxIter);
                    tmp = nan(1, nCols);
                    tmp(idx) = zLocal(:)';
                    localBaselines(i, :) = tmp;
                end
                merged = mergeBaselines(obj, globalBaseline', localBaselines);
            end

            if smoothEnabled
                merged = aslsBaseline(obj, merged', smoothLambda, 0.5, maxIter)';
            end

            meanSpectrum = y';
            baseline = merged;
        end

        function [correctedData, baselineData] = correct(obj, data, lambda, p, intervals, smoothEnabled, smoothLambda, maxIter)
            arguments
                obj
                data (:,:) double {mustBeFinite}
                lambda (1,1) double {mustBePositive}
                p (1,1) double {mustBeGreaterThanOrEqual(p, 0), mustBeLessThanOrEqual(p, 1)}
                intervals struct
                smoothEnabled (1,1) logical
                smoothLambda (1,1) double {mustBePositive}
                maxIter (1,1) double {mustBeInteger, mustBePositive}
            end

            [nRows, nCols] = size(data);
            baselineData = zeros(nRows, nCols);

            % Global baseline for each spectrum.
            for r = 1:nRows
                z = aslsBaseline(obj, data(r, :)', lambda, p, maxIter);
                baselineData(r, :) = z';
            end

            if ~isempty(intervals)
                localBaselines = nan(numel(intervals), nRows, nCols);

                for i = 1:numel(intervals)
                    idx = intervalToIndices(obj, intervals(i).start, intervals(i).end, nCols);
                    if isempty(idx)
                        continue;
                    end

                    for r = 1:nRows
                        ySub = data(r, idx)';
                        zLocal = aslsBaseline(obj, ySub, intervals(i).lambda, intervals(i).p, maxIter);
                        localBaselines(i, r, idx) = zLocal(:)';
                    end
                end

                merged = baselineData;
                for r = 1:nRows
                    rowLocals = squeeze(localBaselines(:, r, :));
                    merged(r, :) = mergeBaselines(obj, baselineData(r, :), rowLocals);
                end
                baselineData = merged;
            end

            if smoothEnabled
                for r = 1:nRows
                    smoothed = aslsBaseline(obj, baselineData(r, :)', smoothLambda, 0.5, maxIter);
                    baselineData(r, :) = smoothed';
                end
            end

            correctedData = data - baselineData;
        end
    end

    methods (Access = private)
        function idx = intervalToIndices(~, startValue, endValue, nCols)
            s = max(1, min(nCols, round(startValue)));
            e = max(1, min(nCols, round(endValue)));
            if e < s
                t = s;
                s = e;
                e = t;
            end
            if e <= s
                idx = [];
                return;
            end
            idx = s:e;
        end

        function rowMerged = mergeBaselines(~, rowGlobal, rowLocals)
            if isempty(rowLocals)
                rowMerged = rowGlobal;
                return;
            end

            if isvector(rowLocals)
                rowLocals = rowLocals(:)';
            end

            localAverage = mean(rowLocals, 1, 'omitnan');
            rowMerged = localAverage;
            noLocal = isnan(localAverage);
            rowMerged(noLocal) = rowGlobal(noLocal);
        end

        function z = aslsBaseline(~, y, lambda, p, maxIter)
            y = y(:);
            n = numel(y);
            if n < 3
                z = y;
                return;
            end

            d2 = diff(speye(n), 2);
            smoothTerm = d2' * d2;
            w = ones(n, 1);

            for k = 1:maxIter
                wDiag = spdiags(w, 0, n, n);
                c = wDiag + lambda * smoothTerm;
                z = c \ (wDiag * y);

                above = y > z;
                belowOrEqual = ~above;
                w = p * above + (1 - p) * belowOrEqual;
            end
        end
    end
end
