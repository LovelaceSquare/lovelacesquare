classdef LASLSCorrector
    % LASLSCorrector - Core algorithm for Local Asymmetric Least Squares (LAsLS) baseline correction.
    %
    % Solves via Iteratively Reweighted Least Squares (IRLS):
    %   (W + D' * diag(lambdaVec) * D + mu * L' * L) * b = W * y
    %
    % where W is the diagonal weight matrix, D and L are second- and
    % first-order sparse difference operators, lambdaVec contains
    % position-dependent smoothing penalties, and mu is the global
    % first-derivative penalty.
    %
    % REFERENCES:
    %   Eilers, Paul H.C., and Hans F.M. Boelens.
    %   "Baseline correction with asymmetric least squares smoothing."
    %   Leiden University Medical Centre Report 1.1 (2005): 5.
    %
    % Authors: Adrián Gómez-Sánchez, Berta Torres-Cobos, Rodrigo Rocha de Oliveira
    % Date Created: 2024-12-16
    % License: MIT
    % Repository: https://github.com/LovelaceSquare/lovelacesquare
    % Reviewed by Lovelace's Square: Yes
    % Version: 1.1

    methods (Access = public)
        function obj = LASLSCorrector()
            % Constructor - no initialization required
        end

        function [baseline, weights] = computeBaseline(~, y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol, globalP)
            % computeBaseline - Compute LALS baseline for a single spectrum.
            %
            % Inputs:
            %   y           - Column vector (n x 1) of spectral intensities
            %   intervals   - (m x 2) matrix of [startIdx, endIdx] per interval
            %   pVals       - (m x 1) vector of local p values per interval
            %   lambdasAsym - (m x 1) vector of local lambda values per interval
            %   lambdaWhit  - Scalar, global lambda (smoothness outside intervals)
            %   mu          - Scalar, first-derivative penalty weight
            %   maxIter     - Scalar, maximum IRLS iterations
            %   tol         - Scalar, convergence tolerance
            %   globalP     - Scalar, global asymmetry parameter (default p)
            %
            % Outputs:
            %   baseline    - (n x 1) estimated baseline
            %   weights     - (n x 1) final IRLS weights

            %% 1) Input checks and setup
            n = length(y);
            if size(y, 2) > 1
                y = y(:);
                n = length(y);
            end

            if isempty(intervals)
                m = 0;
            else
                m = size(intervals, 1);
                if length(pVals) ~= m || length(lambdasAsym) ~= m
                    error('LASLSCorrector:paramMismatch', ...
                        'Number of intervals must match length of pVals and lambdasAsym.');
                end
            end

            %% 2) Build per-index parameter vectors
            % pVec: asymmetry weight per data point (global default, overridden in intervals)
            pVec = globalP * ones(n, 1);
            % lambdaVec: smoothness penalty per difference row (global default, overridden in intervals)
            lambdaVec = lambdaWhit * ones(n - 2, 1);

            for j = 1:m
                iStart = intervals(j, 1);
                iEnd   = intervals(j, 2);
                thisP  = pVals(j);
                thisLam = lambdasAsym(j);

                if iStart < 1 || iEnd > n || iStart > iEnd
                    warning('LASLSCorrector:invalidInterval', ...
                        'Skipping invalid interval [%d, %d].', iStart, iEnd);
                    continue;
                end

                % Override p for data points within this interval
                pVec(iStart:iEnd) = thisP;

                % Override lambda for difference rows that overlap this interval
                % The second-order difference at row r involves indices r, r+1, r+2
                % so rows affected are max(1, iStart-1) to min(n-2, iEnd-1)
                rowStart = max(1, iStart - 1);
                rowEnd   = min(n - 2, iEnd - 1);
                for r = rowStart:rowEnd
                    lambdaVec(r) = thisLam;
                end
            end

            %% 3) Construct sparse difference operators
            % D: second-order difference operator (n-2 x n)
            % D = [1, -2, 1] pattern via spdiags
            e = ones(n, 1);
            D = spdiags([e, -2*e, e], 0:2, n - 2, n);

            % L: first-order difference operator (n-1 x n)
            L = spdiags([-ones(n, 1), ones(n, 1)], [0, 1], n - 1, n);

            %% 4) IRLS initialization
            weights = ones(n, 1);
            baseline = zeros(n, 1);

            %% 5) Iterative Reweighted Least Squares (IRLS)
            LambdaDiag = spdiags(lambdaVec, 0, n - 2, n - 2);
            DtLD = D' * LambdaDiag * D;
            LtL = mu * (L' * L);

            for iter = 1:maxIter
                W = spdiags(weights, 0, n, n);
                A = W + DtLD + LtL;
                rhs = W * y;
                newBaseline = A \ rhs;

                relChange = norm(newBaseline - baseline, 2) / (norm(baseline, 2) + eps);
                baseline = newBaseline;

                if relChange < tol
                    break;
                end

                % Update weights based on asymmetry
                for i = 1:n
                    if y(i) > baseline(i)
                        weights(i) = pVec(i);
                    else
                        weights(i) = 1 - pVec(i);
                    end
                end
            end
        end

        function [corrected, baseline, weights] = correctMatrix(obj, data, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol, globalP)
            % correctMatrix - Apply LALS baseline correction to each row of a matrix.
            %
            % Inputs:
            %   data        - (nRows x nCols) spectral data matrix
            %   intervals   - (m x 2) matrix of [startIdx, endIdx] per interval
            %   pVals       - (m x 1) vector of local p values per interval
            %   lambdasAsym - (m x 1) vector of local lambda values per interval
            %   lambdaWhit  - Scalar, global lambda
            %   mu          - Scalar, first-derivative penalty weight
            %   maxIter     - Scalar, maximum IRLS iterations
            %   tol         - Scalar, convergence tolerance
            %   globalP     - Scalar, global asymmetry parameter
            %
            % Outputs:
            %   corrected   - (nRows x nCols) baseline-corrected data
            %   baseline    - (nRows x nCols) estimated baselines
            %   weights     - (nRows x nCols) final IRLS weights

            [nRows, nCols] = size(data);
            corrected = zeros(nRows, nCols);
            baseline  = zeros(nRows, nCols);
            weights   = zeros(nRows, nCols);

            for i = 1:nRows
                y = data(i, :)';
                [bl, w] = obj.computeBaseline(y, intervals, pVals, ...
                    lambdasAsym, lambdaWhit, mu, maxIter, tol, globalP);
                baseline(i, :)  = bl';
                weights(i, :)   = w';
                corrected(i, :) = data(i, :) - bl';
            end
        end
    end
end

%
%
%                                               AAA       AAA
%                                            AA  AAA      AA
%                                          AA  AAA       A
%                              AAAAAA     A  AA  AAA    A
%                            A   A  AA  AA  AAA        A
%                            A   A  A  AA AAAAAA      A
%                            AAA A  A  AAAAAAAAAAA  AA
%                            A AAA  A  A AAAA      A
%                            A  A A A  A AA     AA      AAAAAAAA
%                            A  A A A  AA     AA   AAA     AAAA   AAA
%                             A  AAAAAAAAAAAA  AA   AAAAAAAAAAAAAAAA  A
%                             AA  AAAA       A   AAAAAAAAAAAAAAAAAAAAA  A
%                              A  AAAAA   AA  AAAAAAAAAAAAAAAAAAAAAAAAAA A
%                               AAAAAAAAAA  AAAAAAAAAAAAAAAAAAAAAAAAAAAAA A
%                                AAA AAA  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA A
%                                AAA AA   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%                         AAAA        AAA AAAAAAAAAAAAAAAAAAAAAA    AAAAAAAA
%                     AAA     AAAAAAAAAAAAAAAAAAAAAAAAA   AAAAAAAAAAAAAA  A
%                   A    AAAAAAAAAAAAAAAAAAAAAAAA  AAAAAAAAAAAAAAAAAAAAAA
%                AA   AAAAAAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAAAAA
%              AA   AAAAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAAA AAA                              A
%             A   A   AAAAAAAAAAAAAAAA AAA  AAAAAAAAAAAAAAAAAAAAAA                               AA
%            A  AA    AAAAAAAAAAAAA AAA     AAAAAAAAAAAAAA AAAA                               AA A
%           A  AA    AAAAAAAAAAA  AAAAAA     AAAAAAAAAAAAAA                                 AA  A
%          A  AAAA AAAAAAAAAAA AAAAAAAAAA     AAAAAAAA                                    AA   A     A
%         A  AAAAAAAAAAAAAA  AAAAAAAAAAAAA     A               AAAAAAA      AAAA        A     A    AAA
%        A  AAAAAAAAAAAAA  AAAAAAAAAAAAAAA     A         AAAA                   AA    A   A  A  AA  A
%        A  AAAAAAAAAAAA AAAAAAAAAAAAAAAA A     A    AAA AAA                      AAA    AA AAAA   A           AAA
%       AA AAAAAAAAAAA AAAAAAAAAAAAA AA    A    AAA   AA                           A    AAAAAA    A       AAA AA
%       A  AAAAAAAAAAAAAAAAAAAAAAAAAA       AAAA    A      AAA                     AA AAAA       A   AAAA    A
%       A  AAAAAAAA AAAAAAAAAAAAAAA       AA      AA    AA  AAA                     A  AA       AAAA       AA
%       AA AAAAAAAAAAAAAAAAAAAAA       AA       AA    A  AA                 A  A    A AA    A            AA    A
%        A  AAAAAAAAAAAA   AAA       AA AAA      A      A                     A A   AA   AAA            AAAAAAA
%        AA AAAAAAAAA   AAA       AA AA     AAAAA    A A       AAA     AAA     A AA  A AAA     AA         AA
%         A  AAAAA   AA         AA AA   AA   AAAA  A A      AA             AA   A AA AAAA      AA      AA      AAAAAA
%          AAAAAAAA            A  A  AA   AA      A A A   AAA           A     A A AAAA A     AAAAA   AAAA         AA
%                            A   A  A   AA          A  AA     AA    A    A     AAAAAAA     AAAAAA               AA
%                          AA   A AA   A            A       A A  AA  AA   A AA  AAAAAA   AA AAA               AA
%                         A    A  A   A     AAA  AAA  A  AAAA A A  AA  A   A AA AAAAAA  A                   AA
%                        A           A   AA     AAA  AAAAAAA  AA      AAA   A  AAAAAAA       A      AAA   AA
%                      AA              A       A   AAAAAAA    AA        AAA  A  AAAAA     AAAA AAA AAAA   AAAAAAAAAAAAAAAAAAAA
%                     AA    A    A   AA       A  AAAAA        A       AA   AAAAA   AA AAAAAAA   AAAAA                    AA
%                    AA       A  A  A      AAA   AAAAAA       A        AA      AAAAAA  AAAAA     AA                  AAA
%                   AAA      AA A  AA    A AA       AA A      A         A          AA   AA     AA                AAA
%        AAAAAAAA     A      A  A  A     A AA   AA      A     A         AA        AA     AAAAAAAAAAAA         AA
%         AAA      A   A     A  A A      A A      A A    AA   A                   AA           AAAAA      AA
%             AAA       AA   A  A A    A AA       AAA     AA   A                   AA       AA AAA        AAA
%         AAAA           AA   A A A    AAAA                AA  A              A A   A       AAA                AAA
%           AA       AAAA AA     AA    AAA                  A   A             AAAAAAA   AAA                        AA
%              AA       AA AAA   AAA   AAA                   A  A    AA    A      AAAA                      AAAA
%         AAAAA      AAAAA  A AA  AA   AAA                    A  A   AAA   AAA    AAAA    A   AA    AAAAA
%  AAAA           AAAAAAAAA     A AAA  AA                      A A  AAA     AAA  AAAA A  A AA  AA      AA
% AAA                    AA A     A  A AA               A  A    A A AAA  A  AAA AAAAA A AAAA    A AAA    AA
%      AAAA             A  A  A     AAAAA             AAAA  A    AAA  A AA  AAA    A   AAAAAAA  AAAAA       AA
%            AAAA     AA AAAAA        AA   A         AAA      A  AAA   AA AA A          AA AA A A  AAAAAAAAA
%              AA    AAA AAAA   A       A   A A     AAAA       AA AAA  AAAAAAAA   A     AA  AAA AAA  AAA
%           AA    AA   AA AA  AA  A AAA  AA AAA      AAA  A AAA    AAA   A  A  A A       AAA   A
%        AA          A   AA  A   A AAA   AAA          AA AAA    A  AAAA  A  A  AAA   AA    A AA
%       AAAA    AA       A AA  AA AAAA AA   A            AA  AAAA   A AA AAAAAA AA  AAAAA   A
%              A       AA A  AAA  A A AA     A      A    AAAA     A  A AA AAAAA     AAA  A   AA
%            A        AAAAAAAAA  A A A   A A AA    A AA   AAA AAA AA A   A   AAAAAA     AA    AA
%          AA          AA  AAA  A AAA   AAA A  A    AAA AA  A    AAA  A   AAAAAAAA A    A  A    A
%        AA               A      A  A  AAAAA    A         A  AAAA A   AA AAAAAAA   A    A AAA    A
%      AA          AAA           A  A AA AA AAAAAA           AAAA  AAAA AAAAAAA    AAAAAAA AAAAA  A
%     A       AAA   A           A  A AAA A A A   AA              AAAAAAAAAAA       A        AAAAAAAA
%   A   AAAA      AA       AA  A  A  AAAA AA A AA  A     AAAA   AA AAA            A            AA       AAAAAAA
% AAAAA          A     AAA  A        AAA   A AAAAA AA    AAAAA  A                A    AAAAAA AAA  AAAAAAA        AAAAAAAAA
%              AA   AAA    A       AA  A AAAA   AAAAA     AA    AA              A       AAAAAAAAAAAA AAAAAAAAA
%             A AAA        A      A   AA AA  AAA     A         AA AAA        AA      A  A   AAAAAA   AA
%           AAA            A    AA    A A  AA     A AA   AAAAAAA  A   AAAAA        AA AAAAA           A
%                         A    A      A  AA          A         AAA                AAAAAAAAA            A
%                         A  A       A A             AA  A AAA  A      A AA AAAA  AAAAAAAAA  AAA     A A
%                         AAA        A                A  A A AAA      A AAA A  AAAAAAAA  AAA  AA     A A
%                        AA                           A   A AAA    A  A AAAAAAAAAAA AA    AA   AA    AAA
%                                                     AA  AAAA    A   A AAAAAAAAAA A A     AA  AA  AAA A
%                                                      A   AA    AAA AAAAAAAAAAA  A         AA AA  AAAAA
%                                                       A A AA   A AA   AAAAAAA   A          AAAA  AAAAA
%                                                       AAA A A  AA   A  AA        A          AA   AAAAA
%                                                       AA A  AAA     A  AAA       AA          A   AAAAA
%                                                      AA  A AA     AAAAA AA        A           A  AAAA
%                                                     AAA AAA     AA   AAAAA        AA           A  AAA
%                                                   AA A AA  A     AA   AAA          A           AAA A
%                                                  AA AAA     AA    AA    A          AA           A A
%                                                 A  A          A     A  AA           A            A
%                                                AAA              A    AAA            A            AA
%                                               AA                  AA  A             AA            A
%                                                                      A               A             A
%                                                                     A                A             A
%                                                                    A                 A             A
%                                                                  AA                  A              A
%                                                                 AA                   AA             A
%                                                                AA                    AA             A
%                                                               AA                     A              A
%                                                               A                      A              A
%                                                              A                       A              AA
%                                                              A                      AA               A
%                                                             A                                       A
%                                                             A                      A                A
%                                                             A                    A                  A
%                                                             A                                      AA
%                                                             A                                      A
%                                                              A                                    A
%                                                              A                                   A
%                                                               A                                 A
%                                                               AAA                              A
%                                                                AAA                        A  AA
%                                                                  AA                 A  A AAAA
%                                                                   AAAAAAAAAA AA AA A  A AAA
%                                                                     AAAAA AAA AA AA AAAA
%                                                                         AAAAA AAAAAAA
%
