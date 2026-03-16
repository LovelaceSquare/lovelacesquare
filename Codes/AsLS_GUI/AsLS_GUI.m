classdef AsLS_GUI < matlab.apps.AppBase
% AsLS_GUI. Interactive baseline correction using Asymmetric Least Squares.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Estimates and removes baselines from spectral datasets using the AsLS
% (Asymmetric Least Squares Smoothing) algorithm. The method fits a smooth
% baseline below the signal using penalized least squares with asymmetric
% weights, controlled by two parameters: lambda (smoothness) and p
% (asymmetry weight).
%
% The GUI provides real-time baseline preview as sliders are dragged, with
% all spectra overlaid, signal-by-signal browsing, and before/after views.
%
% Uses the AppBase + uihtml architecture: HTML/CSS/JS frontend with a
% MATLAB backend. Charts rendered with Canvas 2D (no MATLAB axes).
%
% ALGORITHM:
%   Solves (I + lambda * D' * D) * z = y with asymmetric weights, where:
%   - I is the identity matrix
%   - D is a difference matrix of order d
%   - y is the input signal
%   - Weights are updated iteratively: w_i = p if y_i > z_i, else (1 - p)
%
% PARAMETERS:
%   lambda   - Smoothness penalty (default: 1e6, range: 10^3 -- 10^8)
%   p        - Asymmetry weight (default: 1e-3, range: 10^-4 -- 10^-1)
%   maxIter  - Maximum iterations (default: 10)
%
% GUI FEATURES:
%   - Real-time baseline preview as sliders are dragged
%   - Logarithmic sliders for lambda and p with manual numeric input
%   - Signal-by-signal navigation with bold highlight
%   - Show all spectra overlay with subsampling for large datasets
%   - Batch correction across all spectra with progress bar
%   - Zoom (mouse wheel) and pan (click-drag) on both charts
%   - Large dataset guard with subsampling options (> 1000 samples)
%   - Export corrected data and baselines to MATLAB workspace
%   - Dark mode toggle (Ctrl+D)
%   - Resizable control/chart panels via drag handle
%   - Demo data with synthetic spectra and polynomial baselines
%
% USAGE:
%   app = AsLS_GUI();              % Empty (load from GUI)
%   app = AsLS_GUI(spectra);       % Direct numeric matrix
%   app = AsLS_GUI(s);             % Struct with .data, .wavelength
%
% INPUTS:
%   spectra    - Numeric matrix [nSamples x nChannels]
%   s          - Struct with fields:
%                .data       - Numeric matrix [nSamples x nChannels]
%                .wavelength - 1D vector [1 x nChannels] (optional)
%                .varName    - Character vector (optional)
%
% PUBLIC API:
%   app.setInputData(s)       - Load data programmatically
%   app.getData()             - Retrieve corrected data (struct with
%                               .data, .baselines, .wavelength)
%   app.waitForOrchestrator() - Block until user clicks Next or closes
%
% EXAMPLE:
%   % Generate demo data and launch GUI
%   AsLS_test;
%   s.data = spectra; s.wavelength = wavelength;
%   app = AsLS_GUI(s);
%
%   % Use in a pipeline (orchestrator mode)
%   app = AsLS_GUI(s);
%   app.waitForOrchestrator();
%   result = app.getData();
%   correctedSpectra = result.data;
%
% REFERENCES:
%   Eilers, P. H. C. (2003). A Perfect Smoother. Analytical Chemistry,
%   75(14), 3631-3636.
%
%   Eilers, P. H. C. & Boelens, H. F. M. (2005). Baseline Correction
%   with Asymmetric Least Squares Smoothing.
%
% See also: AsLSCorrector, DataValidator
%
% Disclaimer:
%   Author and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this software.

    properties (Access = public)
        UIFigure        matlab.ui.Figure
        HTMLComponent   matlab.ui.control.HTML
    end

    properties (Access = private)
        Corrector                       % AsLSCorrector instance
        Validator                       % DataValidator instance
        OriginalData    double = []
        CorrectedData   double = []
        BaselineData    double = []
        Wavelength      double = []
        LoadedVarName   char = ''
        DataLoaded      logical = false
        IsClosed        logical = false
        IsWaiting       logical = false
        NextClicked     logical = false
        ResponseCounter double = 0
        PendingPayload  struct = struct.empty
    end

    methods (Access = public)
        function app = AsLS_GUI(inputData)
            createComponents(app);
            initializeBusinessLogic(app);
            registerApp(app, app.UIFigure);

            if nargin >= 1
                try
                    setInputData(app, inputData);
                catch ME
                    warning('AsLS_GUI:inputError', ...
                        'Could not load input data: %s', ME.message);
                end
            end

            runStartupFcn(app, @startupFcn);

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            app.IsClosed = true;
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)
        function createComponents(app)
            screenSize = get(0, 'ScreenSize');
            figWidth = min(1400, round(screenSize(3) * 0.8));
            figHeight = min(900, round(screenSize(4) * 0.85));

            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 figWidth figHeight];
            app.UIFigure.Name = 'AsLS GUI — Baseline Correction';
            app.UIFigure.Color = [0.91 0.92 0.93];
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @closeRequest, true);
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @figureSizeChanged, true);

            modulePath = fileparts(mfilename('fullpath'));
            htmlPath = fullfile(modulePath, 'ui', 'asls_baseline_correction_ui.html');

            app.HTMLComponent = uihtml(app.UIFigure);
            app.HTMLComponent.Position = [1 1 figWidth figHeight];
            app.HTMLComponent.HTMLSource = htmlPath;
            app.HTMLComponent.DataChangedFcn = createCallbackFcn(app, @HTMLDataChanged, true);

            app.UIFigure.Visible = 'on';
        end

        function initializeBusinessLogic(app)
            modulePath = fileparts(mfilename('fullpath'));
            addpath(fullfile(modulePath, 'business_logic'));
            app.Corrector = AsLSCorrector();
            app.Validator = DataValidator();
        end

        function startupFcn(app)
            movegui(app.UIFigure, 'center');
            sendResponse(app, 'ready', struct('message', 'Ready'));
        end

        function figureSizeChanged(app, ~)
            if isempty(app.UIFigure) || ~isvalid(app.UIFigure), return; end
            pos = app.UIFigure.Position;
            if ~isempty(app.HTMLComponent) && isvalid(app.HTMLComponent)
                app.HTMLComponent.Position = [1 1 pos(3) pos(4)];
            end
        end

        function closeRequest(app, ~)
            app.IsClosed = true;
            if app.IsWaiting
                try uiresume(app.UIFigure); catch; end
            end
            delete(app);
        end
    end

    % =====================================================================
    %  ACTION DISPATCHER
    % =====================================================================
    methods (Access = private)
        function HTMLDataChanged(app, ~)
            try
                data = app.HTMLComponent.Data;
                if ~isstruct(data) || ~isfield(data, 'action'), return; end
                if isfield(data, 'source') && strcmp(data.source, 'matlab'), return; end

                action = char(string(data.action));
                switch action
                    case 'uiReady',             handleUiReady(app);
                    case 'loadData',            handleLoadData(app);
                    case 'loadVariable',        handleLoadVariable(app, data);
                    case {'preview', 'previewBaseline'}, handlePreview(app, data);
                    case 'apply',               handleApply(app, data);
                    case 'prepareExport',       handlePrepareExport(app);
                    case 'doExport',            handleDoExport(app, data);
                    case 'getSpectrum',         handleGetSpectrum(app, data);
                    case 'createDemoData',       handleCreateDemoData(app);
                    case 'next',                handleNext(app);
                end
            catch ME
                sendResponse(app, 'error', struct('message', ME.message));
            end
        end

        function handleLoadData(app)
            try
                variables = collectWorkspaceMatrices(app);
                vectors = collectWorkspaceVectors(app);
                if isempty(variables)
                    sendResponse(app, 'error', struct('message', 'No numeric 2D matrices found in workspace.'));
                    return;
                end
                sendResponse(app, 'showVarList', struct('variables', {variables}, 'vectors', {vectors}));
            catch ME
                sendResponse(app, 'error', struct('message', ['Could not list variables: ' ME.message]));
            end
        end

        function handleLoadVariable(app, data)
            if ~isfield(data, 'varName')
                sendResponse(app, 'error', struct('message', 'Missing variable name.'));
                return;
            end

            varName = char(string(data.varName));
            if ~isvarname(varName)
                sendResponse(app, 'error', struct('message', sprintf('Invalid variable name: %s', varName)));
                return;
            end
            if ~evalin('base', sprintf('exist(''%s'', ''var'')', varName))
                sendResponse(app, 'error', struct('message', sprintf('Variable "%s" not found.', varName)));
                return;
            end

            sendResponse(app, 'progress', struct('title', 'Loading variable...', 'percent', 5));
            drawnow limitrate;

            raw = evalin('base', varName);
            if ~isreal(raw), raw = real(raw); end

            sendResponse(app, 'progress', struct('title', 'Validating data...', 'percent', 15));
            drawnow limitrate;

            [ok, msg] = app.Validator.validateData(raw);
            if ~ok
                sendResponse(app, 'error', struct('message', msg));
                return;
            end

            nChannels = size(raw, 2);
            lambda = safeDouble(app, data, 'lambda', 1e6);
            p = safeDouble(app, data, 'p', 1e-3);
            maxIter = max(1, round(safeDouble(app, data, 'maxIter', 10)));

            % Load x-axis vector if provided
            waveAxis = 1:nChannels;
            if isfield(data, 'xAxisVar') && ~isempty(data.xAxisVar)
                xVarName = char(string(data.xAxisVar));
                if isvarname(xVarName) && evalin('base', sprintf('exist(''%s'', ''var'')', xVarName))
                    xRaw = evalin('base', xVarName);
                    xRaw = double(xRaw(:))';
                    if numel(xRaw) == nChannels
                        waveAxis = xRaw;
                    end
                end
            end

            loadDataIntoState(app, double(raw), waveAxis, varName);

            sendResponse(app, 'progress', struct('title', 'Computing baseline preview...', 'percent', 35));
            drawnow limitrate;

            [meanSpectrum, previewBaseline] = app.Corrector.previewBaseline( ...
                app.OriginalData, lambda, p, struct.empty, false, 1e5, maxIter);

            nR = size(raw, 1);
            payload = struct();
            payload.message = sprintf('Loaded "%s" (%d x %d).', varName, nR, nChannels);
            payload.varName = varName;
            payload.nSamples = nR;
            payload.nChannels = nChannels;
            payload.wavelength = app.Wavelength;
            payload.meanSpectrum = meanSpectrum;
            payload.previewBaseline = previewBaseline;
            payload.lambda = lambda;
            payload.p = p;
            payload.maxIter = maxIter;

            if nR * nChannels <= 2000000
                sendResponse(app, 'progress', struct('title', 'Preparing spectra data...', 'percent', 70));
                drawnow limitrate;
                allSpec = cell(1, nR);
                specProgress = max(1, floor(nR / 10));
                for i = 1:nR
                    allSpec{i} = app.OriginalData(i, :);
                    if mod(i, specProgress) == 0
                        pct = 70 + round(25 * i / nR);
                        sendResponse(app, 'progress', struct('title', 'Preparing spectra data...', 'percent', pct, 'current', i, 'total', nR));
                        drawnow limitrate;
                    end
                end
                payload.allSpectra = allSpec;
            end

            sendResponse(app, 'dataLoaded', payload);
        end

        function handlePreview(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'error', struct('message', 'Load data first.'));
                return;
            end

            lambda = safeDouble(app, data, 'lambda', 1e6);
            p = safeDouble(app, data, 'p', 1e-3);
            maxIter = max(1, round(safeDouble(app, data, 'maxIter', 10)));

            signalIndex = round(safeDouble(app, data, 'signalIndex', 0));
            if signalIndex > 0 && signalIndex <= size(app.OriginalData, 1)
                dataForPreview = app.OriginalData(signalIndex, :);
            else
                dataForPreview = app.OriginalData;
            end

            [meanSpectrum, previewBaseline] = app.Corrector.previewBaseline( ...
                dataForPreview, lambda, p, struct.empty, false, 1e5, maxIter);

            payload = struct();
            payload.wavelength = app.Wavelength;
            payload.meanSpectrum = meanSpectrum;
            payload.previewBaseline = previewBaseline;
            payload.lambda = lambda;
            payload.p = p;
            payload.maxIter = maxIter;

            sendResponse(app, 'previewUpdated', payload);
        end

        function handleGetSpectrum(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'error', struct('message', 'Load data first.'));
                return;
            end
            idx = max(1, min(size(app.OriginalData, 1), round(safeDouble(app, data, 'signalIndex', 1))));
            payload = struct();
            payload.signalIndex = idx;
            payload.spectrum = app.OriginalData(idx, :);
            payload.wavelength = app.Wavelength;
            payload.nSamples = size(app.OriginalData, 1);
            if ~isempty(app.CorrectedData)
                payload.correctedSpectrum = app.CorrectedData(idx, :);
            end
            if ~isempty(app.BaselineData)
                payload.baseline = app.BaselineData(idx, :);
            end
            sendResponse(app, 'spectrumLoaded', payload);
        end

        function handleApply(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'error', struct('message', 'Load data first.'));
                return;
            end

            lambda = safeDouble(app, data, 'lambda', 1e6);
            p = safeDouble(app, data, 'p', 1e-3);
            maxIter = max(1, round(safeDouble(app, data, 'maxIter', 10)));

            nR = size(app.OriginalData, 1);
            nC = size(app.OriginalData, 2);
            signalIndex = round(safeDouble(app, data, 'signalIndex', 0));
            singleMode = signalIndex > 0 && signalIndex <= nR;

            % Initialize matrices if empty (first correction)
            if isempty(app.CorrectedData) || size(app.CorrectedData, 1) ~= nR
                app.CorrectedData = app.OriginalData;
                app.BaselineData = zeros(nR, nC);
            end

            if singleMode
                % Signal-by-signal: correct only the current spectrum
                sendResponse(app, 'progress', struct('title', sprintf('Correcting signal %d...', signalIndex), 'current', 0, 'total', 1));
                drawnow limitrate;

                [corrRow, blRow] = app.Corrector.correct( ...
                    app.OriginalData(signalIndex, :), lambda, p, struct.empty, false, 1e5, maxIter);
                app.CorrectedData(signalIndex, :) = corrRow;
                app.BaselineData(signalIndex, :) = blRow;

                sendResponse(app, 'progress', struct('title', 'Done', 'current', 1, 'total', 1));
                drawnow limitrate;

                payload = struct();
                payload.message = sprintf('Correction applied to signal %d.', signalIndex);
                payload.singleMode = true;
                payload.signalIndex = signalIndex;
                payload.nSamples = nR;
                payload.nChannels = nC;
                payload.wavelength = app.Wavelength;
                payload.signalCorrected = app.CorrectedData(signalIndex, :);
                payload.signalBaseline = app.BaselineData(signalIndex, :);
                payload.meanOriginal = mean(app.OriginalData, 1);
                payload.meanCorrected = mean(app.CorrectedData, 1);

                % Preview baseline on current signal (consistent with what user sees)
                [~, previewBl] = app.Corrector.previewBaseline( ...
                    app.OriginalData(signalIndex, :), lambda, p, struct.empty, false, 1e5, maxIter);
                payload.meanBaseline = previewBl;
            else
                % Batch mode: correct all spectra
                app.CorrectedData = zeros(nR, nC);
                app.BaselineData = zeros(nR, nC);

                progressInterval = max(1, floor(nR / 20));
                sendResponse(app, 'progress', struct('title', 'Applying correction...', 'current', 0, 'total', nR));

                for i = 1:nR
                    [corrRow, blRow] = app.Corrector.correct( ...
                        app.OriginalData(i, :), lambda, p, struct.empty, false, 1e5, maxIter);
                    app.CorrectedData(i, :) = corrRow;
                    app.BaselineData(i, :) = blRow;

                    if mod(i, progressInterval) == 0 || i == nR
                        sendResponse(app, 'progress', struct('title', 'Applying correction...', 'current', i, 'total', nR));
                        drawnow limitrate;
                    end
                end

                % Preview baseline on mean spectrum (consistent with preview)
                [~, previewBl] = app.Corrector.previewBaseline( ...
                    app.OriginalData, lambda, p, struct.empty, false, 1e5, maxIter);

                payload = struct();
                payload.message = sprintf('Correction applied to %d spectra.', nR);
                payload.singleMode = false;
                payload.nSamples = nR;
                payload.nChannels = nC;
                payload.wavelength = app.Wavelength;
                payload.meanOriginal = mean(app.OriginalData, 1);
                payload.meanCorrected = mean(app.CorrectedData, 1);
                payload.meanBaseline = previewBl;

                if nR * nC <= 2000000
                    allCorr = cell(1, nR);
                    for i = 1:nR
                        allCorr{i} = app.CorrectedData(i, :);
                    end
                    payload.allCorrected = allCorr;
                end
            end

            sendResponse(app, 'correctionApplied', payload);
        end

        function handlePrepareExport(app)
            if isempty(app.CorrectedData)
                sendResponse(app, 'error', struct('message', 'Apply correction before exporting.'));
                return;
            end
            existingNames = evalin('base', 'who');
            payload = struct();
            payload.suggestedCorrectedName = suggestName(app, existingNames, 'aslsCorrectedData');
            payload.suggestedBaselineName = suggestName(app, existingNames, 'aslsBaselineData');
            sendResponse(app, 'showExportDialog', payload);
        end

        function handleDoExport(app, data)
            if isempty(app.CorrectedData)
                sendResponse(app, 'error', struct('message', 'No corrected data available.'));
                return;
            end

            correctedName = 'aslsCorrectedData';
            baselineName = 'aslsBaselineData';
            if isfield(data, 'correctedName'), correctedName = strtrim(char(string(data.correctedName))); end
            if isfield(data, 'baselineName'), baselineName = strtrim(char(string(data.baselineName))); end

            if ~isvarname(correctedName)
                sendResponse(app, 'error', struct('message', sprintf('Invalid name: %s', correctedName)));
                return;
            end
            if ~isvarname(baselineName)
                sendResponse(app, 'error', struct('message', sprintf('Invalid name: %s', baselineName)));
                return;
            end

            assignin('base', correctedName, app.CorrectedData);
            assignin('base', baselineName, app.BaselineData);

            payload = struct();
            payload.exportedNames = {correctedName, baselineName};
            payload.message = sprintf('Exported "%s" and "%s".', correctedName, baselineName);
            sendResponse(app, 'exportCompleted', payload);
        end

        function handleCreateDemoData(app)
            rng(42);
            nSamples = 30; nChannels = 500;
            wl = linspace(400, 900, nChannels);
            spectra = zeros(nSamples, nChannels);
            for s = 1:nSamples
                bl = polyval([rand*0.002-0.001, rand*0.5-0.25, rand*20-10, rand*200], linspace(-1,1,nChannels));
                bl = bl + 50*exp(-linspace(0,3,nChannels));
                for pk = [480 550 630 720 810]
                    spectra(s,:) = spectra(s,:) + (20+rand*40)*exp(-0.5*((wl-pk)/(8+rand*12)).^2);
                end
                spectra(s,:) = spectra(s,:) + bl + randn(1, nChannels)*3;
            end

            assignin('base', 'demoSpectra', spectra);
            assignin('base', 'demoWavelength', wl);
            loadDataIntoState(app, spectra, wl, 'demoSpectra');

            lambda = 1e6; p = 1e-3;
            [meanS, blS] = app.Corrector.previewBaseline(spectra, lambda, p, struct.empty, false, 1e5, 10);

            payload = struct();
            payload.message = 'Demo data loaded (30 x 500).';
            payload.varName = 'demoSpectra';
            payload.nSamples = 30;
            payload.nChannels = 500;
            payload.wavelength = wl;
            payload.meanSpectrum = meanS;
            payload.previewBaseline = blS;
            payload.lambda = lambda;
            payload.p = p;
            payload.maxIter = 10;

            allSpec = cell(1, nSamples);
            for i = 1:nSamples
                allSpec{i} = spectra(i, :);
            end
            payload.allSpectra = allSpec;

            sendResponse(app, 'dataLoaded', payload);
        end

        function handleUiReady(app)
            % Send any pending data that was queued during construction
            if ~isempty(app.PendingPayload)
                sendResponse(app, 'dataLoaded', app.PendingPayload);
                app.PendingPayload = struct.empty;
            end
        end

        function handleNext(app)
            app.NextClicked = true;
            if app.IsWaiting
                try uiresume(app.UIFigure); catch; end
            end
            sendResponse(app, 'statusUpdate', struct('type', 'success', 'message', 'Continuing.'));
        end
    end

    % =====================================================================
    %  HELPERS
    % =====================================================================
    methods (Access = private)
        function sendResponse(app, responseType, payload)
            if app.IsClosed, return; end
            try
                r = struct();
                r.response = responseType;
                r.payload = payload;
                r.source = 'matlab';
                app.ResponseCounter = app.ResponseCounter + 1;
                r.timestamp = app.ResponseCounter;
                app.HTMLComponent.Data = r;
            catch
            end
        end

        function loadDataIntoState(app, data, wavelength, varName)
            app.OriginalData = double(data);
            app.CorrectedData = [];
            app.BaselineData = [];
            app.DataLoaded = true;
            app.LoadedVarName = char(varName);
            if isempty(wavelength)
                app.Wavelength = 1:size(data, 2);
            else
                app.Wavelength = double(wavelength(:)');
            end
        end

        function variables = collectWorkspaceMatrices(~)
            vars = evalin('base', 'whos');
            variables = {};
            for i = 1:numel(vars)
                if vars(i).global, continue; end
                if ~ismember(vars(i).class, {'double','single','int8','int16','int32','int64','uint8','uint16','uint32','uint64'})
                    continue;
                end
                if numel(vars(i).size) ~= 2, continue; end
                nR = vars(i).size(1); nC = vars(i).size(2);
                if nR < 1 || nC < 4, continue; end
                variables{end+1} = struct('name', vars(i).name, 'rows', nR, 'cols', nC, ...
                    'size', sprintf('%d x %d', nR, nC)); %#ok<AGROW>
            end
        end

        function vectors = collectWorkspaceVectors(~)
            vars = evalin('base', 'whos');
            vectors = {};
            for i = 1:numel(vars)
                if vars(i).global, continue; end
                if ~ismember(vars(i).class, {'double','single','int8','int16','int32','int64','uint8','uint16','uint32','uint64'})
                    continue;
                end
                if numel(vars(i).size) ~= 2, continue; end
                nR = vars(i).size(1); nC = vars(i).size(2);
                if ~(nR == 1 || nC == 1), continue; end
                len = max(nR, nC);
                if len < 2, continue; end
                vectors{end+1} = struct('name', vars(i).name, 'length', len); %#ok<AGROW>
            end
        end

        function val = safeDouble(~, s, fieldName, defaultVal)
            val = defaultVal;
            if ~isstruct(s) || ~isfield(s, fieldName), return; end
            raw = s.(fieldName);
            if ischar(raw) || isstring(raw)
                parsed = str2double(raw);
                if isfinite(parsed), val = parsed; end
            elseif isnumeric(raw) && isscalar(raw) && isfinite(raw)
                val = double(raw);
            elseif islogical(raw) && isscalar(raw)
                val = double(raw);
            end
        end

        function [wl, sOrig, sCorr, mOrig, mCorr] = sampleForPlot(~, original, corrected, wavelength, maxRows, maxCols)
            nRows = size(original, 1);
            nCols = size(original, 2);
            rowStep = max(1, ceil(nRows / maxRows));
            colStep = max(1, ceil(nCols / maxCols));
            rowIdx = 1:rowStep:nRows;
            colIdx = 1:colStep:nCols;
            wl = wavelength(colIdx);
            sOrig = original(rowIdx, colIdx);
            sCorr = corrected(rowIdx, colIdx);
            mOrig = mean(original(:, colIdx), 1);
            mCorr = mean(corrected(:, colIdx), 1);
        end

        function outName = suggestName(~, existingNames, baseName)
            outName = baseName;
            if ~iscell(existingNames), return; end
            if ~any(strcmp(existingNames, outName)), return; end
            suffix = 1;
            while any(strcmp(existingNames, sprintf('%s_%d', baseName, suffix)))
                suffix = suffix + 1;
            end
            outName = sprintf('%s_%d', baseName, suffix);
        end
    end

    % =====================================================================
    %  PUBLIC API
    % =====================================================================
    methods (Access = public)
        function out = getData(app)
            out = struct();
            out.data = app.CorrectedData;
            if isempty(out.data), out.data = app.OriginalData; end
            out.correctedData = app.CorrectedData;
            out.baselineData = app.BaselineData;
            out.wavelength = app.Wavelength;
            out.varName = app.LoadedVarName;
        end

        function setInputData(app, inputData)
            if isnumeric(inputData)
                data = inputData;
                wavelength = 1:size(data, 2);
                varName = 'inputData';
            elseif isstruct(inputData) && isfield(inputData, 'data')
                data = inputData.data;
                wavelength = 1:size(data, 2);
                if isfield(inputData, 'wavelength') && ~isempty(inputData.wavelength)
                    wavelength = inputData.wavelength;
                end
                varName = 'inputData';
                if isfield(inputData, 'varName'), varName = char(string(inputData.varName)); end
            else
                error('Input must be a numeric matrix or struct with "data" field.');
            end

            [ok, msg] = app.Validator.validateData(data);
            if ~ok, error(msg); end

            loadDataIntoState(app, double(data), wavelength, varName);
            [meanS, blS] = app.Corrector.previewBaseline(app.OriginalData, 1e6, 1e-3, struct.empty, false, 1e5, 10);

            payload = struct();
            payload.message = sprintf('Data loaded (%d x %d).', size(data, 1), size(data, 2));
            payload.varName = varName;
            payload.nSamples = size(data, 1);
            payload.nChannels = size(data, 2);
            payload.wavelength = app.Wavelength;
            payload.meanSpectrum = meanS;
            payload.previewBaseline = blS;
            payload.lambda = 1e6;
            payload.p = 1e-3;
            payload.maxIter = 10;

            % Store as pending (sent when JS fires 'uiReady') or send now
            app.PendingPayload = payload;
            if ~isempty(app.HTMLComponent) && isvalid(app.HTMLComponent)
                try
                    sendResponse(app, 'dataLoaded', payload);
                catch
                end
            end
        end

        function waitForOrchestrator(app)
        % waitForOrchestrator  Block until user clicks Next or window closes.
            app.IsWaiting = true;
            app.NextClicked = false;
            while ~app.IsClosed && ~app.NextClicked
                uiwait(app.UIFigure);
            end
            app.IsWaiting = false;
        end
    end
end
