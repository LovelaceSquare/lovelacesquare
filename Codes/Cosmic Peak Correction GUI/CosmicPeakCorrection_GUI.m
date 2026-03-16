classdef CosmicPeakCorrection_GUI < matlab.apps.AppBase
% CosmicPeakCorrection_GUI. Interactive GUI for cosmic ray spike removal.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Detects and removes cosmic ray spikes (sharp, narrow artifacts caused by
% cosmic radiation hitting the CCD detector) from spectral datasets. The
% algorithm computes the absolute derivative of each spectrum, flags
% channels where the derivative exceeds a threshold, expands the flagged
% region by k channels on each side, and replaces the affected values with
% linearly interpolated data from neighbouring clean channels.
%
% The GUI provides two interactive charts:
%   - Spectra chart: all individual spectra overlaid (original + corrected),
%     with signal-by-signal browsing and show-all/corrected-only toggles.
%   - Derivative chart: pointwise max |derivative| across all samples, with
%     sorted rank mode, click-to-set threshold, and per-sample background.
%
% Uses the AppBase + uihtml architecture: HTML/CSS/JS frontend with a
% MATLAB backend. Charts rendered with Canvas 2D (no MATLAB axes).
%
% ALGORITHM (per sample):
%   1. Compute abs(diff(spectrum, derivativeOrder)) along channels.
%   2. Flag channels where |derivative| > threshold.
%   3. Expand each flagged position by channelsToRemove on each side.
%   4. Set flagged channels to NaN.
%   5. Fill NaN gaps with fillmissing('linear', 'EndValues', 'nearest').
%
% PARAMETERS:
%   derivativeOrder  - Order of the finite difference (default: 1)
%   channelsToRemove - Expansion radius around each spike (default: 2)
%   threshold        - Detection threshold on |derivative| (default: auto)
%
% GUI FEATURES:
%   - Load workspace variables or generate built-in demo data
%   - Logarithmic threshold slider with manual numeric input
%   - Click on derivative chart to set threshold visually
%   - Sort derivative by rank (descending) or by channel position
%   - Signal-by-signal navigation with bold highlight
%   - Show all / corrected only toggle badges
%   - Zoom (mouse wheel) and pan (click-drag) on both charts
%   - Progress bar for loading, correction, and export
%   - Large dataset guard with subsampling options (> 200 spectra)
%   - Dark mode toggle (Ctrl+D)
%   - Resizable control/chart panels via drag handle
%   - Export corrected data to MATLAB workspace
%
% USAGE:
%   app = CosmicPeakCorrection_GUI();              % Empty (load from GUI)
%   app = CosmicPeakCorrection_GUI(spectra);       % Direct numeric matrix
%   app = CosmicPeakCorrection_GUI(s);             % Struct with .data, .wavelength
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
%                               .data and .wavelength fields)
%   app.waitForOrchestrator() - Block until user clicks Next or closes
%
% EXAMPLE:
%   % Generate demo data and launch GUI
%   CosmicPeak_test;
%   s.data = spectra; s.wavelength = wavelength;
%   app = CosmicPeakCorrection_GUI(s);
%
%   % Use in a pipeline (orchestrator mode)
%   app = CosmicPeakCorrection_GUI(s);
%   app.waitForOrchestrator();
%   result = app.getData();
%   correctedSpectra = result.data;
%
% See also: CosmicPeakCorrector, DataValidator
%
% Disclaimer:
%   Author and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this software.

    % ------------------------------------------------------------------ %
    %  Properties
    % ------------------------------------------------------------------ %
    properties (Access = public)
        Figure              matlab.ui.Figure
        HTMLComponent       matlab.ui.control.HTML
    end

    properties (Access = private)
        Corrector           CosmicPeakCorrector
        Validator           DataValidator
        OriginalData        double = []
        CorrectedData       double = []
        Wavelength          double = []
        LoadedVarName       char = ''
        DataLoaded          logical = false
        IsClosed            logical = false
        IsWaiting           logical = false
        NextClicked         logical = false
        ResponseCounter     double = 0
        PendingPayload      struct = struct.empty
    end

    % ------------------------------------------------------------------ %
    %  Constructor / Destructor
    % ------------------------------------------------------------------ %
    methods (Access = public)

        function app = CosmicPeakCorrection_GUI(inputData)
            createComponents(app);
            initializeBusinessLogic(app);
            registerApp(app, app.Figure);

            if nargin >= 1
                try
                    setInputData(app, inputData);
                catch ME
                    warning('CosmicPeakCorrection_GUI:inputError', ...
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
            if ~isempty(app.Figure) && isvalid(app.Figure)
                delete(app.Figure);
            end
        end
    end

    % ------------------------------------------------------------------ %
    %  Component initialisation
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function createComponents(app)
            screenSize = get(0, 'ScreenSize');
            figWidth = min(1400, round(screenSize(3) * 0.8));
            figHeight = min(900, round(screenSize(4) * 0.85));

            app.Figure = uifigure('Visible', 'off');
            app.Figure.Position = [100 100 figWidth figHeight];
            app.Figure.Name = 'Cosmic Peak Correction';
            app.Figure.Color = [0.91 0.92 0.93];
            app.Figure.AutoResizeChildren = 'off';
            app.Figure.CloseRequestFcn = createCallbackFcn(app, @closeRequest, true);
            app.Figure.SizeChangedFcn = createCallbackFcn(app, @figureSizeChanged, true);

            modulePath = fileparts(mfilename('fullpath'));
            htmlPath = fullfile(modulePath, 'ui', 'cosmic_peak_correction_ui.html');

            app.HTMLComponent = uihtml(app.Figure);
            app.HTMLComponent.Position = [1 1 figWidth figHeight];
            app.HTMLComponent.HTMLSource = htmlPath;
            app.HTMLComponent.DataChangedFcn = createCallbackFcn(app, @HTMLDataChanged, true);

            app.Figure.Visible = 'on';
        end

        function initializeBusinessLogic(app)
            modulePath = fileparts(mfilename('fullpath'));
            addpath(fullfile(modulePath, 'business_logic'));
            app.Corrector = CosmicPeakCorrector();
            app.Validator = DataValidator();
        end

        function startupFcn(app)
            movegui(app.Figure, 'center');
            sendResponse(app, 'ready', struct('message', 'Ready'));
        end

        function figureSizeChanged(app, ~)
            if isempty(app.Figure) || ~isvalid(app.Figure), return; end
            pos = app.Figure.Position;
            if ~isempty(app.HTMLComponent) && isvalid(app.HTMLComponent)
                app.HTMLComponent.Position = [1 1 pos(3) pos(4)];
            end
        end

        function closeRequest(app, ~)
            app.IsClosed = true;
            if app.IsWaiting
                try uiresume(app.Figure); catch; end
            end
            delete(app);
        end
    end

    % ------------------------------------------------------------------ %
    %  Action Dispatcher
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function HTMLDataChanged(app, ~)
            try
                data = app.HTMLComponent.Data;
                if ~isstruct(data) || ~isfield(data, 'action'), return; end
                if isfield(data, 'source') && strcmp(data.source, 'matlab'), return; end

                action = char(string(data.action));
                switch action
                    case 'uiReady',         handleUiReady(app);
                    case 'loadData',        handleLoadData(app);
                    case 'loadVariable',    handleLoadVariable(app, data);
                    case 'apply',           handleApply(app, data);
                    case 'prepareExport',   handlePrepareExport(app);
                    case 'doExport',        handleDoExport(app, data);
                    case 'createDemoData',  handleCreateDemoData(app);
                    case 'getSpectrum',     handleGetSpectrum(app, data);
                    case 'next',            handleNext(app);
                end
            catch ME
                sendResponse(app, 'error', struct('message', ME.message));
            end
        end
    end

    % ------------------------------------------------------------------ %
    %  Action Handlers
    % ------------------------------------------------------------------ %
    methods (Access = private)

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

            sendResponse(app, 'progress', struct('title', 'Computing auto threshold...', 'percent', 40));
            drawnow limitrate;

            autoThresh = app.Corrector.autoThreshold(app.OriginalData, 1);

            nR = size(raw, 1);
            payload = struct();
            payload.message = sprintf('Loaded "%s" (%d x %d).', varName, nR, nChannels);
            payload.varName = varName;
            payload.nSamples = nR;
            payload.nChannels = nChannels;
            payload.wavelength = app.Wavelength;
            payload.meanSpectrum = mean(app.OriginalData, 1);
            payload.autoThreshold = autoThresh;

            if nR * nChannels <= 2000000
                sendResponse(app, 'progress', struct('title', 'Preparing spectra data...', 'percent', 60));
                drawnow limitrate;
                allSpec = cell(1, nR);
                specProgress = max(1, floor(nR / 10));
                for i = 1:nR
                    allSpec{i} = app.OriginalData(i, :);
                    if mod(i, specProgress) == 0
                        pct = 60 + round(35 * i / nR);
                        sendResponse(app, 'progress', struct('title', 'Preparing spectra data...', 'percent', pct, 'current', i, 'total', nR));
                        drawnow limitrate;
                    end
                end
                payload.allSpectra = allSpec;
            end

            sendResponse(app, 'dataLoaded', payload);
        end

        function handleApply(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'error', struct('message', 'Load data first.'));
                return;
            end

            derivativeOrder = max(1, round(safeDouble(app, data, 'derivativeOrder', 1)));
            channelsToRemove = max(0, round(safeDouble(app, data, 'channelsToRemove', 2)));
            threshold = safeDouble(app, data, 'threshold', 5);

            nR = size(app.OriginalData, 1);
            nC = size(app.OriginalData, 2);

            [valid, msg] = app.Validator.validateParams(derivativeOrder, channelsToRemove, threshold, nC);
            if ~valid
                sendResponse(app, 'error', struct('message', msg));
                return;
            end

            sendResponse(app, 'progress', struct('title', 'Applying correction...', 'percent', 10));
            drawnow limitrate;

            [correctedData, correctionMask] = app.Corrector.correct( ...
                app.OriginalData, derivativeOrder, channelsToRemove, threshold);

            app.CorrectedData = correctedData;

            sendResponse(app, 'progress', struct('title', 'Preparing results...', 'percent', 70));
            drawnow limitrate;

            correctedSamples = find(any(correctionMask, 2));
            nCorrected = length(correctedSamples);

            payload = struct();
            payload.message = sprintf('Correction applied. %d / %d samples corrected.', nCorrected, nR);
            payload.nSamples = nR;
            payload.nChannels = nC;
            payload.nCorrectedSamples = nCorrected;
            payload.wavelength = app.Wavelength;
            payload.meanOriginal = mean(app.OriginalData, 1);
            payload.meanCorrected = mean(app.CorrectedData, 1);

            if nR * nC <= 2000000
                sendResponse(app, 'progress', struct('title', 'Preparing corrected spectra...', 'percent', 80));
                drawnow limitrate;
                allCorr = cell(1, nR);
                for i = 1:nR
                    allCorr{i} = app.CorrectedData(i, :);
                end
                payload.allCorrected = allCorr;
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
            payload.suggestedName = suggestName(app, existingNames, 'cosmicCorrectedData');
            sendResponse(app, 'showExportDialog', payload);
        end

        function handleDoExport(app, data)
            if isempty(app.CorrectedData)
                sendResponse(app, 'error', struct('message', 'No corrected data available.'));
                return;
            end

            exportName = 'cosmicCorrectedData';
            if isfield(data, 'exportName')
                exportName = strtrim(char(string(data.exportName)));
            end

            if ~isvarname(exportName)
                sendResponse(app, 'error', struct('message', sprintf('Invalid variable name: %s', exportName)));
                return;
            end

            sendResponse(app, 'progress', struct('title', 'Exporting...', 'percent', 50));
            drawnow limitrate;

            assignin('base', exportName, app.CorrectedData);

            payload = struct();
            payload.exportedName = exportName;
            payload.message = sprintf('Exported "%s" (%d x %d) to workspace.', ...
                exportName, size(app.CorrectedData, 1), size(app.CorrectedData, 2));
            sendResponse(app, 'exportCompleted', payload);
        end

        function handleCreateDemoData(app)
            sendResponse(app, 'progress', struct('title', 'Generating demo data...', 'percent', 20));
            drawnow limitrate;

            rng(42);
            nS = 40; nC = 800;
            wl = linspace(200, 3500, nC);
            gauss = @(x, mu, h, w) h .* exp(-((x - mu).^2) ./ (2*w.^2));
            spectra = zeros(nS, nC);

            for s = 1:nS
                bg = 500 + 200*randn() + gauss(wl, 1800, 400+50*randn(), 800);
                peaks = zeros(1, nC);
                peakPos = [520, 780, 1000, 1200, 1450, 1600, 2200, 2900, 3100];
                peakH   = [80, 60, 120, 90, 150, 200, 50, 180, 100];
                peakW   = [15, 12, 20, 18, 25, 22, 30, 40, 35];
                for p = 1:length(peakPos)
                    h = peakH(p) * (0.7 + 0.6*rand());
                    peaks = peaks + gauss(wl, peakPos(p)+3*(rand()-0.5), h, peakW(p));
                end
                spectra(s, :) = bg + peaks + 5*randn(1, nC);
            end

            % Insert cosmic ray spikes
            nAffected = round(nS * 0.3);
            affectedIdx = randperm(nS, nAffected);
            for i = 1:nAffected
                s = affectedIdx(i);
                nSpikes = randi([2 5]);
                spikePos = randi([5, nC-4], 1, nSpikes);
                for j = 1:nSpikes
                    spikeHeight = 2000 + 3000*rand();
                    col = spikePos(j);
                    spectra(s, col) = spectra(s, col) + spikeHeight;
                    if rand() > 0.5 && col < nC
                        spectra(s, col+1) = spectra(s, col+1) + spikeHeight*0.3;
                    end
                end
            end

            sendResponse(app, 'progress', struct('title', 'Loading demo data...', 'percent', 60));
            drawnow limitrate;

            assignin('base', 'demoSpectra', spectra);
            assignin('base', 'demoWavelength', wl);
            loadDataIntoState(app, spectra, wl, 'demoSpectra');

            autoThresh = app.Corrector.autoThreshold(spectra, 1);

            payload = struct();
            payload.message = sprintf('Demo data loaded (%d x %d, %d with spikes).', nS, nC, nAffected);
            payload.varName = 'demoSpectra';
            payload.nSamples = nS;
            payload.nChannels = nC;
            payload.wavelength = wl;
            payload.meanSpectrum = mean(spectra, 1);
            payload.autoThreshold = autoThresh;

            allSpec = cell(1, nS);
            for i = 1:nS
                allSpec{i} = spectra(i, :);
            end
            payload.allSpectra = allSpec;

            sendResponse(app, 'dataLoaded', payload);
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
            sendResponse(app, 'spectrumLoaded', payload);
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
                try uiresume(app.Figure); catch; end
            end
            sendResponse(app, 'statusUpdate', struct('type', 'success', 'message', 'Continuing.'));
        end
    end

    % ------------------------------------------------------------------ %
    %  Helpers
    % ------------------------------------------------------------------ %
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
                if nR < 1 || nC < 3, continue; end
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

    % ------------------------------------------------------------------ %
    %  Public API (orchestrator interface)
    % ------------------------------------------------------------------ %
    methods (Access = public)

        function result = getData(app)
        % getData  Return the corrected data (or original if not yet corrected).
            result = struct();
            if ~isempty(app.CorrectedData)
                result.data = app.CorrectedData;
            else
                result.data = app.OriginalData;
            end
            result.wavelength = app.Wavelength;
        end

        function setInputData(app, dataStruct)
        % setInputData  Accept data from the orchestrator.
            if isnumeric(dataStruct)
                rawData = dataStruct;
                wl = 1:size(rawData, 2);
                vName = 'inputData';
            elseif isstruct(dataStruct) && isfield(dataStruct, 'data')
                rawData = dataStruct.data;
                wl = 1:size(rawData, 2);
                if isfield(dataStruct, 'wavelength') && ~isempty(dataStruct.wavelength)
                    wl = dataStruct.wavelength;
                end
                vName = 'inputData';
                if isfield(dataStruct, 'varName')
                    vName = char(string(dataStruct.varName));
                end
            else
                error('Input must be a numeric matrix or struct with "data" field.');
            end

            [valid, msg] = app.Validator.validateData(rawData);
            if ~valid, error(msg); end

            loadDataIntoState(app, double(rawData), wl, vName);

            autoThresh = app.Corrector.autoThreshold(app.OriginalData, 1);
            nR = size(app.OriginalData, 1);
            nC = size(app.OriginalData, 2);

            payload = struct();
            payload.message = sprintf('Data loaded (%d x %d).', nR, nC);
            payload.varName = vName;
            payload.nSamples = nR;
            payload.nChannels = nC;
            payload.wavelength = app.Wavelength;
            payload.meanSpectrum = mean(app.OriginalData, 1);
            payload.autoThreshold = autoThresh;

            if nR * nC <= 2000000
                allSpec = cell(1, nR);
                for i = 1:nR
                    allSpec{i} = app.OriginalData(i, :);
                end
                payload.allSpectra = allSpec;
            end

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
                uiwait(app.Figure);
            end
            app.IsWaiting = false;
        end
    end
end
