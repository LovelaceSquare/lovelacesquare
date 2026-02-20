classdef LASLS < matlab.apps.AppBase
    % LASLS - HTML-based GUI for Local Asymmetric Least Squares (LAsLS)
    % baseline correction.
    %
    % Interactive graphical interface for LAsLS baseline correction with
    % per-interval asymmetry and smoothing parameters, signal-by-signal
    % navigation, per-signal parameter mode, peak detection, interactive
    % interval drawing, and first-derivative penalty (mu).
    %
    % Architecture: MATLAB AppBase + uihtml (HTML/JS frontend, MATLAB backend)
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

    properties (Access = public)
        UIFigure        matlab.ui.Figure
        HTMLComponent   matlab.ui.control.HTML
    end

    properties (Access = private)
        Corrector                           % LASLSCorrector instance
        Validator                           % DataValidator instance
        OriginalData    double = []         % Raw spectra matrix (nRows x nCols)
        CorrectedData   double = []         % Corrected spectra
        BaselineData    double = []         % Estimated baselines
        WeightsData     double = []         % IRLS weights
        Wavelength      double = []         % Channel axis (1:nCols)
        DataLoaded      logical = false
        IsClosed        logical = false
        UIUpdateCounter double = 0
        LoadedVarName   char   = ''
    end

    % ====================================================================
    %  CONSTRUCTOR / DESTRUCTOR
    % ====================================================================
    methods (Access = public)

        function app = LASLS()
            createComponents(app);
            initializeBusinessLogic(app);
            registerApp(app, app.UIFigure);
            runStartupFcn(app, @startupFcn);
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            app.IsClosed = true;
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end

    end

    % ====================================================================
    %  COMPONENT CREATION
    % ====================================================================
    methods (Access = private)

        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1400 900];
            app.UIFigure.Name = 'LASLS Baseline Correction';
            app.UIFigure.Color = [0.91 0.92 0.93];
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.SizeChangedFcn  = createCallbackFcn(app, @UIFigureSizeChanged, true);

            modulePath = fileparts(mfilename('fullpath'));
            htmlPath = fullfile(modulePath, 'ui', 'lasls_baseline_correction_ui.html');

            app.HTMLComponent = uihtml(app.UIFigure);
            app.HTMLComponent.Position = [1 1 1400 900];
            app.HTMLComponent.HTMLSource = htmlPath;
            app.HTMLComponent.DataChangedFcn = createCallbackFcn(app, @HTMLDataChanged, true);

            app.UIFigure.Visible = 'on';
            app.UIFigure.WindowState = 'maximized';
        end

        function initializeBusinessLogic(app)
            modulePath = fileparts(mfilename('fullpath'));
            blPath = fullfile(modulePath, 'business_logic');
            addpath(blPath);
            app.Corrector = LASLSCorrector();
            app.Validator = DataValidator();
        end

        function startupFcn(app)
            movegui(app.UIFigure, 'center');
            pause(0.3);
            sendResponse(app, 'statusUpdate', struct('type', 'idle', 'message', 'Ready'));
        end

    end

    % ====================================================================
    %  UI CALLBACKS
    % ====================================================================
    methods (Access = private)

        function UIFigureCloseRequest(app, ~)
            app.IsClosed = true;
            delete(app);
        end

        function UIFigureSizeChanged(app, ~)
            pos = app.UIFigure.Position;
            app.HTMLComponent.Position = [1 1 pos(3) pos(4)];
        end

        % ----------------------------------------------------------------
        %  HTML DataChanged dispatcher
        %  IMPORTANT: Each handler must send exactly ONE sendResponse call.
        %  Two rapid .Data writes cause the first to be lost (race condition).
        % ----------------------------------------------------------------
        function HTMLDataChanged(app, ~)
            try
                data = app.HTMLComponent.Data;

                if isempty(data) || ~isstruct(data), return; end
                if isfield(data, 'response'), return; end
                if ~isfield(data, 'action') || isempty(data.action), return; end

                action = data.action;
                if ~ischar(action)
                    try, action = char(action); catch, return; end
                end

                switch action
                    case 'loadData',            handleLoadData(app);
                    case 'loadVariable',        handleLoadVariable(app, data);
                    case 'getSpectrum',          handleGetSpectrum(app, data);
                    case 'previewBaseline',      handlePreviewBaseline(app, data);
                    case 'apply',                handleApply(app, data);
                    case 'prepareExport',        handlePrepareExport(app);
                    case 'doExport',             handleDoExport(app, data);
                    case 'prepareImport',        handlePrepareImport(app);
                    case 'doImport',             handleDoImport(app, data);
                    case 'exportIntervals',      handleExportIntervals(app, data);
                    case 'loadIntervals',        handleLoadIntervals(app);
                    case 'loadIntervalTable',    handleLoadIntervalTable(app, data);
                    case 'createDemoData',       handleCreateDemoData(app, data);
                    otherwise
                        sendResponse(app, 'statusUpdate', ...
                            struct('type', 'error', 'message', ['Unknown action: ' action]));
                end
            catch ME
                disp('=== LASLS Action Error ===');
                disp(ME.getReport());
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Error: ' ME.message]));
            end
        end

    end

    % ====================================================================
    %  ACTION HANDLERS
    %  RULE: Each handler sends exactly ONE sendResponse with ALL data
    %  embedded (including status). Never two consecutive sendResponse calls.
    % ====================================================================
    methods (Access = private)

        % ---- Load data from workspace (send variable list to JS) -----------
        function handleLoadData(app)
            try
                vars = evalin('base', 'whos');
                varList = {};
                for i = 1:length(vars)
                    if (strcmp(vars(i).class, 'double') || ...
                        strcmp(vars(i).class, 'single')) && ...
                        length(vars(i).size) == 2 && all(vars(i).size > 0)
                        % Include variable info for display
                        varList{end+1} = struct(...
                            'name', vars(i).name, ...
                            'size', sprintf('%dx%d', vars(i).size(1), vars(i).size(2)), ...
                            'rows', vars(i).size(1), ...
                            'cols', vars(i).size(2)); %#ok<AGROW>
                    end
                end

                % Always show modal, even if empty (so user knows it's working)
                sendResponse(app, 'showVarList', struct('variables', {varList}));

            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Error: ' ME.message]));
            end
        end

        % ---- Load a specific variable by name -----------------------------
        function handleLoadVariable(app, data)
            try
                selectedName = data.varName;

                if ~isvarname(selectedName)
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', ['Invalid variable name: "' selectedName '"']));
                    return;
                end

                if ~evalin('base', ['exist(''' selectedName ''', ''var'')'])
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', ['Variable "' selectedName '" not found']));
                    return;
                end

                rawData = evalin('base', selectedName);

                if ~isreal(rawData), rawData = real(rawData); end
                [isValid, msg] = app.Validator.validateData(rawData);
                if ~isValid
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', msg));
                    return;
                end

                app.OriginalData = double(rawData);
                [nRows, nCols] = size(rawData);
                app.Wavelength = (1:nCols)';
                app.DataLoaded = true;
                app.LoadedVarName = selectedName;
                app.CorrectedData = [];
                app.BaselineData = [];
                app.WeightsData = [];

                % Build payload with all needed data
                payload = struct();
                payload.wavelength = app.Wavelength';
                payload.meanSpectrum = mean(app.OriginalData, 1);
                payload.nSamples = nRows;
                payload.nChannels = nCols;
                payload.varName = selectedName;
                payload.statusType = 'success';
                payload.statusMessage = sprintf('Loaded "%s": %dx%d', selectedName, nRows, nCols);

                % Send first spectrum for signal mode
                if nRows >= 1
                    payload.firstSpectrum = app.OriginalData(1, :);
                end

                % If dataset small enough, send all spectra for "show all data"
                if nRows <= 500 && nRows * nCols <= 500000
                    allData = cell(1, nRows);
                    for i = 1:nRows
                        allData{i} = app.OriginalData(i, :);
                    end
                    payload.allSpectra = allData;
                end

                sendResponse(app, 'dataLoaded', payload);

            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Load error: ' ME.message]));
            end
        end

        % ---- Get specific spectrum ----------------------------------------
        function handleGetSpectrum(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'No data loaded'));
                return;
            end

            idx = round(safeDouble(app, data, 'signalIndex', 1));
            nRows = size(app.OriginalData, 1);
            idx = max(1, min(nRows, idx));

            payload = struct();
            payload.signalIndex = idx;
            payload.spectrum = app.OriginalData(idx, :);
            payload.wavelength = app.Wavelength';
            payload.nSamples = nRows;
            payload.statusType = 'success';
            payload.statusMessage = sprintf('Signal %d/%d', idx, nRows);

            sendResponse(app, 'spectrumLoaded', payload);
        end

        % ---- Preview baseline on displayed spectrum -----------------------
        function handlePreviewBaseline(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'No data loaded'));
                return;
            end

            lambda    = safeDouble(app, data, 'lambda', 10);
            p         = safeDouble(app, data, 'p', 0.5);
            mu        = safeDouble(app, data, 'mu', 10);
            maxIter   = round(safeDouble(app, data, 'maxIter', 50));
            tolerance = safeDouble(app, data, 'tolerance', 1e-6);

            % Get the spectrum to use for preview
            signalMode = false;
            if isfield(data, 'signalMode')
                signalMode = logical(data.signalMode);
            end
            signalIndex = round(safeDouble(app, data, 'signalIndex', 1));

            if signalMode && signalIndex >= 1 && signalIndex <= size(app.OriginalData, 1)
                y = app.OriginalData(signalIndex, :)';
            else
                y = mean(app.OriginalData, 1)';
            end

            % Build interval arrays from JS data (now using flat arrays)
            intervals = [];
            pVals = [];
            lambdasAsym = [];

            if isfield(data, 'intStarts') && isfield(data, 'intEnds')
                starts = toDoubleArray(app, data.intStarts);
                ends = toDoubleArray(app, data.intEnds);
                lambdas = toDoubleArray(app, data.intLambdas);
                ps = toDoubleArray(app, data.intPs);

                nInt = numel(starts);
                if nInt > 0
                    intervals = [starts(:), ends(:)];
                    lambdasAsym = lambdas(:);
                    pVals = ps(:);
                end
            else
                % Fallback to old format
                [intervals, pVals, lambdasAsym] = buildIntervalsFromData(app, data);
            end

            try
                [baseline, ~] = app.Corrector.computeBaseline(y, intervals, pVals, ...
                    lambdasAsym, lambda, mu, maxIter, tolerance, p);

                result = struct();
                result.baseline = baseline';
                result.corrected = y' - baseline';
                result.wavelength = app.Wavelength';
                result.statusType = 'success';
                result.statusMessage = 'Baseline updated';

                sendResponse(app, 'previewUpdated', result);
            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Preview: ' ME.message]));
            end
        end

        % ---- Apply correction to all spectra ------------------------------
        function handleApply(app, data)
            if ~app.DataLoaded
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'No data loaded'));
                return;
            end

            lambda    = safeDouble(app, data, 'lambda', 10);
            p         = safeDouble(app, data, 'p', 0.5);
            mu        = safeDouble(app, data, 'mu', 10);
            maxIter   = round(safeDouble(app, data, 'maxIter', 50));
            tolerance = safeDouble(app, data, 'tolerance', 1e-6);

            perSignalMode = false;
            if isfield(data, 'perSignalMode')
                perSignalMode = logical(data.perSignalMode);
            end

            nRows = size(app.OriginalData, 1);
            nCols = size(app.OriginalData, 2);

            app.CorrectedData = zeros(nRows, nCols);
            app.BaselineData  = zeros(nRows, nCols);
            app.WeightsData   = zeros(nRows, nCols);

            try
                % Send progress updates every N samples (or at least 20 updates total)
                progressInterval = max(1, floor(nRows / 20));

                % Send initial progress (0 of total) so user sees the total
                sendResponse(app, 'correctionProgress', struct(...
                    'current', 0, 'total', nRows));
                drawnow limitrate;

                if perSignalMode && isfield(data, 'psIntervalSignalIdx')
                    % Per-signal mode: flat arrays encode per-signal intervals
                    allSigIdx = toDoubleArray(app, data.psIntervalSignalIdx);
                    allStarts = toDoubleArray(app, data.psIntervalStart);
                    allEnds   = toDoubleArray(app, data.psIntervalEnd);
                    allLam    = toDoubleArray(app, data.psIntervalLambda);
                    allPv     = toDoubleArray(app, data.psIntervalP);

                    % Per-signal global params
                    psLambda = toDoubleArray(app, data.psGlobalLambda);
                    psP      = toDoubleArray(app, data.psGlobalP);
                    psMu     = toDoubleArray(app, data.psGlobalMu);

                    for i = 1:nRows
                        y = app.OriginalData(i, :)';

                        % Get intervals for this signal
                        mask = (allSigIdx == i);
                        sigStarts = allStarts(mask);
                        sigEnds   = allEnds(mask);
                        sigLam    = allLam(mask);
                        sigPv     = allPv(mask);

                        nInt = numel(sigStarts);
                        if nInt > 0
                            sigIntervals = [sigStarts(:), sigEnds(:)];
                        else
                            sigIntervals = [];
                        end

                        % Get global params for this signal
                        if i <= numel(psLambda)
                            sigLambda = psLambda(i);
                            sigGlobalP = psP(i);
                            sigMu = psMu(i);
                        else
                            sigLambda = lambda;
                            sigGlobalP = p;
                            sigMu = mu;
                        end

                        [bl, w] = app.Corrector.computeBaseline(y, sigIntervals, ...
                            sigPv(:), sigLam(:), sigLambda, sigMu, maxIter, tolerance, sigGlobalP);

                        app.BaselineData(i, :)  = bl';
                        app.WeightsData(i, :)   = w';
                        app.CorrectedData(i, :) = app.OriginalData(i, :) - bl';

                        % Send progress update
                        if mod(i, progressInterval) == 0 || i == nRows
                            sendResponse(app, 'correctionProgress', struct(...
                                'current', i, 'total', nRows));
                            drawnow limitrate;
                        end
                    end
                else
                    % Global mode: same intervals/params for all spectra
                    [intervals, pVals, lambdasAsym] = buildIntervalsFromData(app, data);

                    for i = 1:nRows
                        y = app.OriginalData(i, :)';
                        [bl, w] = app.Corrector.computeBaseline(y, intervals, pVals, ...
                            lambdasAsym, lambda, mu, maxIter, tolerance, p);
                        app.BaselineData(i, :)  = bl';
                        app.WeightsData(i, :)   = w';
                        app.CorrectedData(i, :) = app.OriginalData(i, :) - bl';

                        % Send progress update
                        if mod(i, progressInterval) == 0 || i == nRows
                            sendResponse(app, 'correctionProgress', struct(...
                                'current', i, 'total', nRows));
                            drawnow limitrate;
                        end
                    end
                end

                result = struct();
                result.correctedMean = mean(app.CorrectedData, 1);
                result.nCorrected = nRows;
                result.wavelength = app.Wavelength';
                result.statusType = 'success';
                result.statusMessage = sprintf('Corrected %d spectra', nRows);

                % Send all corrected data if small enough for "show all"
                if nRows <= 500 && nRows * nCols <= 500000
                    allCorr = cell(1, nRows);
                    for i = 1:nRows
                        allCorr{i} = app.CorrectedData(i, :);
                    end
                    result.allCorrected = allCorr;
                end

                sendResponse(app, 'correctionApplied', result);
            catch ME
                % Show error in MATLAB command window for debugging
                disp('=== LASLS Apply Error ===');
                disp(ME.message);
                disp(ME.getReport());
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Apply error: ' ME.message]));
            end
        end

        % ---- Prepare export (check existing vars) --------------------------
        function handlePrepareExport(app)
            % Check what data is available
            hasCorrected = ~isempty(app.CorrectedData);
            hasBaseline = ~isempty(app.BaselineData);
            hasWeights = ~isempty(app.WeightsData);

            % Check which default variable names already exist in workspace
            existingVars = {};
            if evalin('base', 'exist(''correctedData'', ''var'')')
                existingVars{end+1} = 'correctedData';
            end
            if evalin('base', 'exist(''baselineData'', ''var'')')
                existingVars{end+1} = 'baselineData';
            end
            if evalin('base', 'exist(''weightsData'', ''var'')')
                existingVars{end+1} = 'weightsData';
            end
            if evalin('base', 'exist(''laslsParams'', ''var'')')
                existingVars{end+1} = 'laslsParams';
            end

            % Send info to JS to show export dialog
            sendResponse(app, 'showExportDialog', struct(...
                'hasCorrected', hasCorrected, ...
                'hasBaseline', hasBaseline, ...
                'hasWeights', hasWeights, ...
                'existingVars', {existingVars}));
        end

        % ---- Do the actual export with user-specified names ----------------
        function handleDoExport(app, data)
            if isempty(app.CorrectedData)
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'No corrected data. Apply first.'));
                return;
            end

            exported = {};

            % Export corrected data
            correctedName = 'correctedData';
            if isfield(data, 'correctedName') && ~isempty(data.correctedName)
                correctedName = data.correctedName;
            end
            if isvarname(correctedName)
                assignin('base', correctedName, app.CorrectedData);
                exported{end+1} = correctedName;
            end

            % Export baseline if available
            if ~isempty(app.BaselineData)
                baselineName = 'baselineData';
                if isfield(data, 'baselineName') && ~isempty(data.baselineName)
                    baselineName = data.baselineName;
                end
                if isvarname(baselineName)
                    assignin('base', baselineName, app.BaselineData);
                    exported{end+1} = baselineName;
                end
            end

            % Export weights if available
            if ~isempty(app.WeightsData)
                weightsName = 'weightsData';
                if isfield(data, 'weightsName') && ~isempty(data.weightsName)
                    weightsName = data.weightsName;
                end
                if isvarname(weightsName)
                    assignin('base', weightsName, app.WeightsData);
                    exported{end+1} = weightsName;
                end
            end

            % Export parameters struct for reproducibility
            if isfield(data, 'paramsName') && ~isempty(data.paramsName) && isfield(data, 'params')
                paramsName = data.paramsName;
                if isvarname(paramsName)
                    % Convert JS params to MATLAB struct
                    params = struct();
                    jsParams = data.params;
                    if isfield(jsParams, 'lambda'), params.lambda = jsParams.lambda; end
                    if isfield(jsParams, 'p'), params.p = jsParams.p; end
                    if isfield(jsParams, 'mu'), params.mu = jsParams.mu; end
                    if isfield(jsParams, 'maxIter'), params.maxIter = jsParams.maxIter; end
                    if isfield(jsParams, 'tolerance'), params.tolerance = jsParams.tolerance; end
                    if isfield(jsParams, 'perSignalMode'), params.perSignalMode = jsParams.perSignalMode; end

                    % Convert intervals
                    if isfield(jsParams, 'intervals')
                        params.intervals = jsParams.intervals;
                    end

                    % Per-signal data
                    if isfield(jsParams, 'perSignalIntervals')
                        params.perSignalIntervals = jsParams.perSignalIntervals;
                    end
                    if isfield(jsParams, 'perSignalGlobalParams')
                        params.perSignalGlobalParams = jsParams.perSignalGlobalParams;
                    end

                    assignin('base', paramsName, params);
                    exported{end+1} = paramsName;
                end
            end

            if isempty(exported)
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'No valid variable names'));
            else
                sendResponse(app, 'exportCompleted', struct('exportedNames', {exported}));
            end
        end

        % ---- Prepare import (list param structs in workspace) --------------
        function handlePrepareImport(app) %#ok<INUSL>
            try
                vars = evalin('base', 'whos');
                varList = {};
                for i = 1:length(vars)
                    if strcmp(vars(i).class, 'struct')
                        % Check if it looks like a LASLS params struct
                        try
                            val = evalin('base', vars(i).name);
                            if isfield(val, 'lambda') || isfield(val, 'intervals')
                                nInt = 0;
                                perSig = false;
                                if isfield(val, 'intervals')
                                    if isstruct(val.intervals)
                                        nInt = numel(val.intervals);
                                    elseif iscell(val.intervals)
                                        nInt = numel(val.intervals);
                                    end
                                end
                                if isfield(val, 'perSignalMode')
                                    perSig = logical(val.perSignalMode);
                                end
                                varList{end+1} = struct(...
                                    'name', vars(i).name, ...
                                    'nIntervals', nInt, ...
                                    'perSignalMode', perSig); %#ok<AGROW>
                            end
                        catch
                            % Skip if can't read
                        end
                    end
                end

                sendResponse(app, 'showImportDialog', struct('variables', {varList}));
            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Import error: ' ME.message]));
            end
        end

        % ---- Do the actual import ------------------------------------------
        function handleDoImport(app, data) %#ok<INUSL>
            try
                varName = data.varName;
                if ~isvarname(varName)
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', ['Invalid variable name: "' varName '"']));
                    return;
                end
                if ~evalin('base', ['exist(''' varName ''', ''var'')'])
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', ['Variable "' varName '" not found']));
                    return;
                end

                params = evalin('base', varName);
                sendResponse(app, 'paramsImported', params);
            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Import error: ' ME.message]));
            end
        end

        % ---- Export interval table to workspace ---------------------------
        function handleExportIntervals(app, data)
            if ~isfield(data, 'tableStarts') || isempty(data.tableStarts)
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'No intervals to export'));
                return;
            end

            try
                starts = toDoubleArray(app, data.tableStarts);
                ends   = toDoubleArray(app, data.tableEnds);
                lams   = toDoubleArray(app, data.tableLambdas);
                ps     = toDoubleArray(app, data.tablePs);
                nInt   = numel(starts);

                T = table();

                hasPerSignal = isfield(data, 'tableSignalIdx') && ~isempty(data.tableSignalIdx);
                if hasPerSignal
                    sigIdx = toDoubleArray(app, data.tableSignalIdx);
                    T.SignalIndex = sigIdx(:);
                    if isfield(data, 'tableLambdaOut')
                        T.LambdaOut = toDoubleArray(app, data.tableLambdaOut);
                        T.LambdaOut = T.LambdaOut(:);
                    end
                    if isfield(data, 'tablePOut')
                        T.pOut = toDoubleArray(app, data.tablePOut);
                        T.pOut = T.pOut(:);
                    end
                    if isfield(data, 'tableMu')
                        T.Mu = toDoubleArray(app, data.tableMu);
                        T.Mu = T.Mu(:);
                    end
                else
                    T.ID = (1:nInt)';
                end

                T.Start  = starts(:);
                T.End    = ends(:);
                T.Lambda = lams(:);
                T.p      = ps(:);

                assignin('base', 'intervalTable', T);
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'success', ...
                    'message', sprintf('Exported %d intervals as "intervalTable"', nInt)));
            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Export intervals: ' ME.message]));
            end
        end

        % ---- Load interval table - send table list to JS -------------------
        function handleLoadIntervals(app)
            if ~app.DataLoaded
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', 'Load data first'));
                return;
            end

            try
                vars = evalin('base', 'whos');
                tableList = {};
                for i = 1:length(vars)
                    if strcmp(vars(i).class, 'table')
                        tableList{end+1} = struct(...
                            'name', vars(i).name, ...
                            'size', sprintf('%dx%d', vars(i).size(1), vars(i).size(2)), ...
                            'rows', vars(i).size(1)); %#ok<AGROW>
                    end
                end

                if isempty(tableList)
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', 'No tables in workspace'));
                    return;
                end

                % Send table list to JavaScript for modal display
                sendResponse(app, 'showTableList', struct('tables', {tableList}));

            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Error: ' ME.message]));
            end
        end

        % ---- Load a specific interval table by name -----------------------
        function handleLoadIntervalTable(app, data)
            try
                tableName = data.tableName;

                if ~isvarname(tableName)
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', ['Invalid table name: "' tableName '"']));
                    return;
                end

                if ~evalin('base', ['exist(''' tableName ''', ''var'')'])
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', ['Table "' tableName '" not found']));
                    return;
                end

                T = evalin('base', tableName);
                if ~istable(T) || ~ismember('Start', T.Properties.VariableNames) || ...
                        ~ismember('End', T.Properties.VariableNames)
                    sendResponse(app, 'statusUpdate', ...
                        struct('type', 'error', 'message', 'Table must have Start and End columns'));
                    return;
                end

                % Convert table to struct arrays for HTML
                hasPerSignal = ismember('SignalIndex', T.Properties.VariableNames);
                hasLam = ismember('Lambda', T.Properties.VariableNames);
                hasP   = ismember('p', T.Properties.VariableNames);

                nRows = height(T);
                nCh = size(app.OriginalData, 2);

                intervals = cell(1, nRows);
                for i = 1:nRows
                    intv = struct();
                    s = round(double(T.Start(i)));
                    e = round(double(T.End(i)));
                    s = max(1, min(nCh, s));
                    e = max(1, min(nCh, e));
                    if e < s, tmp=s; s=e; e=tmp; end

                    intv.startIdx = s;
                    intv.endIdx = e;
                    intv.lambda = 1000;
                    intv.p = 0.01;

                    if hasLam, intv.lambda = double(T.Lambda(i)); end
                    if hasP,   intv.p = double(T.p(i)); end

                    if hasPerSignal
                        intv.signalIndex = double(T.SignalIndex(i));
                    end
                    if ismember('LambdaOut', T.Properties.VariableNames)
                        intv.lambdaOut = double(T.LambdaOut(i));
                    end
                    if ismember('pOut', T.Properties.VariableNames)
                        intv.pOut = double(T.pOut(i));
                    end
                    if ismember('Mu', T.Properties.VariableNames)
                        intv.mu = double(T.Mu(i));
                    end

                    intervals{i} = intv;
                end

                payload = struct();
                payload.intervals = intervals;
                payload.perSignalFormat = hasPerSignal;
                payload.statusType = 'success';
                payload.statusMessage = sprintf('Loaded %d intervals from "%s"', nRows, tableName);

                sendResponse(app, 'intervalsLoaded', payload);
            catch ME
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', ['Load intervals: ' ME.message]));
            end
        end

    end

    % ====================================================================
    %  COMMUNICATION HELPERS
    % ====================================================================
    methods (Access = private)

        function sendResponse(app, responseName, payload)
            if app.IsClosed, return; end
            try
                app.UIUpdateCounter = app.UIUpdateCounter + 1;
                out = struct();
                out.response  = responseName;
                out.payload   = payload;
                out.timestamp = app.UIUpdateCounter;
                app.HTMLComponent.Data = out;
            catch
            end
        end

    end

    % ====================================================================
    %  INTERVAL BUILDING HELPER
    % ====================================================================
    methods (Access = private)

        function [intervals, pVals, lambdasAsym] = buildIntervalsFromData(app, data) %#ok<INUSL>
            intervals   = [];
            pVals       = [];
            lambdasAsym = [];

            if ~isfield(data, 'intervals') || isempty(data.intervals)
                return;
            end

            raw = data.intervals;

            if iscell(raw)
                nInt = numel(raw);
                intervals   = zeros(nInt, 2);
                pVals       = zeros(nInt, 1);
                lambdasAsym = zeros(nInt, 1);
                for i = 1:nInt
                    item = raw{i};
                    if isstruct(item)
                        intervals(i, :)  = [double(item.startIdx), double(item.endIdx)];
                        pVals(i)         = double(item.p);
                        lambdasAsym(i)   = double(item.lambda);
                    end
                end
            elseif isstruct(raw)
                if isfield(raw, 'startIdx')
                    % Could be single interval or struct array
                    if numel(raw) == 1
                        % Single interval
                        intervals   = [double(raw.startIdx), double(raw.endIdx)];
                        pVals       = double(raw.p);
                        lambdasAsym = double(raw.lambda);
                    else
                        % Struct array - multiple intervals
                        nInt = numel(raw);
                        intervals   = [[raw.startIdx]', [raw.endIdx]'];
                        pVals       = [raw.p]';
                        lambdasAsym = [raw.lambda]';
                    end
                else
                    % Multiple intervals with numeric field names
                    fns = fieldnames(raw);
                    nInt = numel(fns);
                    intervals   = zeros(nInt, 2);
                    pVals       = zeros(nInt, 1);
                    lambdasAsym = zeros(nInt, 1);
                    for i = 1:nInt
                        item = raw.(fns{i});
                        if isstruct(item)
                            intervals(i, :)  = [double(item.startIdx), double(item.endIdx)];
                            pVals(i)         = double(item.p);
                            lambdasAsym(i)   = double(item.lambda);
                        end
                    end
                end
            end
        end

    end

    % ====================================================================
    %  UTILITY
    % ====================================================================
    methods (Access = private)

        function val = safeDouble(~, data, fieldName, defaultVal)
            if isfield(data, fieldName)
                raw = data.(fieldName);
                if isnumeric(raw)
                    val = double(raw);
                else
                    try
                        val = str2double(char(raw));
                        if isnan(val), val = defaultVal; end
                    catch
                        val = defaultVal;
                    end
                end
            else
                val = defaultVal;
            end
        end

        function arr = toDoubleArray(~, raw)
            if isnumeric(raw)
                arr = double(raw(:)');
            elseif iscell(raw)
                arr = cellfun(@double, raw);
            elseif isstruct(raw)
                fns = fieldnames(raw);
                arr = zeros(1, numel(fns));
                for i = 1:numel(fns)
                    arr(i) = double(raw.(fns{i}));
                end
            else
                arr = [];
            end
        end

    end

    % ====================================================================
    %  PUBLIC API
    % ====================================================================
    methods (Access = public)

        function data = getData(app)
            if ~isempty(app.CorrectedData)
                data = app.CorrectedData;
            else
                data = app.OriginalData;
            end
        end

        function setInputData(app, data)
            if ~isreal(data), data = real(data); end
            [isValid, msg] = app.Validator.validateData(data);
            if ~isValid
                sendResponse(app, 'statusUpdate', ...
                    struct('type', 'error', 'message', msg));
                return;
            end

            app.OriginalData = double(data);
            [nRows, nCols] = size(data);
            app.Wavelength = (1:nCols)';
            app.DataLoaded = true;
            app.LoadedVarName = 'external';
            app.CorrectedData = [];
            app.BaselineData = [];
            app.WeightsData = [];

            payload = struct();
            payload.wavelength = app.Wavelength';
            payload.meanSpectrum = mean(app.OriginalData, 1);
            payload.nSamples = nRows;
            payload.nChannels = nCols;
            payload.varName = 'external';
            payload.statusType = 'success';
            payload.statusMessage = sprintf('Loaded: %dx%d matrix', nRows, nCols);

            if nRows >= 1
                payload.firstSpectrum = app.OriginalData(1, :);
            end
            if nRows <= 500 && nRows * nCols <= 500000
                allData = cell(1, nRows);
                for i = 1:nRows
                    allData{i} = app.OriginalData(i, :);
                end
                payload.allSpectra = allData;
            end

            sendResponse(app, 'dataLoaded', payload);
        end

    end

    % ====================================================================
    %  DEMO DATA HANDLER
    % ====================================================================
    methods (Access = private)

        function handleCreateDemoData(app, data)
            % Create demo data in MATLAB workspace for tutorial.
            % JavaScript sends spectra as array of arrays [sample][channel].

            if iscell(data.spectra)
                nSamples = numel(data.spectra);
                spectra = zeros(nSamples, numel(data.spectra{1}));
                for i = 1:nSamples
                    spectra(i, :) = data.spectra{i};
                end
            else
                spectra = data.spectra;
            end

            if iscell(data.wavelength)
                wl = cell2mat(data.wavelength(:))';
            else
                wl = data.wavelength(:)';
            end

            assignin('base', 'tutorialData', spectra);
            assignin('base', 'tutorialWavelength', wl);

            sendResponse(app, 'statusUpdate', struct(...
                'type', 'success', ...
                'message', 'Demo data "tutorialData" created in workspace'));
        end

    end

end
