classdef CropBackground < matlab.apps.AppBase
% CROPBACKGROUND Crop background pixels from a 3D image cube.
%
%   Interactive GUI for removing background pixels from hyperspectral or
%   multichannel image cubes. Supports three cropping methods:
%     - Threshold: min/max intensity sliders with live preview
%     - Selection: draw rectangle, polygon, or freehand regions
%     - Auto-Detect: automatic particle detection (Otsu/adaptive + watershed)
%
%   Preprocessing transforms (square, log, sqrt) can be applied before
%   any method. Returns retained pixels as a 2D matrix.
%
%   Uses the AppBase + uihtml architecture: HTML/CSS/JS frontend with
%   MATLAB backend.
%
%   Usage:
%       app = CropBackground();
%       app = CropBackground(cube);
%       app = CropBackground(struct('cube', cube));
%
%   Inputs:
%       cube or dataStruct.cube - 3D numeric array [rows x cols x channels]
%
%   Public API:
%       app.setInputData(s) - Load data programmatically
%       app.getData()       - Retrieve cropped matrix and indices
%
%   See also: BackgroundCropper, ParticleDetector, DataValidator
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.1

    % ------------------------------------------------------------------ %
    %  Properties
    % ------------------------------------------------------------------ %
    properties (Access = public)
        Figure              matlab.ui.Figure
        HTMLComponent       matlab.ui.control.HTML
    end

    properties (Access = public)
        ImageCube           double = []
        GlobalIntensity     double = []
        CroppedMatrix       double = []
        RetainedIdx         double = []
        DiscardedIdx        double = []
        MinThreshold        double = 0
        MaxThreshold        double = 0
        DataLoaded          logical = false
        IsClosed            logical = false
    end

    properties (Access = private)
        Cropper
        Validator
        Detector
        ImageRows           double = 0
        ImageCols           double = 0
        ImageChannels       double = 0
        TimestampCounter    double = 0
        RawIntensity        double = []
        PreprocessMethod    char   = 'none'
        SelectionMode       char   = 'threshold'
        SelectedBands       double = []   % empty = all bands
    end

    % ------------------------------------------------------------------ %
    %  Component initialisation (AppBase pattern)
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function createComponents(app)
            screenSize = get(0, 'ScreenSize');
            figWidth  = min(1200, screenSize(3) * 0.8);
            figHeight = min(800,  screenSize(4) * 0.85);

            app.Figure = uifigure( ...
                'Visible', 'off', ...
                'Name', 'Crop Background', ...
                'Position', [100 100 figWidth figHeight], ...
                'Color', [0.91 0.92 0.93], ...
                'CloseRequestFcn', createCallbackFcn(app, @closeRequest, true), ...
                'AutoResizeChildren', 'off', ...
                'SizeChangedFcn', createCallbackFcn(app, @figureResized, true));

            htmlPath = fullfile(fileparts(mfilename('fullpath')), ...
                'ui', 'crop_background_ui.html');

            app.HTMLComponent = uihtml(app.Figure, ...
                'HTMLSource', htmlPath, ...
                'Position', [1 1 figWidth figHeight], ...
                'DataChangedFcn', createCallbackFcn(app, @HTMLDataChanged, true));

            app.Figure.Visible = 'on';
        end

        function figureResized(app, ~)
            if app.IsClosed, return; end
            pos = app.Figure.Position;
            app.HTMLComponent.Position = [1 1 pos(3) pos(4)];
        end

        function initializeBusinessLogic(app)
            app.Cropper   = BackgroundCropper();
            app.Validator  = DataValidator();
            app.Detector   = ParticleDetector();
        end
    end

    % ------------------------------------------------------------------ %
    %  Constructor / startup
    % ------------------------------------------------------------------ %
    methods (Access = public)

        function app = CropBackground(dataStruct)
            modulePath = fileparts(mfilename('fullpath'));
            addpath(fullfile(modulePath, 'business_logic'));

            app.createComponents();
            app.initializeBusinessLogic();
            registerApp(app, app.Figure);

            if nargin >= 1 && ~isempty(dataStruct)
                if isnumeric(dataStruct)
                    dataStruct = struct('cube', dataStruct);
                end
            end

            runStartupFcn(app, @(~) app.onStartup());

            if nargin >= 1 && ~isempty(dataStruct)
                app.setInputData(dataStruct);
            end

            if nargout == 0
                clear app
            end
        end
    end

    % ------------------------------------------------------------------ %
    %  Startup
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function onStartup(app)
            movegui(app.Figure, 'center');
            figure(app.Figure);
            app.updateHTMLStatus('Ready. Load a 3-D image cube to begin.');
        end
    end

    methods (Access = public)
        function delete(app)
            app.IsClosed = true;
            delete(app.Figure);
        end
    end

    % ------------------------------------------------------------------ %
    %  Public API
    % ------------------------------------------------------------------ %
    methods (Access = public)

        function setInputData(app, dataStruct)
        % setInputData  Load a 3-D image cube into the GUI.
        %   dataStruct must contain field: cube (3-D numeric array)
            if ~isstruct(dataStruct) || ~isfield(dataStruct, 'cube')
                app.updateHTMLStatus('Invalid input: expected struct with "cube" field.');
                return;
            end

            cube = dataStruct.cube;

            [valid, msg] = DataValidator.validate3DArray(cube);
            if ~valid
                app.updateHTMLStatus(msg);
                return;
            end

            app.ImageCube = double(cube);
            [app.ImageRows, app.ImageCols, app.ImageChannels] = size(app.ImageCube);

            % Compute global intensity and store raw copy for preprocessing
            app.RawIntensity    = app.Cropper.computeIntensity(app.ImageCube);
            app.GlobalIntensity = app.RawIntensity;
            app.PreprocessMethod = 'none';
            app.SelectionMode    = 'threshold';
            app.SelectedBands    = [];

            app.DataLoaded   = true;
            app.MinThreshold = min(app.GlobalIntensity(:));
            app.MaxThreshold = max(app.GlobalIntensity(:));
            app.CroppedMatrix  = [];
            app.RetainedIdx    = [];
            app.DiscardedIdx   = [];

            app.sendDataToHTML();
        end

        function result = getData(app)
        % getData  Return the cropped data.
            result = struct();
            if ~isempty(app.CroppedMatrix)
                result.croppedMatrix = app.CroppedMatrix;
                result.retainedIdx   = app.RetainedIdx;
                result.discardedIdx  = app.DiscardedIdx;
            else
                result.croppedMatrix = [];
                result.retainedIdx   = [];
                result.discardedIdx  = [];
            end
            result.imageCube       = app.ImageCube;
            result.globalIntensity = app.GlobalIntensity;
        end

    end

    % ------------------------------------------------------------------ %
    %  HTML event dispatcher
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function HTMLDataChanged(app, ~, ~)
            eventData = app.HTMLComponent.Data;
            if ~isstruct(eventData)
                return;
            end
            if isfield(eventData, 'source') && strcmp(eventData.source, 'matlab')
                return;
            end
            if ~isfield(eventData, 'action')
                return;
            end

            try
                switch eventData.action
                    case 'load_data'
                        app.handleLoadData();
                    case 'load_variable'
                        app.handleLoadVariable(eventData);
                    case 'update_threshold'
                        app.handleUpdateThreshold(eventData);
                    case 'apply'
                        app.handleApply(eventData);
                    case 'export'
                        app.handleExport();
                    case 'preprocess'
                        app.handlePreprocess(eventData);
                    case 'apply_mask'
                        app.handleApplyMask(eventData);
                    case 'auto_detect'
                        app.handleAutoDetect(eventData);
                    case 'get_spectrum'
                        app.handleGetSpectrum();
                    case 'select_bands'
                        app.handleSelectBands(eventData);
                    otherwise
                        app.updateHTMLStatus(['Unknown action: ' eventData.action]);
                end
            catch ex
                app.updateHTMLStatus(['Error: ' ex.message]);
            end
        end
    end

    % ------------------------------------------------------------------ %
    %  Action handlers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function handleLoadData(app)
        % handleLoadData  Send 3-D workspace variables to the HTML modal.
            vars = evalin('base', 'whos');
            variables = {};
            for k = 1:numel(vars)
                v = vars(k);
                numericTypes = {'double','single','uint8','uint16','int16','int32','uint32','int64','uint64'};
                if numel(v.size) >= 3 && v.size(3) > 1 && any(strcmp(v.class, numericTypes))
                    variables{end+1} = struct( ...
                        'name', v.name, ...
                        'class', v.class, ...
                        'size', sprintf('%dx%dx%d', v.size(1), v.size(2), v.size(3))); %#ok<AGROW>
                end
            end

            if isempty(variables)
                app.updateHTMLStatus('No 3-D numeric arrays found in workspace.');
                return;
            end

            payload = struct();
            payload.kind = 'showVarList';
            payload.variables = {variables};
            app.sendUIUpdate(payload);
        end

        function handleLoadVariable(app, evt)
        % handleLoadVariable  Load the selected variable from workspace.
            if ~isfield(evt, 'varName')
                app.updateHTMLStatus('Missing variable name.');
                return;
            end

            varName = char(string(evt.varName));
            if ~evalin('base', sprintf('exist(''%s'', ''var'')', varName))
                app.updateHTMLStatus(sprintf('Variable "%s" not found.', varName));
                return;
            end

            cube = evalin('base', varName);
            ds = struct('cube', cube);
            app.setInputData(ds);
        end


        function handleUpdateThreshold(app, evt)
        % handleUpdateThreshold  Update min/max thresholds and recompute preview.
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            minT = evt.minThreshold;
            maxT = evt.maxThreshold;
            app.SelectionMode = 'threshold';

            [valid, msg] = DataValidator.validateThresholds(minT, maxT);
            if ~valid
                app.updateHTMLStatus(msg);
                return;
            end

            app.MinThreshold = minT;
            app.MaxThreshold = maxT;

            [retIdx, disIdx] = app.Cropper.applyThreshold( ...
                app.GlobalIntensity, minT, maxT);

            % Build retained mask for the heatmap overlay
            retainedMask = true(app.ImageRows, app.ImageCols);
            retainedMask(disIdx) = false;

            payload = struct();
            payload.kind           = 'threshold_result';
            payload.retainedCount  = numel(retIdx);
            payload.discardedCount = numel(disIdx);
            payload.retainedMask   = reshape(retainedMask', 1, []);
            payload.minThreshold   = minT;
            payload.maxThreshold   = maxT;
            payload.statusMessage  = sprintf('Preview: %d retained, %d discarded.', numel(retIdx), numel(disIdx));
            app.sendUIUpdate(payload);
        end

        function handleApply(app, evt)
        % handleApply  Finalise the cropping selection (mode-aware).
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            % Read mode and thresholds from JS (single message, no race)
            if nargin >= 2 && isstruct(evt)
                if isfield(evt, 'mode')
                    app.SelectionMode = char(string(evt.mode));
                end
                if isfield(evt, 'minThreshold')
                    app.MinThreshold = double(evt.minThreshold);
                end
                if isfield(evt, 'maxThreshold')
                    app.MaxThreshold = double(evt.maxThreshold);
                end
            end

            % In threshold mode, recompute from current thresholds
            if strcmp(app.SelectionMode, 'threshold')
                [app.RetainedIdx, app.DiscardedIdx] = ...
                    app.Cropper.applyThreshold(app.GlobalIntensity, ...
                        app.MinThreshold, app.MaxThreshold);
            end

            % For selection / autodetect, RetainedIdx was set by the preview
            if isempty(app.RetainedIdx)
                app.updateHTMLStatus('No selection made. Adjust parameters first.');
                return;
            end

            app.CroppedMatrix = app.Cropper.extractPixels( ...
                app.ImageCube, app.RetainedIdx);

            payload = struct();
            payload.kind          = 'apply_result';
            payload.retainedCount = numel(app.RetainedIdx);
            payload.croppedRows   = size(app.CroppedMatrix, 1);
            payload.croppedCols   = size(app.CroppedMatrix, 2);
            payload.statusMessage = sprintf('Applied. Cropped matrix: %d x %d.', ...
                size(app.CroppedMatrix, 1), size(app.CroppedMatrix, 2));
            app.sendUIUpdate(payload);
        end

        function handleExport(app)
        % handleExport  Export cropped matrix + indices to the base workspace.
            if isempty(app.CroppedMatrix)
                app.updateHTMLStatus('Nothing to export. Apply thresholds first.');
                return;
            end

            assignin('base', 'croppedMatrix', app.CroppedMatrix);
            assignin('base', 'retainedIdx',   app.RetainedIdx);
            assignin('base', 'discardedIdx',  app.DiscardedIdx);

            % Single message — no overwrite risk
            app.updateHTMLStatus('Exported croppedMatrix, retainedIdx, discardedIdx to workspace.');
        end


        function handlePreprocess(app, evt)
        % handlePreprocess  Apply intensity transform and resend data.
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            method = 'none';
            if isfield(evt, 'transform')
                method = char(string(evt.transform));
            end

            app.PreprocessMethod = method;
            app.GlobalIntensity  = app.Cropper.preprocessIntensity( ...
                app.RawIntensity, method);

            app.MinThreshold   = min(app.GlobalIntensity(:));
            app.MaxThreshold   = max(app.GlobalIntensity(:));
            app.CroppedMatrix  = [];
            app.RetainedIdx    = [];
            app.DiscardedIdx   = [];

            app.sendDataToHTML(sprintf('Preprocessing: %s applied.', method));
        end

        function handleApplyMask(app, evt)
        % handleApplyMask  Receive a binary mask from the JS drawing tools.
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            if ~isfield(evt, 'mask')
                app.updateHTMLStatus('No mask received.');
                return;
            end

            % JS sends row-major flat array; reshape to MATLAB column-major
            maskFlat = logical(evt.mask(:));
            nPix = app.ImageRows * app.ImageCols;
            if numel(maskFlat) ~= nPix
                app.updateHTMLStatus('Mask size mismatch.');
                return;
            end

            % Row-major (JS) to column-major (MATLAB)
            maskMatrix = reshape(maskFlat, app.ImageCols, app.ImageRows)';

            app.SelectionMode = 'selection';
            app.RetainedIdx   = find(maskMatrix(:));
            app.DiscardedIdx  = find(~maskMatrix(:));

            % Send preview back (column-major to row-major for JS)
            retainedMaskRM = reshape(maskMatrix', 1, []);

            payload = struct();
            payload.kind           = 'threshold_result';
            payload.retainedCount  = numel(app.RetainedIdx);
            payload.discardedCount = numel(app.DiscardedIdx);
            payload.retainedMask   = retainedMaskRM;
            payload.minThreshold   = app.MinThreshold;
            payload.maxThreshold   = app.MaxThreshold;
            payload.statusMessage  = sprintf('Selection: %d retained, %d discarded.', ...
                numel(app.RetainedIdx), numel(app.DiscardedIdx));
            app.sendUIUpdate(payload);
        end

        function handleAutoDetect(app, evt)
        % handleAutoDetect  Run automatic particle/cell detection.
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            params = struct();
            if isfield(evt, 'method'),      params.method      = char(string(evt.method)); end
            if isfield(evt, 'sensitivity'),  params.sensitivity = double(evt.sensitivity);  end
            if isfield(evt, 'minArea'),      params.minArea     = double(evt.minArea);      end
            if isfield(evt, 'watershed'),    params.watershed   = logical(evt.watershed);   end
            if isfield(evt, 'smoothSigma'),  params.smoothSigma = double(evt.smoothSigma);  end

            try
                [mask, numObj] = ParticleDetector.detect( ...
                    app.GlobalIntensity, params);
            catch ex
                app.updateHTMLStatus(['Auto-detect: ' ex.message]);
                return;
            end

            app.SelectionMode = 'autodetect';
            app.RetainedIdx   = find(mask(:));
            app.DiscardedIdx  = find(~mask(:));

            retainedMaskRM = reshape(mask', 1, []);

            payload = struct();
            payload.kind           = 'threshold_result';
            payload.retainedCount  = numel(app.RetainedIdx);
            payload.discardedCount = numel(app.DiscardedIdx);
            payload.retainedMask   = retainedMaskRM;
            payload.minThreshold   = app.MinThreshold;
            payload.maxThreshold   = app.MaxThreshold;
            payload.statusMessage  = sprintf('Detected %d objects (%d pixels retained).', ...
                numObj, numel(app.RetainedIdx));
            app.sendUIUpdate(payload);
        end

        function handleGetSpectrum(app)
        % handleGetSpectrum  Compute and send the mean spectrum across all pixels.
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            % Reshape cube to [nPixels x nChannels], compute mean across pixels
            nPix = app.ImageRows * app.ImageCols;
            reshaped = reshape(app.ImageCube, nPix, app.ImageChannels);
            meanSpec = mean(reshaped, 1);

            payload = struct();
            payload.kind          = 'mean_spectrum';
            payload.spectrum      = meanSpec;
            payload.nBands        = app.ImageChannels;
            payload.statusMessage = sprintf('Mean spectrum: %d bands.', app.ImageChannels);
            app.sendUIUpdate(payload);
        end

        function handleSelectBands(app, evt)
        % handleSelectBands  Recompute intensity using only selected spectral bands.
            if ~app.DataLoaded
                app.updateHTMLStatus('No data loaded.');
                return;
            end

            if isfield(evt, 'bands') && ~isempty(evt.bands)
                bands = round(double(evt.bands(:)));
                bands = bands(bands >= 1 & bands <= app.ImageChannels);
                app.SelectedBands = bands;
                app.RawIntensity = sum(app.ImageCube(:,:,bands), 3);
            else
                % Empty = use all bands
                app.SelectedBands = [];
                app.RawIntensity = app.Cropper.computeIntensity(app.ImageCube);
            end

            % Reapply current preprocessing on the new band selection
            app.GlobalIntensity = app.Cropper.preprocessIntensity( ...
                app.RawIntensity, app.PreprocessMethod);

            app.MinThreshold  = min(app.GlobalIntensity(:));
            app.MaxThreshold  = max(app.GlobalIntensity(:));
            app.CroppedMatrix = [];
            app.RetainedIdx   = [];
            app.DiscardedIdx  = [];

            nSel = numel(app.SelectedBands);
            if nSel == 0
                msg = 'Using all bands.';
            else
                msg = sprintf('Using %d of %d bands.', nSel, app.ImageChannels);
            end
            app.sendDataToHTML(msg);
        end
    end

    % ------------------------------------------------------------------ %
    %  Data transmission helpers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function sendDataToHTML(app, statusMsg)
        % sendDataToHTML  Send intensity map data and distribution to the HTML UI.

            if nargin < 2, statusMsg = ''; end

            % Build sorted distribution (down-sampled for performance)
            sorted = sort(app.GlobalIntensity(:));
            nPix   = numel(sorted);
            if nPix > 5000
                idx      = round(linspace(1, nPix, 5000));
                sortedDS = sorted(idx);
            else
                sortedDS = sorted;
            end

            % Normalise intensity to 0-1 for the heatmap
            minI = min(app.GlobalIntensity(:));
            maxI = max(app.GlobalIntensity(:));
            if maxI == minI
                normImg = zeros(app.ImageRows, app.ImageCols);
            else
                normImg = (app.GlobalIntensity - minI) / (maxI - minI);
            end

            % Initial retained mask: all retained
            retainedMask = true(1, nPix);

            payload = struct();
            payload.kind          = 'data_loaded';
            payload.rows          = app.ImageRows;
            payload.cols          = app.ImageCols;
            payload.channels      = app.ImageChannels;
            payload.intensityNorm = reshape(normImg', 1, []);
            payload.sortedDist    = sortedDS(:)';
            payload.minVal        = minI;
            payload.maxVal        = maxI;
            payload.totalPixels   = nPix;
            payload.retainedMask  = retainedMask;
            payload.statusMessage = statusMsg;
            app.sendUIUpdate(payload);
        end
    end

    % ------------------------------------------------------------------ %
    %  UI communication helpers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function sendUIUpdate(app, payload)
        % sendUIUpdate  Send a payload to the HTML component with source and timestamp.
            payload.source    = 'matlab';
            app.TimestampCounter = app.TimestampCounter + 1;
            payload.timestamp = app.TimestampCounter;
            app.HTMLComponent.Data = payload;
        end

        function updateHTMLStatus(app, msg)
        % updateHTMLStatus  Send a status message to the HTML UI.
            payload = struct();
            payload.kind    = 'status';
            payload.message = msg;
            app.sendUIUpdate(payload);
        end
    end

    % ------------------------------------------------------------------ %
    %  Close request
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function closeRequest(app, ~, ~)
            app.IsClosed = true;
            delete(app.Figure);
        end
    end
end
