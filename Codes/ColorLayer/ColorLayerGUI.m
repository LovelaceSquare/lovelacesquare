classdef ColorLayerGUI < matlab.apps.AppBase
    % ColorLayerGUI. Interactive false-color composite image editor
    %
    %   Creates composite images from multiple 2D maps using additive
    %   color blending (6 channels: R, G, B, Y, C, M) with per-channel
    %   brightness/saturation controls.
    %
    %   Uses uihtml + AppBase architecture with JS-side compositing for
    %   real-time slider preview.
    %
    %   Usage:
    %       app = ColorLayerGUI();
    %       app = ColorLayerGUI(mdis);
    %       app = ColorLayerGUI(mdis, Original);
    %
    %   Inputs:
    %       mdis     - Cell array of 2D maps {map1, map2, ..., mapN},
    %                  3D numeric array (H x W x N), or 2D matrix (H x W).
    %                  Any numeric type (converted to double internally).
    %       Original - Background image (H x W) or (H x W x 3)
    %
    % Author: Adrian Gomez-Sanchez
    % Date Created: 2026-03-08
    % License: MIT
    % Reviewed by Lovelace's Square: Yes
    % Version: v 1.0

    properties (Access = public)
        UIFigure        matlab.ui.Figure
        HTMLComponent   matlab.ui.control.HTML
    end

    properties (Access = private)
        MapsData        cell   = {}       % mdis cell array (full resolution)
        OriginalBg      double = []       % Background image
        Channels        struct            % 6-channel settings
        DataLoaded      logical = false
        IsClosed        logical = false
        UIUpdateCounter double  = 0
        SliderLimits    struct            % Custom slider ranges
    end

    % =====================================================================
    %  CONSTRUCTOR / DESTRUCTOR
    % =====================================================================

    methods (Access = public)
        function app = ColorLayerGUI(mdis, Original)
            % Initialize default channel settings
            app.initChannels();
            app.SliderLimits = struct( ...
                'brightnessMin', 0.1, 'brightnessMax', 3.0, ...
                'saturationMin', 0.1, 'saturationMax', 2.0);

            % Create UI components
            createComponents(app);
            registerApp(app, app.UIFigure);

            % Store input data if provided
            if nargin >= 1 && ~isempty(mdis)
                if isnumeric(mdis)
                    mdis = double(mdis);
                    if ndims(mdis) == 3
                        nCh = size(mdis, 3);
                        maps = cell(1, nCh);
                        for k = 1:nCh
                            maps{1, k} = mdis(:, :, k);
                        end
                        mdis = maps;
                    elseif ismatrix(mdis)
                        mdis = {mdis};
                    end
                end
                app.MapsData = mdis;
                app.DataLoaded = true;
            end
            if nargin >= 2 && ~isempty(Original)
                app.OriginalBg = Original;
            end

            runStartupFcn(app, @startupFcn);

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            app.IsClosed = true;
            delete(app.UIFigure);
        end
    end

    % =====================================================================
    %  UI CREATION
    % =====================================================================

    methods (Access = private)
        function createComponents(app)
            screenSize = get(0, 'ScreenSize');
            figWidth  = min(1500, screenSize(3) * 0.8);
            figHeight = min(900,  screenSize(4) * 0.85);

            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100, 100, figWidth, figHeight];
            app.UIFigure.Name = 'ColorLayer Editor';
            app.UIFigure.Color = [0.91, 0.92, 0.93];
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @closeRequest, true);
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @figureResized, true);

            % Set icon if available
            modulePath = fileparts(mfilename('fullpath'));
            iconPath = fullfile(modulePath, 'favicon.png');
            if isfile(iconPath)
                app.UIFigure.Icon = iconPath;
            end

            % uihtml fills entire figure
            app.HTMLComponent = uihtml(app.UIFigure);
            app.HTMLComponent.Position = [1, 1, figWidth, figHeight];
            app.HTMLComponent.HTMLSource = fullfile(modulePath, 'ui', 'colorlayer_ui.html');
            app.HTMLComponent.DataChangedFcn = createCallbackFcn(app, @HTMLDataChanged, true);

            app.UIFigure.Visible = 'on';
        end

        function startupFcn(app)
            movegui(app.UIFigure, 'center');

            % If data was passed to constructor, send it to JS
            if app.DataLoaded
                autoEnableChannels(app);
                sendMapsToJS(app);
            else
                sendResponse(app, struct('type', 'ready', 'message', 'Ready'));
            end
        end

        function closeRequest(app, ~)
            app.IsClosed = true;
            try, uiresume(app.UIFigure); catch, end
            delete(app);
        end

        function figureResized(app, ~)
            if app.IsClosed, return; end
            pos = app.UIFigure.Position;
            app.HTMLComponent.Position = [1, 1, pos(3), pos(4)];
        end
    end

    % =====================================================================
    %  HTML COMMUNICATION
    % =====================================================================

    methods (Access = private)
        function HTMLDataChanged(app, ~)
            if app.IsClosed, return; end
            try
                data = app.HTMLComponent.Data;
                if ~isstruct(data) || ~isfield(data, 'action'), return; end
                if isfield(data, 'source') && strcmp(data.source, 'matlab'), return; end

                switch data.action
                    case 'loadData',           handleLoadData(app, data);
                    case 'loadVariable',       handleLoadVariable(app, data);
                    case 'loadFromFile',       handleLoadFromFile(app, data);
                    case 'exportImage',        handleExportImage(app, data);
                    case 'savePreset',         handleSavePreset(app, data);
                    case 'loadPreset',         handleLoadPreset(app, data);
                    case 'exportToWorkspace',  handleExportToWorkspace(app, data);
                    case 'updateLimits',       handleUpdateLimits(app, data);
                    otherwise
                        sendResponse(app, struct('type', 'error', ...
                            'message', sprintf('Unknown action: %s', data.action)));
                end
            catch ME
                sendResponse(app, struct('type', 'error', 'message', ME.message));
            end
        end

        function sendResponse(app, response)
            if app.IsClosed, return; end
            response.source = 'matlab';
            response.timestamp = posixtime(datetime('now'));
            app.UIUpdateCounter = app.UIUpdateCounter + 1;
            response.counter = app.UIUpdateCounter;
            app.HTMLComponent.Data = response;
        end
    end

    % =====================================================================
    %  ACTION HANDLERS
    % =====================================================================

    methods (Access = private)
        function handleLoadData(app, ~)
            % List workspace cell arrays and numeric matrices/images
            vars = evalin('base', 'whos');
            varList = {};
            numericTypes = {'double','single','uint8','uint16','int16','int32'};
            for i = 1:length(vars)
                v = vars(i);
                isCell = strcmp(v.class, 'cell');
                isNumeric = ismember(v.class, numericTypes);
                is2D = isNumeric && length(v.size) == 2;
                is3D = isNumeric && length(v.size) == 3;
                if isCell || is2D || is3D
                    sizeStr = strjoin(arrayfun(@num2str, v.size, 'UniformOutput', false), 'x');
                    varList{end+1} = struct( ...
                        'name', v.name, ...
                        'class', v.class, ...
                        'size', sizeStr); %#ok<AGROW>
                end
            end
            sendResponse(app, struct('type', 'showVarList', 'variables', {varList}));
        end

        function handleLoadVariable(app, data)
            try
                mdisName = data.mdisVar;
                mdisVar = evalin('base', mdisName);

                % Convert numeric arrays to cell array of maps
                if isnumeric(mdisVar)
                    mdisVar = double(mdisVar);
                    if ndims(mdisVar) == 3
                        nCh = size(mdisVar, 3);
                        maps = cell(1, nCh);
                        for k = 1:nCh
                            maps{1, k} = mdisVar(:, :, k);
                        end
                        mdisVar = maps;
                    elseif ismatrix(mdisVar)
                        mdisVar = {mdisVar};
                    end
                end

                [isValid, msg] = utils.validateMaps(mdisVar);
                if ~isValid
                    sendResponse(app, struct('type', 'error', 'message', msg));
                    return;
                end

                app.MapsData = mdisVar;
                app.DataLoaded = true;

                % Load optional background
                if isfield(data, 'bgVar') && ~isempty(data.bgVar)
                    try
                        app.OriginalBg = evalin('base', data.bgVar);
                    catch
                        app.OriginalBg = [];
                    end
                end

                autoEnableChannels(app);
                sendMapsToJS(app);

            catch ME
                sendResponse(app, struct('type', 'error', 'message', ME.message));
            end
        end

        function handleLoadFromFile(app, ~)
            [file, filepath] = uigetfile( ...
                {'*.mat', 'MAT files'; ...
                 '*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Image files'; ...
                 '*.*', 'All files'}, ...
                'Select file');
            if file == 0
                sendResponse(app, struct('type', 'statusUpdate', ...
                    'message', 'File load cancelled'));
                return;
            end

            fullPath = fullfile(filepath, file);
            [~, ~, ext] = fileparts(file);

            try
                if strcmpi(ext, '.mat')
                    % MAT file: look for mdis cell array
                    loadedData = load(fullPath);
                    fields = fieldnames(loadedData);

                    if isfield(loadedData, 'mdis')
                        mdisVar = loadedData.mdis;
                    elseif ~isempty(fields)
                        mdisVar = loadedData.(fields{1});
                    else
                        sendResponse(app, struct('type', 'error', ...
                            'message', 'No data found in file'));
                        return;
                    end

                    [isValid, msg] = utils.validateMaps(mdisVar);
                    if ~isValid
                        sendResponse(app, struct('type', 'error', 'message', msg));
                        return;
                    end

                    app.MapsData = mdisVar;

                    if isfield(loadedData, 'Original')
                        app.OriginalBg = loadedData.Original;
                    end
                else
                    % Image file: split channels into maps
                    img = double(imread(fullPath));
                    nCh = size(img, 3);

                    if nCh == 1
                        % Grayscale: single map
                        mdisVar = {img};
                    else
                        % RGB or multi-channel: one map per channel
                        mdisVar = cell(1, nCh);
                        for k = 1:nCh
                            mdisVar{1, k} = img(:, :, k);
                        end
                    end

                    [isValid, msg] = utils.validateMaps(mdisVar);
                    if ~isValid
                        sendResponse(app, struct('type', 'error', 'message', msg));
                        return;
                    end

                    app.MapsData = mdisVar;
                    app.OriginalBg = [];
                end

                app.DataLoaded = true;
                autoEnableChannels(app);
                sendMapsToJS(app);

            catch ME
                sendResponse(app, struct('type', 'error', 'message', ME.message));
            end
        end

        function handleExportImage(app, data)
            if ~app.DataLoaded
                sendResponse(app, struct('type', 'error', ...
                    'message', 'No data loaded'));
                return;
            end

            % Get channel settings from JS
            if isfield(data, 'channels')
                app.Channels = parseChannelsFromJS(app, data.channels);
            end

            options = struct();
            if isfield(data, 'normMode')
                options.normalizationMode = data.normMode;
            else
                options.normalizationMode = 'Global mat2gray';
            end
            if isfield(data, 'perMapNorm')
                options.perMapNormalize = logical(data.perMapNorm);
            else
                options.perMapNormalize = false;
            end

            % Parse gamma and contrast limits from JS channels
            gammaVals = ones(1, 6);
            contrastMins = zeros(1, 6);
            contrastMaxs = ones(1, 6);
            if isfield(data, 'channels')
                for ci = 1:min(length(data.channels), 6)
                    chJS = data.channels{ci};
                    if isstruct(chJS)
                        if isfield(chJS, 'gamma'), gammaVals(ci) = chJS.gamma; end
                        if isfield(chJS, 'contrastMin'), contrastMins(ci) = chJS.contrastMin; end
                        if isfield(chJS, 'contrastMax'), contrastMaxs(ci) = chJS.contrastMax; end
                    end
                end
            end

            % Compute composite at full resolution
            composite = core.createCompositeImage(app.MapsData, app.Channels, options);

            % Apply per-channel gamma and contrast (post-processing)
            % Note: For full-res export, we apply gamma/contrast to the
            % composite. Since createCompositeImage already produced the
            % additive blend, gamma is applied as a global post-process
            % on the composite. For precise per-channel gamma, the user
            % should tune visually and the JS preview matches the export.
            hasGammaOrContrast = false;
            for ci = 1:6
                if gammaVals(ci) ~= 1.0 || contrastMins(ci) > 0 || contrastMaxs(ci) < 1
                    hasGammaOrContrast = true;
                    break;
                end
            end
            if hasGammaOrContrast
                % Recompute with per-channel gamma/contrast at full res
                composite = computeWithGamma(app, options, gammaVals, contrastMins, contrastMaxs);
            end

            % Blend with background if present and enabled
            bgEnabled = true;
            if isfield(data, 'backgroundEnabled')
                bgEnabled = logical(data.backgroundEnabled);
            end
            if ~isempty(app.OriginalBg) && bgEnabled && isfield(data, 'bgAlpha')
                composite = core.blendWithBackground(composite, app.OriginalBg, data.bgAlpha);
            end

            % Determine export format
            exportFormat = 'png';
            if isfield(data, 'format')
                exportFormat = lower(data.format);
            end

            if strcmp(exportFormat, 'emf')
                % EMF export via print -dmeta
                [file, filepath] = uiputfile( ...
                    {'*.emf', 'Enhanced Metafile'}, ...
                    'Export as EMF');
                if file == 0
                    sendResponse(app, struct('type', 'statusUpdate', ...
                        'message', 'Export cancelled'));
                    return;
                end
                [H, W, ~] = size(composite);
                fig = figure('Visible', 'off', 'Color', 'w');
                ax = axes(fig, 'Position', [0 0 1 1]);
                image(ax, im2uint8(composite));
                axis(ax, 'image', 'off');
                fig.PaperPositionMode = 'auto';
                fig.Units = 'pixels';
                fig.Position = [100, 100, W, H];
                print(fig, fullfile(filepath, file), '-dmeta');
                close(fig);
            else
                % PNG/TIFF export via imwrite
                [file, filepath] = uiputfile( ...
                    {'*.png', 'PNG Image'; '*.tif', 'TIFF Image'}, ...
                    'Export Composite Image');
                if file == 0
                    sendResponse(app, struct('type', 'statusUpdate', ...
                        'message', 'Export cancelled'));
                    return;
                end
                imwrite(composite, fullfile(filepath, file));
            end

            sendResponse(app, struct('type', 'exportCompleted', ...
                'message', sprintf('Image saved to %s', file)));
        end

        function handleSavePreset(app, data)
            [file, filepath] = uiputfile('*.mat', 'Save Preset');
            if file == 0
                sendResponse(app, struct('type', 'statusUpdate', ...
                    'message', 'Save cancelled'));
                return;
            end

            try
                if isfield(data, 'channels')
                    app.Channels = parseChannelsFromJS(app, data.channels);
                end

                options = struct();
                if isfield(data, 'normMode')
                    options.normalizationMode = data.normMode;
                else
                    options.normalizationMode = 'Global mat2gray';
                end
                if isfield(data, 'perMapNorm')
                    options.perMapNormalize = logical(data.perMapNorm);
                else
                    options.perMapNormalize = false;
                end

                utils.savePreset(fullfile(filepath, file), app.Channels, options);
                sendResponse(app, struct('type', 'statusUpdate', ...
                    'message', 'Preset saved successfully'));
            catch ME
                sendResponse(app, struct('type', 'error', 'message', ME.message));
            end
        end

        function handleLoadPreset(app, ~)
            [file, filepath] = uigetfile('*.mat', 'Load Preset');
            if file == 0
                sendResponse(app, struct('type', 'statusUpdate', ...
                    'message', 'Load cancelled'));
                return;
            end

            try
                [channels_loaded, options_loaded] = utils.loadPreset(fullfile(filepath, file));

                app.Channels = channels_loaded;

                % Send preset data to JS
                chData = cell(1, 6);
                for i = 1:6
                    ch = channels_loaded(i);
                    satVal = ch.saturation;
                    if isinf(satVal), satVal = -1; end
                    chData{i} = struct( ...
                        'enabled', ch.enabled, ...
                        'mapIndex', ch.mapIndex, ...
                        'brightness', ch.brightness, ...
                        'saturation', satVal);
                end

                sendResponse(app, struct('type', 'presetLoaded', ...
                    'channels', {chData}, ...
                    'normMode', options_loaded.normalizationMode, ...
                    'perMapNorm', options_loaded.perMapNormalize));
            catch ME
                sendResponse(app, struct('type', 'error', 'message', ME.message));
            end
        end

        function handleExportToWorkspace(app, data)
            if ~app.DataLoaded
                sendResponse(app, struct('type', 'error', ...
                    'message', 'No data loaded'));
                return;
            end

            if isfield(data, 'channels')
                app.Channels = parseChannelsFromJS(app, data.channels);
            end

            options = struct();
            if isfield(data, 'normMode')
                options.normalizationMode = data.normMode;
            else
                options.normalizationMode = 'Global mat2gray';
            end
            if isfield(data, 'perMapNorm')
                options.perMapNormalize = logical(data.perMapNorm);
            else
                options.perMapNormalize = false;
            end

            composite = core.createCompositeImage(app.MapsData, app.Channels, options);

            if ~isempty(app.OriginalBg) && isfield(data, 'bgAlpha')
                composite = core.blendWithBackground(composite, app.OriginalBg, data.bgAlpha);
            end

            varName = 'compositeImage';
            if isfield(data, 'varName') && ~isempty(data.varName)
                varName = data.varName;
            end

            assignin('base', varName, composite);
            sendResponse(app, struct('type', 'statusUpdate', ...
                'message', sprintf('Composite exported as ''%s''', varName)));
        end

        function handleUpdateLimits(app, data)
            if isfield(data, 'brightnessMin')
                app.SliderLimits.brightnessMin = data.brightnessMin;
            end
            if isfield(data, 'brightnessMax')
                app.SliderLimits.brightnessMax = data.brightnessMax;
            end
            if isfield(data, 'saturationMin')
                app.SliderLimits.saturationMin = data.saturationMin;
            end
            if isfield(data, 'saturationMax')
                app.SliderLimits.saturationMax = data.saturationMax;
            end
            sendResponse(app, struct('type', 'limitsUpdated', ...
                'brightnessMin', app.SliderLimits.brightnessMin, ...
                'brightnessMax', app.SliderLimits.brightnessMax, ...
                'saturationMin', app.SliderLimits.saturationMin, ...
                'saturationMax', app.SliderLimits.saturationMax));
        end
    end

    % =====================================================================
    %  UTILITY METHODS
    % =====================================================================

    methods (Access = private)
        function initChannels(app)
            channelNames = {'Red', 'Green', 'Blue', 'Yellow', 'Cyan', 'Magenta'};
            channelRGB   = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [0 1 1], [1 0 1]};

            for i = 1:6
                app.Channels(i).name       = channelNames{i};
                app.Channels(i).enabled    = false;
                app.Channels(i).mapIndex   = min(i, 1);
                app.Channels(i).brightness = 1.0;
                app.Channels(i).saturation = Inf;
                app.Channels(i).lumination = 1.0;
                app.Channels(i).rgb        = channelRGB{i};
            end
        end

        function autoEnableChannels(app)
            nMaps = length(app.MapsData);
            nEnable = min(nMaps, 6);
            for i = 1:6
                if i <= nEnable
                    app.Channels(i).enabled  = true;
                    app.Channels(i).mapIndex = i;
                else
                    app.Channels(i).enabled  = false;
                    app.Channels(i).mapIndex = 1;
                end
            end
        end

        function sendMapsToJS(app)
            % Prepare maps payload for JS-side compositing
            nMaps = length(app.MapsData);

            % Get reference size
            if size(app.MapsData, 1) == 1
                refMap = app.MapsData{1,1};
            else
                refMap = app.MapsData{1};
            end
            [H, W] = size(refMap);

            % Downsample if too large (>500K total pixels per map)
            maxPixels = 500000;
            scale = 1;
            if H * W > maxPixels
                scale = sqrt(maxPixels / (H * W));
                H = round(H * scale);
                W = round(W * scale);
            end

            % Convert maps to flat row-major arrays
            mapsFlat = cell(1, nMaps);
            for i = 1:nMaps
                if size(app.MapsData, 1) == 1
                    M = double(app.MapsData{1, i});
                else
                    M = double(app.MapsData{i});
                end
                if scale < 1
                    M = utils.imresize(M, [H, W], 'bicubic');
                end
                mapsFlat{i} = reshape(M', 1, []);  % Row-major: transpose then flatten
            end

            % Prepare channel info
            chData = cell(1, 6);
            for i = 1:6
                ch = app.Channels(i);
                satVal = ch.saturation;
                if isinf(satVal), satVal = -1; end
                chData{i} = struct( ...
                    'enabled', ch.enabled, ...
                    'mapIndex', ch.mapIndex, ...
                    'brightness', ch.brightness, ...
                    'saturation', satVal, ...
                    'rgb', ch.rgb);
            end

            % Prepare background if available
            payload = struct();
            payload.type = 'dataLoaded';
            payload.maps = mapsFlat;
            payload.mapHeight = H;
            payload.mapWidth = W;
            payload.nMaps = nMaps;
            payload.channels = chData;

            if ~isempty(app.OriginalBg)
                bg = double(app.OriginalBg);
                if scale < 1
                    bg = utils.imresize(bg, [H, W], 'bicubic');
                end
                if size(bg, 3) == 1
                    bg = repmat(bg, [1, 1, 3]);
                end
                bg = utils.mat2gray(bg);
                % Flatten row-major: [R1,G1,B1, R2,G2,B2, ...] per row
                bgFlat = zeros(1, H * W * 3);
                for r = 1:H
                    for c = 1:W
                        idx = ((r-1)*W + (c-1)) * 3;
                        bgFlat(idx+1) = bg(r, c, 1);
                        bgFlat(idx+2) = bg(r, c, 2);
                        bgFlat(idx+3) = bg(r, c, 3);
                    end
                end
                payload.hasBackground = true;
                payload.background = bgFlat;
            else
                payload.hasBackground = false;
            end

            payload.message = sprintf('Loaded %d maps (%dx%d)', nMaps, H, W);
            sendResponse(app, payload);
        end

        function composite = computeWithGamma(app, options, gammaVals, contrastMins, contrastMaxs)
            % Recompute composite with per-channel gamma and contrast limits
            % Mirrors the JS compositing engine logic
            mdis = app.MapsData;
            chs = app.Channels;

            if size(mdis, 1) == 1
                [H, W] = size(mdis{1,1});
            else
                [H, W] = size(mdis{1});
            end

            R = zeros(H, W);
            G = zeros(H, W);
            B = zeros(H, W);

            for i = 1:length(chs)
                ch = chs(i);
                if ~ch.enabled || ch.mapIndex < 1 || ch.mapIndex > length(mdis)
                    continue;
                end

                if size(mdis, 1) == 1
                    M = double(mdis{1, ch.mapIndex});
                else
                    M = double(mdis{ch.mapIndex});
                end

                if ~isequal(size(M), [H, W])
                    M = utils.imresize(M, [H, W], 'bicubic');
                end

                if options.perMapNormalize
                    M = utils.mat2gray(M);
                end

                if isfinite(ch.saturation)
                    M(M > ch.saturation) = ch.saturation;
                end

                if ~options.perMapNormalize
                    M = utils.mat2gray(M);
                end

                % Apply contrast limits
                cMin = contrastMins(i);
                cMax = contrastMaxs(i);
                if cMax > cMin && (cMin > 0 || cMax < 1)
                    M = max(0, min(1, (M - cMin) / (cMax - cMin)));
                end

                % Apply gamma
                if gammaVals(i) ~= 1.0 && gammaVals(i) > 0
                    M = M .^ gammaVals(i);
                end

                M = M * ch.brightness;

                R = R + M * ch.rgb(1);
                G = G + M * ch.rgb(2);
                B = B + M * ch.rgb(3);
            end

            % Final normalization
            switch options.normalizationMode
                case 'None'
                    R = max(0, min(1, R));
                    G = max(0, min(1, G));
                    B = max(0, min(1, B));
                case 'Global mat2gray'
                    comp = cat(3, R, G, B);
                    comp = utils.mat2gray(comp);
                    R = comp(:,:,1); G = comp(:,:,2); B = comp(:,:,3);
                case 'Per-channel normalize'
                    R = utils.mat2gray(R);
                    G = utils.mat2gray(G);
                    B = utils.mat2gray(B);
                otherwise
                    R = max(0, min(1, R));
                    G = max(0, min(1, G));
                    B = max(0, min(1, B));
            end

            composite = cat(3, R, G, B);
        end

        function channels = parseChannelsFromJS(app, jsChannels)
            channels = app.Channels;
            channelRGB = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [0 1 1], [1 0 1]};

            for i = 1:min(length(jsChannels), 6)
                ch = jsChannels{i};
                if isstruct(ch)
                    if isfield(ch, 'enabled'),    channels(i).enabled    = logical(ch.enabled); end
                    if isfield(ch, 'mapIndex'),   channels(i).mapIndex   = ch.mapIndex; end
                    if isfield(ch, 'brightness'), channels(i).brightness = ch.brightness; end
                    if isfield(ch, 'saturation')
                        if ch.saturation < 0
                            channels(i).saturation = Inf;
                        else
                            channels(i).saturation = ch.saturation;
                        end
                    end
                    channels(i).lumination = 1.0;
                    if isfield(ch, 'rgb') && ~isempty(ch.rgb)
                        channels(i).rgb = double(ch.rgb);
                    else
                        channels(i).rgb = channelRGB{i};
                    end
                    % Gamma and contrast limits are JS-side only (for preview)
                    % Full-res export uses MATLAB's +core/ pipeline
                end
            end
        end
    end

    % =====================================================================
    %  PUBLIC API
    % =====================================================================

    methods (Access = public)
        function setData(app, mdis, Original)
            app.MapsData = mdis;
            app.DataLoaded = true;
            if nargin >= 3
                app.OriginalBg = Original;
            end
            autoEnableChannels(app);
            sendMapsToJS(app);
        end

        function composite = getComposite(app, options)
            if nargin < 2
                options = struct('normalizationMode', 'Global mat2gray', ...
                                 'perMapNormalize', false);
            end
            composite = core.createCompositeImage(app.MapsData, app.Channels, options);
            if ~isempty(app.OriginalBg)
                composite = core.blendWithBackground(composite, app.OriginalBg, 0.7);
            end
        end
    end
end
