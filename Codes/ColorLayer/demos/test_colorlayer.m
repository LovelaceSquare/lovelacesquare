function test_colorlayer()
% test_colorlayer. Automated test suite for ColorLayer toolbox.
%
% Runs pass/fail assertions on all core functions and utilities.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

close all;
thisFile = mfilename('fullpath');
colorLayerPath = fileparts(fileparts(thisFile));
addpath(colorLayerPath);

nPass = 0; nFail = 0; failures = {};

%% mat2gray
try
    out = utils.mat2gray([1 2 3; 4 5 6]);
    assert(out(1,1) == 0 && out(2,3) == 1);
    nPass = nPass + 1; fprintf('PASS: mat2gray\n');
catch ME, nFail = nFail + 1; failures{end+1} = ME.message; fprintf('FAIL: mat2gray\n'); end

%% validateMaps
try
    [ok, ~] = utils.validateMaps({rand(10,10), rand(10,10)});
    assert(ok);
    nPass = nPass + 1; fprintf('PASS: validateMaps\n');
catch ME, nFail = nFail + 1; failures{end+1} = ME.message; fprintf('FAIL: validateMaps\n'); end

%% imresize
try
    out = utils.imresize(rand(100,100), 0.5, 'bilinear');
    assert(size(out,1) == 50 && size(out,2) == 50);
    nPass = nPass + 1; fprintf('PASS: imresize\n');
catch ME, nFail = nFail + 1; failures{end+1} = ME.message; fprintf('FAIL: imresize\n'); end

%% createCompositeImage
try
    mdis = {rand(50,50), rand(50,50), rand(50,50)};
    ch(1) = struct('name','Red','enabled',true,'mapIndex',1,'brightness',1,'saturation',Inf,'lumination',1,'rgb',[1 0 0]);
    ch(2) = struct('name','Green','enabled',true,'mapIndex',2,'brightness',1,'saturation',Inf,'lumination',1,'rgb',[0 1 0]);
    ch(3) = struct('name','Blue','enabled',true,'mapIndex',3,'brightness',1,'saturation',Inf,'lumination',1,'rgb',[0 0 1]);
    opts.normalizationMode = 'Global mat2gray'; opts.perMapNormalize = false;
    comp = core.createCompositeImage(mdis, ch, opts);
    assert(size(comp,3) == 3 && max(comp(:)) <= 1);
    nPass = nPass + 1; fprintf('PASS: createCompositeImage\n');
catch ME, nFail = nFail + 1; failures{end+1} = ME.message; fprintf('FAIL: createCompositeImage\n'); end

%% blendWithBackground
try
    blended = core.blendWithBackground(rand(50,50,3), rand(50,50), 0.7);
    assert(isequal(size(blended), [50 50 3]));
    nPass = nPass + 1; fprintf('PASS: blendWithBackground\n');
catch ME, nFail = nFail + 1; failures{end+1} = ME.message; fprintf('FAIL: blendWithBackground\n'); end

%% Summary
fprintf('\n%d passed, %d failed\n', nPass, nFail);
if nFail > 0
    for i = 1:length(failures), fprintf('  - %s\n', failures{i}); end
end
end
