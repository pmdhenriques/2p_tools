function [maskrois,roizs] = gcROIextractZBBIdxMask(datadir,maskstr,maskreference,ZBBMaskDatabaseNames)
% Extracts the index of ROIs that belong to the input ZBB mask
%
% Pedro Henriques, Aug 2017

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
end

disp('Loading `gmranat` structure...')
load(fullfile(datadir,'gmranat.mat'),'gmranat');

if nargin < 4
    disp('Loading MaskDatabase...')
    try
        load('ZBBMaskDatabase','ZBBMaskDatabaseNames');
    catch
        disp('ZBBMaskDatabase not found on path')
        [maskfile, maskdir] = uigetfile('C:\Users\Pedro\Documents\MATLAB', ...
            'Please select ZBBMaskDatabase.mat file');
        disp('Loading ZBBMaskDatabase')
        load(fullfile(maskdir,maskfile),'ZBBMaskDatabaseNames');
    end
end
    
nzs = size(gmranat.z,2);

% Which mask you're looking for
if nargin < 2
    usrinput = 0;
    maskchoice = questdlg('How do you want to select a mask?', ...
        'Mask Selection', ...
        'List','User input','List');
    switch maskchoice
        case 'List'
            ok = 0;
            while ~ok
                [MaskIdx,ok] = listdlg('PromptString','Select a mask:',...
                    'SelectionMode','single',...
                    'ListString',ZBBMaskDatabaseNames(:));
            end
        case 'User input'
            usrinput = 1;
            maskstr = inputdlg('What''s the name of the mask?');
            maskstr = maskstr{1};
    end
else
    usrinput = 1;
end

% Which Masks index to use
if nargin < 3
    statsfields = fieldnames(gmranat.z(1).STATScrop);
    statsfields = statsfields(contains(statsfields,'Mask'));
    while isempty(statsfields)
        addmasks = questdlg('Mask indexes not found... Run gcROIaddZBBIdxMasks?', ...
            'Run gcROIaddZBBIdxMasks?', ...
            'Yes','No','Yes');
        switch addmasks
            case 'Yes'
                gcROIaddZBBIdxMasks(datadir);
                disp('Loading `gmranat` structure again...')
                load([datadir filesep 'gmranat.mat']);
                statsfields = fieldnames(gmranat.z(1).STATScrop);
                statsfields = statsfields(contains(statsfields,'Mask'));
            otherwise
                return
        end
    end
    
    ok = 0;
    while ~ok
        [maskreferenceix,ok] = listdlg('PromptString','Select the mask reference to use:',...
            'SelectionMode','single',...
            'ListString',statsfields);
    end
    maskreference = statsfields{maskreferenceix};
end

%%

% Get mask idx
if usrinput
    MaskIdx = find(contains(ZBBMaskDatabaseNames,maskstr));
    if ~isempty(MaskIdx)
        if length(MaskIdx) > 1
            ok = 0;
            disp('More than one mask found...')
            while ~ok
                [maskselect,ok] = listdlg('PromptString','Select a mask:',...
                    'SelectionMode','single',...
                    'ListString',ZBBMaskDatabaseNames(MaskIdx));
            end
            MaskIdx = MaskIdx(maskselect);
        end
    else
        errordlg(sprintf('No masks found with name %s',maskstr));
    end
end

% Get roi idxs

while ~isfield(gmranat.z(1).STATScrop,'tix')
    addtix = questdlg('tix indexes not found... Run gcAddtixToStruct?', ...
        'Run gcAddtixToStruct?', ...
        'Yes','No','Yes');
    switch addtix
        case 'Yes'
            gcAddtixToStruct(datadir);
            disp('Loading `gmranat` structure again...')
            load([datadir filesep 'gmranat.mat']);
        otherwise
            return
    end
end

maskrois = [];
roizs = [];
for z = 1:nzs
    zmasksix = cat(1,gmranat.z(z).STATScrop.(maskreference));
    ix = zmasksix(:,MaskIdx) == 1;
    maskrois = [maskrois; cat(1,gmranat.z(z).STATScrop(ix).tix)];
    roizs = [roizs; ones(sum(ix),1).*z];
end

fprintf('Found %d ROIs\n',length(maskrois));

end