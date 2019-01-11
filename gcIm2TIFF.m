function [] = gcIm2TIFF(datadir,rot)
% Saves averaged stacks for the imaged planes onto TIFF 16bit format which
% can be used for registration
%
% Pedro Henriques, June 2017

if nargin < 2
    rot = 0;
elseif nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P');
end

disp(['This is gcIm2TIFF, processing ' datadir '...'])
h = waitbar(0.5,'Loading `gm` structure...');
load(fullfile(datadir,'gm.mat'))

dim=[gm.yDef,gm.xDef, size(gm.zindices,2)];

stk = zeros(dim);
for z = 1:dim(3)
    if isfield(gm,'aligntwice')
        stk(:,:,z) = gm.zimg(z).aligngood2;
    else
        stk(:,:,z) = gm.zimg(z).aligngood;
    end
end
stk = uint16(stk);

% The stack is rotated to obey to the standard fish orientation (Anterior
% to left)

if rot
    stk = rot90(stk,2);
end

waitbar(0.5,h,'Saving tiff...');
if ~exist(fullfile(datadir,'registration'),'dir')
    mkdir(fullfile(datadir,'registration'));
end

% Save
% save_nii(nii,fullfile(datadir,'registration',[gm.name '_01.nii.gz']));

for z = 1:size(stk,3)
    if z == 1
        imwrite(stk(:,:,z),fullfile(datadir,'registration',[gm.name '_01.tiff']));
    else
        imwrite(stk(:,:,z),fullfile(datadir,'registration',[gm.name '_01.tiff']), ...
            'WriteMode','append');
    end
end

close(h)
end