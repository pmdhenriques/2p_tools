function [] = gcIm2NifTI(datadir,rot)
% Saves averaged stacks for the imaged planes onto NifTI 32bit format which
% can be used for registration
%
% Pedro Henriques, June 2017

if nargin < 2
    rot = 0;
elseif nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P');
end

disp(['This is gcIm2NifTI, processing ' datadir '...'])
h = waitbar(0.5,'Loading `gm` structure...');
load(fullfile(datadir,'gm'),'gm')

nz = size(gm.zindices,2);
dim = [gm.yDef,gm.xDef, nz];
if nz > 1
    pxdim = [gm.xpx,gm.ypx,diff(gm.zrange(1:2))];
else
    pxdim = [gm.xpx,gm.ypx,1];
end

stk = zeros(dim);
for z = 1:dim(3)
    if isfield(gm,'aligntwice')
        stk(:,:,z) = gm.zimg(z).aligngood2;
    else
        stk(:,:,z) = gm.zimg(z).aligngood;
    end
end

% The stack is rotated to obey to the standard fish orientation (Anterior
% to left)

if rot
    stk = flip(rot90(stk),2);
else
    stk = flip(rot90(stk,3),2);
end

nii = make_nii(stk, ...
    pxdim, ...
    [0 0 0], ...
    512);

waitbar(0.5,h,'Saving nii...');
if ~exist(fullfile(datadir,'registration'),'dir')
    mkdir(fullfile(datadir,'registration'));
end

% Save
save_nii(nii,fullfile(datadir,'registration',[gm.name '_01.nii']));

close(h)
end