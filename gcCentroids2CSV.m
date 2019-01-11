function [] = gcCentroids2CSV(datadir,rot,zcorr)
% Writes csv file of original roi centroids in a format that can be used by
% ANTs to transform the points to a new reference frame
% (antsApplyTransformsToPoint)
%
% Pedro Henriques, Aug 2017

if nargin < 2
    rot = 0;
elseif nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
end

disp('Loading...')
load(fullfile(datadir,'gmranat.mat'))
load(fullfile(datadir,'gm.mat'))

nzs = size(gmranat.z,2);

if nargin < 3 || isempty(zcorr)
    zcorr = 1:nzs;
end

xpx = gm.xpx;
ypx = gm.ypx;
if nzs > 2
    zpx = diff(gm.zrange(1:2));
else
    zpx = 1;
end

XYZ = [];
for z = 1:nzs    
    XYZ = cat(1,XYZ,[cat(1,gmranat.z(z).STATScrop.Centroid), ...
        ones(size(gmranat.z(z).STATScrop,1),1)*zcorr(z)]);
end

% perform a 180 deg counter-clockwise rotation to all point to be in sync
% with image rotation performed before CMTK registration

if rot
    disp('Performing 180 deg rotation!')
    t = [ ...
        -1 0 0 0; ...
        0 -1 0 0; ...
        0 0 1 0; ...
        0 0 0 1];
    tform = affine3d(t);
    XYZ = transformPointsForward(invert(tform),XYZ);
    XYZ = XYZ+[gm.xDef+1 gm.yDef+1 0];
end
XYZ = XYZ.*[xpx ypx zpx];
XYZ(:,4:6) = [ones(size(XYZ,1),2) NaN(size(XYZ,1),1)];

% Build a table with proper labels
T = array2table(XYZ,'VariableNames',{'x','y','z','t','label','comment'});

if ~exist(fullfile(datadir,'registration'),'dir')
    mkdir(fullfile(datadir,'registration'));
    disp('Creating registraiton folder')
end

% Save

savedir = fullfile(datadir,'registration','centroids2ref.csv');
if exist(savedir,'file')
    tdstr = datestr(datetime('Today'),'yymmdd');
    writetable(T,strrep(savedir,'.csv',['_' tdstr '.csv']));
else
    writetable(T,savedir);
end
disp('File saved!')

end