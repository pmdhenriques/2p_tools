function [zcorr] = gcCentroidsRegZCorr(datadir,stackdir)
% Computes the Z correction needed to apply to centroid positions for
% reistration purposes. Used if a larger brain stack has been performed and
% used for registration, instead of the average aligned stack.
% 
% Assumes that the experiment was performed on z positions included in the
% larger stack
%
% stackdir - directory with the single tiff files obtained in the larger
% stack
%
% Pedro Henriques, Fev 2018

load(fullfile(datadir,'gm'),'gm')
files = dir2(fullfile(stackdir,'*.tif'));

if isempty(files)
    error('No TIFF files fount on stackdir')
    return
end

nfiles = size(files,1);
zrange = gm.zrange;

zpos = zeros(nfiles,1);
for i = 1:nfiles
    namespp = strsplit(files(i).name,'.');
    zpos(i) = str2num(namespp{4});
end
zpos = sort(zpos,'ascend');

zcorr = find(ismember(zpos,zrange,'rows'));
end
