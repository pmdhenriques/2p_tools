function [] = gcROIaddZBBIdxMasks(datadir)
% Adds ZBB masks to each transformed centroids.
% Requires centroids to have been previously registered onto the ZBB
% reference
%
% Pedro Henriques, Aug 2017

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
end

h = waitbar(0,'Loading `gmranat` structure...');
load([datadir filesep 'gmranat.mat']);
waitbar(0,h,'Loading MaskDatabase...');
try
    load('ZBBMaskDatabase');
catch
    [ZBBmaskfile,ZBBmaskdir] = uigetfile('\\128.40.168.141\bdata2\Registration\*.mat', ...
        'Select ZBBMaskDatabase.mat file');
    load(fullfile(ZBBmaskdir,ZBBmaskfile));
end
    

zbbdim = [616 1030 420];
nzs = size(gmranat.z,2);
nmsks = size(ZBBMaskDatabase,2);

for z = 1:nzs
    nzcnt = size(gmranat.z(z).STATScrop,1);
    zbbxyz = round(cat(1,gmranat.z(z).STATScrop.Centroid_ZBB));
    
    % Extract mask ids for each centroid
    M = false(nzcnt,nmsks);    
    for m = 1:nmsks
        waitbar(m/nmsks,h, ...
            sprintf('Attributing masks Z = %d/%d',z,nzs))
        Mask = reshape(full(ZBBMaskDatabase(:,m)),zbbdim);
        M(:,m) = Mask(sub2ind(zbbdim,zbbxyz(:,2),zbbxyz(:,1),zbbxyz(:,3)));
    end    
    
    % Write to structure
    for i = 1:nzcnt
        gmranat.z(z).STATScrop(i).Masks_ZBB = M(i,:);
    end
end

clear ZBBMaskDatabase

waitbar(1,h,'Saving...')
save(fullfile(datadir,'gmranat.mat'),'gmranat');
close(h);
end