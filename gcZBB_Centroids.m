function [C] = gcZBB_Centroids(datadir)
% Gets ZBB centroids

C = [];
load(fullfile(datadir,'gmranat'));

% Get ZBB labels
nzs = size(gmranat.z,2);
for z = 1:nzs
    C = cat(1,C, ...
        cat(1,gmranat.z(z).STATScrop.Centroid_ZBB));
end
end