function [] = gcCentroids2NifTI(csvdir,dim)
% Writes a binanty stack with points in the input centroid positions
%
% Pedro Henriques, Fev 2018

if nargin < 2
    dim = [1030 616 420];   % ZBB reference
end

[filepath,~,~] = getfilepath(csvdir);

T = readtable(csvdir);
X = zeros(dim,'uint8');

x = round(T.x);
y = round(T.y);
z = round(T.z);

inds = sub2ind(dim,x,y,z);
X(inds) = 255;

nii = make_nii(X);
save_nii(nii,fullfile(filepath,'Centroids.nii.gz'));
end