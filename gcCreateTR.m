function [] = gcCreateTR(datadir)
% Create matrix of entire response profiles for each individual ROI
%
% Pedro Henriques, June 2017

h = waitbar(0.5,'Loading structures...');
load(fullfile(datadir,'gm'));
load(fullfile(datadir,'gmv'));
load(fullfile(datadir,'gmranat'));
load(fullfile(datadir,'gmrxanat'));

nstim = size(gmrxanat.roi(1).Vprofiles,2);
nreps = size(gmrxanat.roi(1).Vprofiles(1).zProfiles,1);
nzs = size(gm.zindices,2);

visp = zeros(nstim,nreps,nzs);
nzstim = nstim*nreps;

for i = 1:size(gmv.visstim,1)
    [~,r] = ismember(gmv.visstim(i,1:5),gmv.vistypz,'rows');
    visp(r,nnz(visp(r,:,ceil(i/nzstim)))+1,ceil(i/nzstim)) = i-(nzstim*(ceil(i/nzstim)-1));
end

TR = NaN(size(gmrxanat.roi,2),nstim*nreps*gm.nfr);
for z = 1:nzs
    waitbar(z/nzs,h,sprintf('Building TR for Z = %d',z))
    for i = [gmranat.z(z).STATScrop.tix]
        R = [];
        for j = 1:nzstim
            [r,c] = find(visp(:,:,z) == j);
            if isempty(c)
                R = [R zeros(1,gm.nfr)];
            else
                R = [R gmrxanat.roi(i).Vprofiles(r).zProfiles(c,:)];
            end
        end
        TR(i,:) = R;
    end
end

waitbar(1,h,'Saving...')
save(fullfile(datadir,'TR'),'TR');
close(h)
end