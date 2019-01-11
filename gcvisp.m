function [visp] = gcvisp(gmv,gmrxanat)
% Returns the trial number array organized by stimulus type, rep, z

nstim = size(gmrxanat.roi(1).Vprofiles,2);
nreps = size(gmrxanat.roi(1).Vprofiles(1).zProfiles,1);
nzs = length(unique([gmrxanat.roi.z]));

visp = zeros(nstim,nreps,nzs);
nzstim = nstim*nreps;
for i = 1:size(gmv.visstim,1)
    z = ceil(i/nzstim);
    [~,r] = ismember(gmv.visstim(i,1:5),gmv.vistypz,'rows');
    visp(r,nnz(visp(r,:,z))+1,z) = i;
    % = i-(nzstim*(z-1)); % counting from 1 each z
end
end