function [ResponseRate,ResponseNum] = gcResponseRate(datadir,visstim)
% Calculated convergence response rates for input visstim
% 
% Pedro Henriques, Feb 2018

disp('Loading structures')
load(fullfile(datadir,'gmv.mat'),'gmv');
load(fullfile(datadir,'gm.mat'),'gm');
load(fullfile(datadir,'gmb.mat'),'gmb');

nstim = size(visstim,1);
ResponseRate = NaN(nstim,1);
ResponseNum = NaN(nstim,2);

for s = 1:nstim
    pcstim = find(ismember(gmv.visstim(:,1:4),visstim(s,1:4),'rows'));
    [~, ~, vSTt, vEDt] = gcStimInfo(visstim(s,1:5), gm.trfr, gm.frtime, 0);
    pcconv = find(ismember(gmb.convergences(:,1),pcstim) & ...
        gmb.convergences(:,2) >= vSTt &  ...
        gmb.convergences(:,2) <= vEDt);
    ResponseRate(s) = length(pcconv)/length(pcstim);
    ResponseNum(s,1) = length(pcconv);
    ResponseNum(s,2) = length(pcstim);
end

end