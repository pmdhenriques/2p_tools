function [] = gcPlotConvRate(datadir)
% Plots the convergence rate (Stim-ON and Spontaneous for each visual
% stimulus presented
%
% Pedro Henriques, Mar 2018

load(fullfile(datadir,'gm'),'gm');
load(fullfile(datadir,'gmv'),'gmv');
load(fullfile(datadir,'gmb'),'gmb');


nvistypes = size(gmv.vistypz,1);
CR = NaN(nvistypes,2);
for v = 1:nvistypes
    % Stimulus on times
    [~, ~, vSTt, vEDt] = gcStimInfo(gmv.vistypz(v,1:5),gm.trfr,gm.frtime,0);
    
    % Stimulus epochs
    sp = find(ismember(gmv.visstim(:,1:5),gmv.vistypz(v,1:5),'rows'));
    
    % Good convergences
    convix = gmb.convergences(:,3) == 1;
    pconv = intersect(sp,gmb.convergences(convix,1));
    
    % No convergence epochs
    npconv = setdiff(sp,pconv);
    
    % Convergences in stim on time
    convonix = gmb.convergences(:,2) >= vSTt & ...
        gmb.convergences(:,2) <= vEDt;
    pconvon = intersect(sp,gmb.convergences(convix & convonix,1));
    
    % Spontaneous convergences
    pnconvon = intersect(sp,gmb.convergences(convix & ~convonix,1));
    
    crate = length(pconvon)/length(sp);
    sprate = length(pnconvon)/length(sp);
    
    CR(v,1) = crate;
    CR(v,2) = sprate;
end

figure
bar(CR)
legend({'Stim','Spont'})
grid on
xlabel('Visual Stimulus')
ylabel('Convergence rate')
end