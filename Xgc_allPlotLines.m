function Xgc_allPlotLines(gm,gmv,ylims)
% Plots the stimuli start, end and epoch lines
%
% Pedro Henriques, June 2017

for v = 1:size(gmv.vistypz,1)
    [vST, vED] = gcStimInfo(gmv.vistypz(v,:), gm.trfr, gm.frtime, 0);
    line([vST+(gm.nfr*(v-1)) vST+(gm.nfr*(v-1))],ylims,'Color','b','LineWidth',2);
    line([vED+(gm.nfr*(v-1)) vED+(gm.nfr*(v-1))],ylims,'Color','r','LineWidth',2);
    line([gm.nfr*(v-1) gm.nfr*(v-1)],ylims,'Color','k','LineWidth',3);
end
end