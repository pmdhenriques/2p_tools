function [gcbf] = gcbft_verg(datadir)
% Adds vergence information to gmbf structure, based of vergence fitting
% function
%
% Pedro Henriques, April 2018

disp('Loading gmb')
load(fullfile(datadir,'gmb'),'gmb');
disp('Loading gmbf')
load(fullfile(datadir,'gmbf'),'gmbf');

if ~isfield(gmb,'vfit')
    disp('Vergence fit not detected... Running gcVergenceFit')
    gmb = gcVergenceFit(datadir);
end
   
disp('Processing')
for i = 1:size(gmbf.b,2)
    t = gmbf.b(i).ed;
    p = gmbf.b(i).p;
    tix = findnearest(t,gmb.p(p).tt);
    verg = gmb.p(p).Langles(tix)-gmb.p(p).Rangles(tix);
    if verg > gmb.vfit.gintersect
        gmbf.b(i).Vergbout = 1;
    else
        gmbf.b(i).Vergbout = 0;
    end
end

disp('Saving gmbf')
save(fullfile(datadir,'gmbf.mat'),'gmbf')
end