function [] = cluster_k_ihb_ph_batch(filepath)

corr_list = [.7 .75 .8]'; % list of correlation thresholds to test
minoccupancy = 8; % min number of cells per cluster
seeded = 0; % not recommended to use seeding
clustermethod = 3;

[filedir,~,~] = fileparts(filepath);

load(filepath,'X');

for c = 1:size(corr_list,1)
    
    thiscorr = corr_list(c,1);
    
    [Idx, Cent, nclust, clsizesorder, corrthr_set] = ....
        ihbcluster_v3(X', thiscorr, minoccupancy, seeded, clustermethod);
    
    corr_list(c,2) = size(Cent,1);
    corr_list(c,3) = nclust;
    
    clust.corr(c).corr = thiscorr;
    clust.corr(c).Idx = Idx;
    clust.corr(c).Cent = Cent;
    clust.corr(c).nclust = nclust;
    clust.corr(c).clsizesorder = clsizesorder;
    clust.corr(c).corrthr_set = corrthr_set;
    
    save(fullfile(filedir,'clust'),'clust','-v7.3');
end
end