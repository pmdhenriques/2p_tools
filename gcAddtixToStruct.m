function [] = gcAddtixToStruct(datadir)
% Adds total indexes and label to all ROIs in the gmranat structure (helps
% with indexing in some cases...)
%
% Pedro Henriques, June 2017

load(fullfile(datadir,'gmranat'));

nzs = size(gmranat.z,2);
    
roiid = [];
for z = 1:nzs
    roiid = cat(1, roiid, ...
        cat(2, ...
        setdiff(gmranat.z(z).allIdx,gmranat.z(z).excludeIdx), ...
        ones(size(gmranat.z(z).STATScrop,1),1).*z));
end

for z = 1:size(gmranat.z,2)
    gmranat.z(z).STATScrop(1).Label = [];
    labels = setdiff(gmranat.z(z).allIdx,gmranat.z(z).excludeIdx);
    for i = 1:size(gmranat.z(z).STATScrop,1)
        gmranat.z(z).STATScrop(i).Label = labels(i);
        gmranat.z(z).STATScrop(i).tix = find(roiid(:,2) == z & ...
            roiid(:,1) == labels(i));
    end
end

save(fullfile(datadir,'gmranat'),'gmranat');
end