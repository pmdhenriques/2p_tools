function [] = gcROIplotcentroids(rois,zs,gm,gmranat)
% Plots the centroid coordenates of input ROIs onto the average image for
% their given z slice
%
% Pedro Henriques, June 2017

nroitype = size(rois,2);

if isempty(zs)
    zs = 1:size(gm.zindices,2);
end

for i = 1:length(zs)
    z = zs(i);
    
    C = cat(1,gmranat.z(z).STATScrop.Centroid);
    tix = cat(1,gmranat.z(z).STATScrop.tix);
    
    figure('name',sprintf('ROIs for Z = %d',z))
    colormap gray
    imagesc(gm.zimg(z).aligngood2,[0 20])
    hold on    
    
    for j = 1:nroitype
        ix = ismember(tix,rois(:,j));
        scatter(C(ix,1),C(ix,2),'filled')
    end
    axis('square','off')
end
end
