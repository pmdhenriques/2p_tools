function [roitix] = gcROIlab2tix(roilab,z,gmVCNV)
% Transform label roi indexes to total roi indexes for gm structures
%
% Pedro Henriques, June 2017

zix = find(gmVCNV.ROIinfo(:,2) == z);
lia = ismember(gmVCNV.ROIinfo(zix,1),roilab,'rows');
roitix = zix(lia);

end
