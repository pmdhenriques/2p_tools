function [TPr,FNr] = gcConv_error(datadir)
% Computes the True positive and False negative rate of the online
% convergence detector VS offline convergence detector

load(fullfile(datadir,'gmb'),'gmb');

tconv = gmb.convergences(gmb.convergences(:,3) == 1,1);
C = dlmread(fullfile(datadir,'b\convlog.xls'));
C(:,2) = C(:,2)+1;

TPr = sum(ismember(tconv,C(:,2)))/length(tconv);
FNr = sum(~ismember(tconv,C(:,2)))/length(tconv);

end

