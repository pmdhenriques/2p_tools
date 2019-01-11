function [output,C] = imMaxCorr2(regim,floatim,showplot)

if nargin < 3
    showplot = 0;
end

C = NaN(size(regim,3),1);
for i = 1:size(regim,3)
    C(i) = corr2(floatim,regim(:,:,i));
end
[maxC,maxCix] = max(C);

output = [maxC,maxCix];

if showplot
    figure
    plot(C)
    hold on
    plot(maxCix,maxC,'ro')
    ylabel('Correlation coeff')
    xlabel('Ref Z-position')
    grid on
end

end
