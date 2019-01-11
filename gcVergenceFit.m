function [gmb] = gcVergenceFit(datadir,bounds)
% Finds vergence thresholds for convergences by fitting a mixture of two
% gaussians to the eye vergence angle (only on convergence trials
%
% Pedro Henriques, April 2018

if nargin < 2
    bounds = [15,52,52,70]; % Gaussian fit bounds
end

histrng = -40:2.5:120;

disp('Loadig gmb')
load(fullfile(datadir,'gmb'),'gmb');

convp = gmb.convergences(gmb.convergences(:,3) == 1,1);
if length(convp) < 10
    fprintf(2,'Not enought convergences to proceed with fit! \n')
    return
end

% nconvp = setdiff(1:size(gmb.p,2),convp);
% if length(nconvp) < length(convp)
%     nsample = length(nconvp);
% else
%     nsample = length(convp);
% end
% nconvp = datasample(setdiff(1:size(gmb.p,2),convp),nsample, ...
%     'Replace',false)';
% C = [convp; nconvp];

C = convp;
V = cat(1,gmb.p(C).Langles)-cat(1,gmb.p(C).Rangles);
[counts, centres] = hist(V, histrng);

%% Fit

disp('Fitting Gaussian mixture')

aFitType = fittype('gauss2');
aFitOptions = fitoptions('gauss2');
aFitOptions.Lower(1:6) = -Inf;
aFitOptions.Upper(1:6) = Inf;
aFitOptions.Lower(2) = bounds(1);
aFitOptions.Lower(5) = bounds(2);
aFitOptions.Upper(2) = bounds(3);
aFitOptions.Upper(5) = bounds(4);
aFitOptions.Robust = 'LAR';
gfit = fit(centres', counts', aFitType, aFitOptions);

coeffvals = coeffvalues(gfit);
a1 = coeffvals(1);
a2 = coeffvals(4);

b1 = coeffvals(2);
b2 = coeffvals(5);

c1 = coeffvals(3)./sqrt(2);
c2 = coeffvals(6)./sqrt(2);

A = 1/(2*c1.^2) - 1/(2*c2.^2);
B = b2/c2.^2 - b1/c1.^2;
C = b1.^2/(2*c1.^2) - b2.^2/(2*c2.^2) - log((c2*a1)/(c1*a2));

grootz = roots([A; B; C]);
gintersect = findnearest(bounds(2), grootz);
gintersect = grootz(gintersect);
pthrON = b2 - 1.*c2; % ie vergnce centre minus 1 sd (was gintersect)
pthrOFF = gintersect - 5; % ie 5 degrees less than intersect

%%  Write to structure

gmb.vfit.gfit = gfit;
gmb.vfit.pthrON = pthrON;
gmb.vfit.pthrOFF = pthrOFF;
gmb.vfit.gintersect = gintersect;

%%  Save

disp('Saving')
save(fullfile(datadir,'gmb'),'gmb');

%%  Plot

figure('Name','Vergence Histogram')
hold on
hist(V, histrng)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'w';
xlabel('Vergence [degrees]')
plot(gfit, centres, counts)
line([pthrON pthrON],[0 .25.*a2], 'Color', 'r', 'LineWidth', 3)
line([pthrOFF pthrOFF],[0 .25.*a2], 'Color', 'b', 'LineWidth', 3)
text(gintersect, 1.1.*a2, num2str(gintersect,'%.1f'))
text(b1, 1.1.*a1, num2str(b1,'%.1f'))
text(b2, 1.1.*a2, num2str(b2,'%.1f'))
end