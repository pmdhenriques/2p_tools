function [] = gcROIplotconv(rois,visstimix,gm,gmv,gmb,gmrxanat)
% Plots ROI responses of convergence VS non-convergence epochs
%
% Pedro Henriques, June 2017

if nargin < 3
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
    disp('Loadinf structures')
    load(fullfile(datadir,'gm'),'gm');
    load(fullfile(datadir,'gmv'),'gmv');
    load(fullfile(datadir,'gmb'),'gmb');
    load(fullfile(datadir,'gmrxanat'),'gmrxanat');
    if nargin < 2
        visstimix = 1:size(gmv.vistypz,1);
    end
end

if size(rois,1) > size(rois,2)
    rois = rois';
end

for v = visstimix
    for roi = rois
        z = gmrxanat.roi(roi).z;
        e = gm.zindices(z).e;
        
        vepoch = find(ismember(gmv.visstim(:,1:5),gmv.vistypz(v,1:5),'rows'));  % stimulus epochs
        zvepoch = e(ismember(e,vepoch));    % epochs with current z
        
        [vST, vED, vSTt, vEDt] = gcStimInfo(gmv.vistypz(v,1:5), gm.trfr, gm.frtime, 0);        
        
        czvepochix = ismember(gmb.convergences(:,1),zvepoch) & ...
            gmb.convergences(:,2) > vSTt & ...
            gmb.convergences(:,2) < vEDt & ...
            gmb.convergences(:,3) == 1; % Epochs in which convergence is coincidental with stimulus on time
        czvepoch = gmb.convergences(czvepochix,1);
        
        if isempty(czvepoch)
            fprintf('No convergence for ROI %d\n',roi);
            continue
        end
        
        cits = round(gmb.convergences(czvepochix,2)./gm.frtime.*1000);  % Convergence iteration (in 2P time)        
        czvepochid = find(ismember(zvepoch,czvepoch));  % Convergence epochs zs
        nzvepochid = setdiff(1:length(zvepoch),czvepochid); % Non convergence ephochs
        
        X = gmrxanat.roi(roi).Vprofiles(v).zProfiles';
        Xz = (X-nanmean(reshape(X,numel(X),1)))/ ...
            nanstd(reshape(X,numel(X),1));
        Xzf = filtfilt(ones(1,3)/3,1,double(Xz));
        
        lims = [-2 5];  % Ylim in std
%         lims(1) = min(min(Xf))-0.2;
%         lims(2) = max(max(Xf))+0.2;
        
        figure('name',sprintf('ROI # %d; V = %d',roi,gmv.vistypz(v,1)), ...
            'Position',[400 400 1000 400])
        subplot(1,2,1)
        plot(Xzf(:,czvepochid),'k');
        hold on
        for i = 1:length(cits)
            p1 = plot(cits(i),Xzf(cits(i),czvepochid(i)),'b*');
        end
        xlim([0 gm.nfr]); ylim([lims(1) lims(2)]);
        ylabel('F (zscore)'); xlabel('Frames');
        title('Convergence')
        p2 = line([vST vED; vST vED],[lims(1) lims(1); lims(2) lims(2)], ...
            'LineStyle','--');
        legend([p1(1),p2(1),p2(2)],{'Conv','vST','vED'})
        
        subplot(1,2,2)
        plot(Xzf(:,nzvepochid),'k');
        xlim([0 gm.nfr]); ylim([lims(1) lims(2)]);
        xlabel('Frames');
        title('Non-Convergence')
        line([vST vED; vST vED],[lims(1) lims(1); lims(2) lims(2)], ...
            'LineStyle','--')
    end
end

