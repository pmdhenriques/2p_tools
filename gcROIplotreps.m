function [] = gcROIplotreps(rois,visstimix,gmrxanat,gm,gmv,gmb,gmbt)
% Plots, for every ROI and visstim index parsed, a figure containing the
% zscored activity along eye and tail behaviour traces for each rep
%
% Pedro Henriques, June 2017

if nargin < 3
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
    disp('Loadinf structures')
    load(fullfile(datadir,'gm'),'gm');
    load(fullfile(datadir,'gmv'),'gmv');
    load(fullfile(datadir,'gmb'),'gmb');
    load(fullfile(datadir,'gmbt'),'gmbt');
    load(fullfile(datadir,'gmrxanat'),'gmrxanat');
    if nargin < 2
        visstimix = 1:size(gmv.vistypz,1);
    end
end

for roi = rois'
    for v = visstimix
        % Get the F profiles for all reps of the ROI
        X = gmrxanat.roi(roi).Vprofiles(v).zProfiles';
        %
        if gmrxanat.nfr == 4*gm.nfr
            nreps = size(X,2);
            X = movsum(X(:),[0 3]);
            X = X(1:4:end);
            X = reshape(X,gm.nfr,nreps);
        end
        % Z score
        X = (X-nanmean(reshape(X,numel(X),1)))/ ...
            nanstd(reshape(X,numel(X),1));
        % Light zero-phase filter
        Xf = filtfilt(ones(1,3)/3,1,double(X));
        
        reps = size(X,2);
        nfrm = size(X,1);
        
        mlims = [-2 5];     % plot limits of the zscore
        blims = [-40 40];   % plot limits of the eyes    
        tlims = [0 70];     % plot limits of the tail
        
        % Get start and end frames for current visstim
        [vST, vED] = gcStimInfo(gmv.vistypz(v,:), gm.trfr, gm.frtime, 0);
        
        z = gmrxanat.roi(roi).z;    % z pos
        % Indices of the vistim reps
        e = [gm.zindices(z).e(1) gm.zindices(z).e(end)];    
        Lia = find(ismember(gmv.visstim(:,1:5),gmv.vistypz(v,:),'rows'));
        Lia = Lia(Lia >= e(1) & Lia <= e(2));
        
        % Tail vigour resampled tp 2P acquisition
        t = cat(1,gmbt.p(Lia).resampled);
        t = reshape(t(:,3),size(X));
        
        % Eye angles resampled tp 2P acquisition
        b = cat(1,gmb.p(Lia).resampled);
        b = reshape(b(:,2:3)',[2 size(X)]);
        
        %% Plot Stuff
        figure('name',sprintf('ROI # %d; V = %d',roi,gmv.vistypz(v,1)), ...
            'Position',[200 200 1500 700])
        for r = 1:reps
            % zscore
            subplot(4,reps,[r r+reps])
            plot(Xf(:,r),'k'); title(sprintf('rep %d',r));
            xlim([0 nfrm]); ylim(mlims);
            line([vST vED; vST vED],[mlims(1) mlims(1); mlims(2) mlims(2)], ...
                'Color','k','LineStyle','--')
            % eyes
            subplot(4,reps,r+reps*2)
            plot(b(1,:,r)); hold on
            plot(b(2,:,r)); hold off
            xlim([0 nfrm]); ylim(blims);
            line([vST vED; vST vED],[blims(1) blims(1); blims(2) blims(2)], ...
                'Color','k','LineStyle','--')      
            % tail
            subplot(4,reps,r+reps*3)
            plot(t(:,r));
            xlim([0 nfrm]); ylim(tlims);
            line([vST vED; vST vED],[tlims(1) tlims(1); tlims(2) tlims(2)], ...
                'Color','k','LineStyle','--')
        end
    end
end
end