function [] = gcConvMod_at(datadir)
% Computes the convergence modulation for each ROI at AT flightback points
%
% Pedro Henriques, Oct 2018

%% Load structures

h = waitbar(0,'Loading structures');

load(fullfile(datadir,'gm'),'gm');
load(fullfile(datadir,'gmv'),'gmv');
load(fullfile(datadir,'gmb'),'gmb');
load(fullfile(datadir,'gmrxanat'),'gmrxanat');

%%  Setup variables

ffiltwidth = 5; % Fluorescence filter width
ffilter = ones(1,ffiltwidth)./ffiltwidth;   % Fluorescent filter
conv_offset_s = 2;    % Convergence offset for convergece window (s)
conv_offset_fr = ceil(conv_offset_s/(gm.frtime/1000));  % Convergence offset in frames
nshuff = 1000;
nfr = gm.nfr;   % Number of epoch frames
isi = 0.01;

%%

% Visual stimuli types
vtypes = gmv.vistypz(:,1) >= 80006 & gmv.vistypz(:,1) <= 81007;
vstims = gmv.vistypz(vtypes,1:4);   % Visual stimuli codes

nrois = size(gmrxanat.roi,2);   % Number of rois
zs = cat(1,gmrxanat.roi.z);

% Get actual stimulus times
vSTt = (gm.trfr+1)*gm.frtime./1000;
vttime = gmv.visstim(:,6);
vEDst = vSTt + vttime;

%%  Run the loop

for roi = 1:nrois
    waitbar(roi/nrois,h,'Processing')
    
    z = zs(roi);    % z pos
    
    % stimulus epochs to consider
    se = intersect(gm.zindices(z).e, ...
        find(ismember(gmv.visstim(:,1:4),vstims,'rows')));
    
    %% Get convergence epochs
    
    % Good convergences
    gconvix = gmb.convergences(:,3) == 1 & ismember(gmb.convergences(:,1),se);
    % Convergences in stim ON time
    convonix = gmb.convergences(:,2) >= vSTt & ...
        gmb.convergences(:,2) <= vEDst(gmb.convergences(:,1));
    % Stim ON convergence epochs
    econv_on = unique(intersect(se,gmb.convergences(gconvix & convonix,1)));
    % No convergence stim ON epochs
    noeconv = setdiff(se,econv_on);    
    
    % Epoch info
    [v_on_n,~,r_on_n] = gcPresentationInfo(noeconv,gm,gmv);    
    % AT epochs
    v_at_n_at = find(ismember(v_on_n,[5 6]));
    nv_at_n_at = length(v_at_n_at);
    
    %%
    
    % ROI Ca2+ vector of entire experiment
    [v0,~,r0] = gcPresentationInfo(gm.zindices(z).e,gm,gmv);
    F0 = [];
    for i = 1:length(v0)
        F0 = cat(1,F0,gmrxanat.roi(roi).Vprofiles(v0(i)).zProfiles(r0(i),:)');
    end
    F0f = filtfilt(ffilter,1,double(F0));   % zero-phase filter
    F0fz = (F0f-nanmedian(F0f))./nanstd(F0f); % z-score
    
    % Null distribution
    Null = NaN(nshuff,1);
    nullstfr = randi(length(F0)-(2*conv_offset_fr),nshuff,1);
    for i = 1:nshuff
        stfr = nullstfr(i);
        Null(i) = nanmean(F0fz(stfr:stfr+(2*conv_offset_fr)-1));  % Conv centered z F
    end
    nullmean = nanmean(Null);
        
    Fz_at = NaN(nv_at_n_at,nfr,'single');
    F_at = NaN(nv_at_n_at,1); % Mean fluorescence list
    Z_at = NaN(nv_at_n_at,1); % Z list (distance to no response trials distribution)
    Z_null_at = NaN(nv_at_n_at,1);   % Null distribution means
    F_at_nrsp_all = cell(nv_at_n_at,1);
    Fz_at_nrsp = cell(nv_at_n_at,1);
    eyeDelta_at = NaN(nv_at_n_at,1);
    
    for i = 1:nv_at_n_at
        ixx = v_at_n_at(i);
        vis = v_on_n(ixx);  % visual stimulus
        rep = r_on_n(ixx);  % rep
        e = noeconv(ixx);    % epoch
        
        ze = find(v0 == vis & r0 == rep);   % z-wise epoch
        
        cix = find(abs(diff(gmb.p(e).vis_traj(:,3))) > 100);
        cix = cix(1); % use only first stim flighback
        
        cst = gmb.p(e).vis_traj(cix,7);  % Convergence start (s)
        cstfr = round(cst/(gm.frtime/1000));    % Convergence start (fr)        
        cst_z = nfr*(ze-1)+cstfr; % Convergence frame (from start of z plane)
        
        % Eyes
        
        L = gmb.p(e).Langles;
        R = gmb.p(e).Rangles;
        tt = gmb.p(e).tt;
        [tt, ind] = unique(tt);
        timebase = 0:isi:tt(end);
        Ri = interp1(tt, R(ind), timebase);
        Li = interp1(tt, L(ind), timebase);
        eyecix = findnearest(cst,timebase);

        preL = nanmean(Li(eyecix-gmb.tailsearch:eyecix-gmb.tailsearch+gmb.av_int));
        postL = nanmean(Li(eyecix+gmb.tailsearch-gmb.av_int:eyecix+gmb.tailsearch));
        preR = nanmean(Ri(eyecix-gmb.tailsearch:eyecix-gmb.tailsearch+gmb.av_int));
        postR = nanmean(Ri(eyecix+gmb.tailsearch-gmb.av_int:eyecix+gmb.tailsearch));
        
        eyeDelta_at(i,1) = postL-preL; % Left eye delta angle
        eyeDelta_at(i,2) = postR-preR;   % Right eye delta angle
        
        % 
        
        itst = cst_z-conv_offset_fr+1;    % Start frame to analyse
        ited = cst_z+conv_offset_fr;      % End frame to analyse
        if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
            Fz_at(i,:) = F0fz(cst_z-(nfr/2)+1: ...
                cst_z+(nfr/2));
            F_at(i) = nanmean(F0fz(itst:ited));  % Conv centered z F
        else
            continue
        end
        
        % No at distribution
        v_on_nix = find(v_on_n == vis-2); % No convergence epochs of same visual stimulus
        F_on_nrsp = NaN(length(v_on_nix),1); % No response vector
        for j = 1:length(v_on_nix)
            vix = v_on_nix(j);
            rep2 = r_on_n(vix);
            ze = find(v0 == vis & r0 == rep2);
            
            if ~isempty(ze)
                itst = nfr*(ze-1)+cstfr-conv_offset_fr+1;
                ited = nfr*(ze-1)+cstfr+conv_offset_fr;
                if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
                    Fz_at_nrsp{i}(j,:) = F0fz(nfr*(ze-1)+cstfr-(nfr/2)+1: ...
                        nfr*(ze-1)+cstfr+(nfr/2));
                    F_on_nrsp(j) = nanmean(F0fz(itst:ited));
                end
            end
        end
        F_at_nrsp_all{i} = F_on_nrsp;
        Z_at(i) = (F_at(i)-nanmean(F_on_nrsp))/nanstd(F_on_nrsp);   % Z-value
        Z_null_at(i) = (nullmean-nanmean(F_on_nrsp))/nanstd(F_on_nrsp);   % Z-value
    end
        
    %% Write to structure 
    
    lconvat = diff(eyeDelta_at,1,2) < 0;    % Left convergences       
    
%     gmConvMod_at.rois(roi).Fz_at = Fz_at;    
%     gmConvMod_at.rois(roi).Fz_nrsp_at = Fz_at_nrsp;
%     gmConvMod_at.rois(roi).F_at = F_at;
%     gmConvMod_at.rois(roi).F_at_nrsp = F_at_nrsp_all;
%     gmConvMod_at.rois(roi).Z_at = Z_at;
%     gmConvMod_at.rois(roi).Z_null_at = Z_null_at;
%     gmConvMod_at.rois(roi).Lconv_at = lconvat;
    gmConvMod_at.rois(roi).Z_at_Lconv_Mn = nanmean(Z_at(lconvat));
    gmConvMod_at.rois(roi).Z_at_Rconv_Mn = nanmean(Z_at(~lconvat));
    gmConvMod_at.rois(roi).Z_null_at_Lconv_Mn = nanmean(Z_null_at(lconvat));
    gmConvMod_at.rois(roi).Z_null_at_Rconv_Mn = nanmean(Z_null_at(~lconvat));
     
end

%% Save

waitbar(1,h,'Saving')
save(fullfile(datadir,'gmConvMod_at'),'gmConvMod_at');
close(h)

end