function [] = gcConvMod_profile(datadir)
% Computes the convergence modulation profile for each ROI by computing a list of
% Z-values for each convergence during stim ON epochs with a distribution
% during stim OFF epochs.
%
%
% Pedro Henriques, Sep 2018

%% Load structures

h = waitbar(0,'Loading structures');

load(fullfile(datadir,'gm'),'gm');
load(fullfile(datadir,'gmv'),'gmv');
load(fullfile(datadir,'gmb'),'gmb');
load(fullfile(datadir,'gmbf'),'gmbf');
load(fullfile(datadir,'gmrxanat'),'gmrxanat');

%%  Setup variables

restrictToPrey = 1; % Restrict visual stimuli to prey-like

ffiltwidth = 5; % Fluorescence filter width
ffilter = ones(1,ffiltwidth)./ffiltwidth;   % Fluorescent filter
conv_offset_s = 2;    % Convergence offset for convergece window (s)
conv_offset_fr = ceil(conv_offset_s/(gm.frtime/1000));  % Convergence offset in frames
nfr = gm.nfr;   % Number of epoch frames

gmConvMod_profile.settings.restrictToPrey = restrictToPrey;
gmConvMod_profile.settings.ffiltwidth = ffiltwidth;
gmConvMod_profile.settings.conv_offset_fr = conv_offset_fr;
gmConvMod_profile.settings.rundate = str2double(datestr(datetime('today'),'yymmdd'));

%%

% Visual stimuli types
if restrictToPrey
    vtypes = gmv.vistypz(:,1) >= 70000 & gmv.vistypz(:,1) <= 90000;
else
    vtypes = 1:size(gmv.vistypz,1);
end
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
    neconv_on = length(econv_on);
    % Stim OFF (spontaneous) convergence epochs
    econv_off = unique(intersect(se,gmb.convergences(gconvix & ~convonix,1)));
    neconv_off = length(econv_off);
    % No convergence stim ON epochs
    noeconv = setdiff(se,econv_on);
    
    % Epoch info
    [v_on,~,r_on] = gcPresentationInfo(econv_on,gm,gmv);
    [v_off,~,r_off] = gcPresentationInfo(econv_off,gm,gmv);
    [v_on_n,~,r_on_n] = gcPresentationInfo(noeconv,gm,gmv);
    
    % Write to structure
    gmConvMod_profile.rois(roi).z = z;
    % Epochs
    gmConvMod_profile.rois(roi).pconv_on = econv_on;
    gmConvMod_profile.rois(roi).pconv_off = econv_off;
    gmConvMod_profile.rois(roi).nopconv = noeconv;
    % Vistim + reps
    gmConvMod_profile.rois(roi).vis_on = [v_on r_on];
    gmConvMod_profile.rois(roi).vis_off = [v_off r_off];
    gmConvMod_profile.rois(roi).vis_on_n = [v_on_n r_on_n];
    
    %%
    
    % ROI Ca2+ vector of entire experiment
    [v0,~,r0] = gcPresentationInfo(gm.zindices(z).e,gm,gmv);
    F0 = [];
    for i = 1:length(v0)
        F0 = cat(1,F0,gmrxanat.roi(roi).Vprofiles(v0(i)).zProfiles(r0(i),:)');
    end
    F0f = filtfilt(ffilter,1,double(F0));   % zero-phase filter
    F0fz = (F0f-mode(F0f))./std(F0f); % z-score
        
    % Stim on convergeces
    F_on = NaN(neconv_on,nfr); % Mean fluorescence list
    Z_on = NaN(neconv_on,nfr); % Z list (distance to no response trials distribution)
    eyeDelta_on = NaN(neconv_on,1);
    eyePrePost_on = NaN(neconv_on,4);
    Cit_on = NaN(neconv_on,1);
    Cst_on = NaN(neconv_on,1);
    for i = 1:neconv_on
        vis = v_on(i);  % visual stimulus
        rep = r_on(i);  % rep
        e = econv_on(i);    % epoch
        
        ze = find(v0 == vis & r0 == rep);   % z-wise epoch
        
        cix = find(gmb.convergences(:,1) == e & gconvix & convonix); % Convergence index
        cix = cix(1); % use only first convergence        
        
        cst = gmb.convergences(cix,2);  % Convergence start (s)
        Cst_on(i) = cst;
        cstfr = round(cst/(gm.frtime/1000));    % Convergence start (fr)        
        
        eyePrePost_on(i,:) = gmb.convergences(cix,7:10);
        eyeDelta_on(i,1) = diff(gmb.convergences(cix,7:8)); % Left eye delta angle
        eyeDelta_on(i,2) = diff(-gmb.convergences(cix,9:10));   % Right eye delta angle
        
        Cit_on(i) = nfr*(ze-1)+cstfr; % Convergence frame (from start of z plane)
        itst = Cit_on(i)-(nfr/2)+1;    % Start frame to analyse
        ited = Cit_on(i)+(nfr/2);      % End frame to analyse
        if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
            F_on(i,:) = F0fz(itst:ited);
        end
        
        % No response distribution
        v_on_nix = find(v_on_n == vis); % No convergence epochs of same visual stimulus
        F_on_nrsp = NaN(length(v_on_nix),nfr); % No response vector
        for j = 1:length(v_on_nix)
            vix = v_on_nix(j);
            rep2 = r_on_n(vix);
            ze = find(v0 == vis & r0 == rep2);
            
            if ~isempty(ze)
                itst = nfr*(ze-1)+cstfr-(nfr/2)+1;
                ited = nfr*(ze-1)+cstfr+(nfr/2);
                if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
                    F_on_nrsp(j,:) = F0fz(itst:ited);
                end
            end
        end
        Z_on(i,:) = (F_on(i,:)-nanmean(F_on_nrsp))./nanstd(F_on_nrsp);   % Z-value
    end
    
    % Spontaneous convergences    
    F_off = NaN(neconv_off,nfr); % Mean fluorescence list
    Z_off = NaN(neconv_off,nfr); % Z list (distance to no response trials distribution)
    eyeDelta_off = NaN(neconv_off,1);
    eyePrePost_off = NaN(neconv_off,4);
    Cit_off = NaN(neconv_off,1);
    Cst_off = NaN(neconv_off,1);
    for i = 1:neconv_off
        vis = v_off(i);  % visual stimulus
        rep = r_off(i);  % rep
        e = econv_off(i);    % epoch
        
        ze = find(v0 == vis & r0 == rep);   % z-wise epoch
                       
        cix = find(gmb.convergences(:,1) == e & gconvix & ~convonix); % Convergence index
        cix = cix(1); % use only first convergence
        cst = gmb.convergences(cix,2);  % Convergence start (s)
        Cst_off(i) = cst;
        cstfr = round(cst/(gm.frtime/1000));    % Convergence start (fr)       
        
        eyePrePost_off(i,:) = gmb.convergences(cix,7:10);
        eyeDelta_off(i,1) = diff(gmb.convergences(cix,7:8)); % Left eye delta angle
        eyeDelta_off(i,2) = diff(-gmb.convergences(cix,9:10));   % Right eye delta angle
        
        Cit_off(i) = nfr*(ze-1)+cstfr;
        itst = Cit_off(i)-(nfr/2)+1;
        ited = Cit_off(i)+(nfr/2);
        if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
            F_off(i,:) = F0fz(itst:ited);
        end
        
        % No response distribution
        v_nix = find(v0 == vis & r0 ~= rep);
        F_off_nrsp = NaN(length(v_nix),nfr); % No response vector        
        for j = 1:length(v_nix)
            ze = v_nix(j);            
            if ~isempty(ze)
                itst = nfr*(ze-1)+cstfr-(nfr/2)+1;
                ited = nfr*(ze-1)+cstfr+(nfr/2);
                if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
                    F_off_nrsp(j,:) = F0fz(itst:ited);
                end
            end            
        end
        Z_off(i,:) = (F_off(i,:)-nanmean(F_off_nrsp))./nanstd(F_off_nrsp);   % Z-value
    end
    
    %% Compare Z distributions
    
    lconvon = diff(eyeDelta_on,1,2) < 0;    % Left convergences
    lconvoff = diff(eyeDelta_off,1,2) < 0;  % Right convergences      
    
    gmConvMod_profile.rois(roi).Z_on = Z_on;
    gmConvMod_profile.rois(roi).Z_off = Z_off;
    gmConvMod_profile.rois(roi).Lconv_on = lconvon;
    gmConvMod_profile.rois(roi).Lconv_off = lconvoff;
    gmConvMod_profile.rois(roi).Z_on_Lconv_Mn = nanmean(Z_on(lconvon,:));
    gmConvMod_profile.rois(roi).Z_on_Rconv_Mn = nanmean(Z_on(~lconvon,:));
    gmConvMod_profile.rois(roi).Z_off_Lconv_Mn = nanmean(Z_off(lconvoff,:));
    gmConvMod_profile.rois(roi).Z_off_Rconv_Mn = nanmean(Z_off(~lconvoff,:));
    gmConvMod_profile.rois(roi).Cit_on = Cit_on;
    gmConvMod_profile.rois(roi).Cit_off = Cit_off;
    gmConvMod_profile.rois(roi).Cst_on = Cst_on;
    gmConvMod_profile.rois(roi).Cst_off = Cst_off;
    gmConvMod_profile.rois(roi).eyeDelta_on = eyeDelta_on;
    gmConvMod_profile.rois(roi).eyeDelta_off = eyeDelta_off;
    gmConvMod_profile.rois(roi).eyePrePost_on = eyePrePost_on;
    gmConvMod_profile.rois(roi).eyePrePost_off = eyePrePost_off;  
    
end

%%

waitbar(1,h,'Saving')
save(fullfile(datadir,'gmConvMod_profile'),'gmConvMod_profile','-v7.3');
close(h)

end