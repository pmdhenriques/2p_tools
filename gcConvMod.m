function [] = gcConvMod(datadir)
% Computes the convergence modulation for each ROI by computing a list of
% Z-values for each convergence during stim ON epochs with a distribution
% during stim OFF epochs, and testing against a null distribution.
%
% Runs seperate tests for Left and Right convergences.
%
% Pedro Henriques, Jun 2018

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
nshuff = 1000;    % Number of shuffled iteration to compute hr null dist for each ROI
nfr = gm.nfr;   % Number of epoch frames

gmConvMod.settings.restrictToPrey = restrictToPrey;
gmConvMod.settings.ffiltwidth = ffiltwidth;
gmConvMod.settings.conv_offset_fr = conv_offset_fr;
gmConvMod.settings.nshuff = nshuff;
gmConvMod.settings.rundate = str2double(datestr(datetime('today'),'yymmdd'));

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

% Bouts during convergences
btps = [gmbf.b.p];
btst = [gmbf.b.st];
bted = [gmbf.b.ed];
vbts = [gmbf.b.Vergbout];
cbts = [gmbf.b.Convbout];

%%  Run the loop

Htt = NaN(nrois,4);
Httpv = NaN(nrois,4);
Htts = NaN(nrois,4);
Hmw = NaN(nrois,4);
Hmwpv = NaN(nrois,4);
Hmws = NaN(nrois,4);

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
    gmConvMod.rois(roi).z = z;
    % Epochs
    gmConvMod.rois(roi).pconv_on = econv_on;
    gmConvMod.rois(roi).pconv_off = econv_off;
    gmConvMod.rois(roi).nopconv = noeconv;
    % Vistim + reps
    gmConvMod.rois(roi).vis_on = [v_on r_on];
    gmConvMod.rois(roi).vis_off = [v_off r_off];
    gmConvMod.rois(roi).vis_on_n = [v_on_n r_on_n];
    
    %%
    
    % ROI Ca2+ vector of entire experiment
    [v0,~,r0] = gcPresentationInfo(gm.zindices(z).e,gm,gmv);
    F0 = [];
    for i = 1:length(v0)
        F0 = cat(1,F0,gmrxanat.roi(roi).Vprofiles(v0(i)).zProfiles(r0(i),:)');
    end
    F0f = filtfilt(ffilter,1,double(F0));   % zero-phase filter
    F0fz = (F0f-mode(F0f))./std(F0f); % z-score
    F0z = (F0-mode(F0))./std(F0); % z-score
    
    % Null distribution
    Null = NaN(nshuff,1);
    nullstfr = randi(length(F0)-(2*conv_offset_fr),nshuff,1);
    for i = 1:nshuff
        stfr = nullstfr(i);
        Null(i) = nanmean(F0fz(stfr:stfr+(2*conv_offset_fr)-1));  % Conv centered z F
    end
    nullmean = nanmean(Null);
    
    % Stim on convergeces
    Fz_on = NaN(neconv_on,nfr,'single');
    Fz_on_nrsp = cell(neconv_on,1);
    Fz_on_nof = NaN(neconv_on,nfr,'single');
    Fz_on_nrsp_nof = cell(neconv_on,1);
    F_on_nrsp_all = cell(neconv_on,1);
    F_on = NaN(neconv_on,1); % Mean fluorescence list
    Z_on = NaN(neconv_on,1); % Z list (distance to no response trials distribution)
    Z_null_on = NaN(neconv_on,1);   % Null distribution means
    eyeDelta_on = NaN(neconv_on,1);
    eyePrePost_on = NaN(neconv_on,4);
    Cit_on = NaN(neconv_on,1);
    Cst_on = NaN(neconv_on,1);
    nbts_on = NaN(neconv_on,1);
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
        
        % Number of bouts
        btix = find(btps == e & vbts == 1 & cbts == 1 & ...
            btst >= cst-0.1 & bted <= vEDst(e)+0.1);
        if length(btix) > 1
            btvixx = diff(vbts(btix(1):btix(end)));
            if any(btvixx < 0)
                btvixxx = find(btvixx < 0);
                nbts_on(i) = btvixxx(1);
            else
                nbts_on(i) = length(btix);
            end
        else
            nbts_on(i) = length(btix);
        end
        
        eyePrePost_on(i,:) = gmb.convergences(cix,7:10);
        eyeDelta_on(i,1) = diff(gmb.convergences(cix,7:8)); % Left eye delta angle
        eyeDelta_on(i,2) = diff(-gmb.convergences(cix,9:10));   % Right eye delta angle
        
        Cit_on(i) = nfr*(ze-1)+cstfr; % Convergence frame (from start of z plane)
        itst = Cit_on(i)-conv_offset_fr+1;    % Start frame to analyse
        ited = Cit_on(i)+conv_offset_fr;      % End frame to analyse
        if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
            Fz_on(i,:) = F0fz(Cit_on(i)-(nfr/2)+1: ...
                Cit_on(i)+(nfr/2));
            Fz_on_nof(i,:) = F0z(Cit_on(i)-(nfr/2)+1: ...
                Cit_on(i)+(nfr/2));
            F_on(i) = nanmean(F0fz(itst:ited));  % Conv centered z F
        end
        
        % No response distribution
        v_on_nix = find(v_on_n == vis); % No convergence epochs of same visual stimulus
        F_on_nrsp = NaN(length(v_on_nix),1); % No response vector
        Fz_on_nrsp{i} = NaN(length(v_on_nix),nfr,'single');
        for j = 1:length(v_on_nix)
            vix = v_on_nix(j);
            rep2 = r_on_n(vix);
            ze = find(v0 == vis & r0 == rep2);
            
            if ~isempty(ze)
                itst = nfr*(ze-1)+cstfr-conv_offset_fr+1;
                ited = nfr*(ze-1)+cstfr+conv_offset_fr;
                if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
                    Fz_on_nrsp{i}(j,:) = F0fz(nfr*(ze-1)+cstfr-(nfr/2)+1: ...
                        nfr*(ze-1)+cstfr+(nfr/2));
                    Fz_on_nrsp_nof{i}(j,:) = F0z(nfr*(ze-1)+cstfr-(nfr/2)+1: ...
                        nfr*(ze-1)+cstfr+(nfr/2));
                    F_on_nrsp(j) = nanmean(F0fz(itst:ited));
                end
            end
        end
        F_on_nrsp_all{i} = F_on_nrsp;
        Z_on(i) = (F_on(i)-nanmean(F_on_nrsp))/nanstd(F_on_nrsp);   % Z-value
        Z_null_on(i) = (nullmean-nanmean(F_on_nrsp))/nanstd(F_on_nrsp);   % Z-value
    end
    
    % Spontaneous convergences
    Fz_off = NaN(neconv_off,nfr,'single');
    Fz_off_nof = NaN(neconv_off,nfr,'single');
    %     Fz_off_nrsp = cell(neconv_off,1);
    F_off_nrsp_all = cell(neconv_off,1);
    F_off_nrsp_all_nof = cell(neconv_off,1);
    F_off = NaN(neconv_off,1); % Mean fluorescence list
    Z_off = NaN(neconv_off,1); % Z list (distance to no response trials distribution)
    Z_null_off = NaN(neconv_off,1);   % Null distribution means
    eyeDelta_off = NaN(neconv_off,1);
    eyePrePost_off = NaN(neconv_off,4);
    Cit_off = NaN(neconv_off,1);
    Cst_off = NaN(neconv_off,1);
    nbts_off = NaN(neconv_off,1);
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
       
        % Number of bouts
        btix = find(btps == e & vbts == 1 & cbts == 1 & ...
            btst >= cst-0.1);
        if length(btix) > 1
            btvixx = diff(vbts(btix(1):btix(end)));
            if any(btvixx < 0)
                btvixxx = find(btvixx < 0);
                nbts_off(i) = btvixxx(1);
            else
                nbts_off(i) = length(btix);
            end
        else
            nbts_off(i) = length(btix);
        end
        
        eyePrePost_off(i,:) = gmb.convergences(cix,7:10);
        eyeDelta_off(i,1) = diff(gmb.convergences(cix,7:8)); % Left eye delta angle
        eyeDelta_off(i,2) = diff(-gmb.convergences(cix,9:10));   % Right eye delta angle
        
        Cit_off(i) = nfr*(ze-1)+cstfr;
        itst = Cit_off(i)-conv_offset_fr+1;
        ited = Cit_off(i)+conv_offset_fr;
        if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
            Fz_off(i,:) = F0fz(Cit_off(i)-(nfr/2)+1: ...
                Cit_off(i)+(nfr/2));
            Fz_off_nof(i,:) = F0z(Cit_off(i)-(nfr/2)+1: ...
                Cit_off(i)+(nfr/2));
            F_off(i) = nanmean(F0fz(itst:ited));
        end
        
        % No response distribution
        v_off_nix = find(r_off ~= rep); % all off response epochs excluding current
        F_off_nrsp = NaN(length(v_off_nix),1); % No response vector
        F_off_nrsp_nof = NaN(length(v_off_nix),1); % No response vector
        for j = 1:length(v_off_nix)
            vix = v_off_nix(j);
            rep2 = r_off(vix);
            ze = find(v0 == vis & r0 == rep2);
            
            if ~isempty(ze)
                itst = nfr*(ze-1)+cstfr-conv_offset_fr+1;
                ited = nfr*(ze-1)+cstfr+conv_offset_fr;
                if itst >= nfr/2 && ited <= length(F0fz)-nfr/2
                    %                     Fz_off_nrsp{i}(j,:) = F0fz(nfr*(ze-1)+cstfr-(nfr/2)+1: ...
                    %                         nfr*(ze-1)+cstfr+(nfr/2));
                    F_off_nrsp(j) = nanmean(F0fz(itst:ited));
                    F_off_nrsp_nof(j) = nanmean(F0z(itst:ited));
                end
            end
            
        end
        F_off_nrsp_all{i} = F_off_nrsp;
        F_off_nrsp_all_nof{i} = F_off_nrsp_nof;
        Z_off(i) = (F_off(i)-nanmean(F_off_nrsp))/nanstd(F_off_nrsp);   % Z-value
        Z_null_off(i) = (nullmean-nanmean(F_off_nrsp))/nanstd(F_off_nrsp);   % Z-value
    end
    
    %% Compare Z distributions
    
    lconvon = diff(eyeDelta_on,1,2) < 0;    % Left convergences
    lconvoff = diff(eyeDelta_off,1,2) < 0;  % Right convergences
    
    % T-tests for each convergence type
    if length(Z_on) > 1 && length(Z_null_on) > 1
        if sum(lconvon) > 1
            [Htt(roi,1),Httpv(roi,1),~,ttstats] = ttest2(Z_on(lconvon),Z_null_on);
            [Hmwpv(roi,1),Hmw(roi,1),mwstats] = ranksum(Z_on(lconvon),Z_null_on);
            Htts(roi,1) = ttstats.tstat;
            Hmws(roi,1) = mwstats.ranksum;
        end
        if sum(~lconvon) > 1
            [Htt(roi,2),Httpv(roi,2),~,ttstats] = ttest2(Z_on(~lconvon),Z_null_on);
            [Hmwpv(roi,2),Hmw(roi,2),mwstats] = ranksum(Z_on(~lconvon),Z_null_on);
            Htts(roi,2) = ttstats.tstat;
            Hmws(roi,2) = mwstats.ranksum;
        end
    end
    if length(Z_off) > 1 && length(Z_null_off) > 1
        if sum(lconvoff) > 1
            [Htt(roi,3),Httpv(roi,3),~,ttstats] = ttest2(Z_off(lconvoff),Z_null_off);
            [Hmwpv(roi,3),Hmw(roi,3),mwstats] = ranksum(Z_off(lconvoff),Z_null_off);
            Htts(roi,3) = ttstats.tstat;
            Hmws(roi,3) = mwstats.ranksum;
        end
        if sum(~lconvoff) > 1
            [Htt(roi,4),Httpv(roi,4),~,ttstats] = ttest2(Z_off(~lconvoff),Z_null_off);
            [Hmwpv(roi,4),Hmw(roi,4),mwstats] = ranksum(Z_off(~lconvoff),Z_null_off);
            Htts(roi,4) = ttstats.tstat;
            Hmws(roi,4) = mwstats.ranksum;
        end
    end
    
    gmConvMod.rois(roi).Fz_on = Fz_on;
    gmConvMod.rois(roi).Fz_on_nof = Fz_on_nof;
    gmConvMod.rois(roi).Fz_on_nrsp = Fz_on_nrsp;
    gmConvMod.rois(roi).Fz_on_nrsp_nof = Fz_on_nrsp_nof;
    gmConvMod.rois(roi).F_on = F_on;
    gmConvMod.rois(roi).F_on_nrsp = F_on_nrsp_all;
    gmConvMod.rois(roi).Z_on = Z_on;
    gmConvMod.rois(roi).Z_null_on = Z_null_on;
    gmConvMod.rois(roi).Fz_off = Fz_off;
    gmConvMod.rois(roi).Fz_off_nof = Fz_off_nof;
    %     gmConvMod.rois(roi).Fz_off_nrsp = Fz_off_nrsp;
    gmConvMod.rois(roi).F_off = F_off;
    gmConvMod.rois(roi).F_off_nrsp = F_off_nrsp_all;
    gmConvMod.rois(roi).Z_off = Z_off;
    gmConvMod.rois(roi).Z_null_off = Z_null_off;
    gmConvMod.rois(roi).Lconv_on = lconvon;
    gmConvMod.rois(roi).Lconv_off = lconvoff;
    gmConvMod.rois(roi).Z_on_Lconv_Mn = nanmean(Z_on(lconvon));
    gmConvMod.rois(roi).Z_on_Rconv_Mn = nanmean(Z_on(~lconvon));
    gmConvMod.rois(roi).Z_off_Lconv_Mn = nanmean(Z_off(lconvoff));
    gmConvMod.rois(roi).Z_off_Rconv_Mn = nanmean(Z_off(~lconvoff));
    if any(lconvon)
        gmConvMod.rois(roi).Z_on_Lconv_Mx = max(Z_on(lconvon)); else
        gmConvMod.rois(roi).Z_on_Lconv_Mx = nan; end
    if any(~lconvon)
        gmConvMod.rois(roi).Z_on_Rconv_Mx = max(Z_on(~lconvon)); else
        gmConvMod.rois(roi).Z_on_Rconv_Mx = nan;end
    if any(lconvoff)
        gmConvMod.rois(roi).Z_off_Lconv_Mx = max(Z_off(lconvoff)); else
        gmConvMod.rois(roi).Z_off_Lconv_Mx = nan; end
    if any(~lconvoff)
        gmConvMod.rois(roi).Z_off_Rconv_Mx = max(Z_off(~lconvoff)); else
        gmConvMod.rois(roi).Z_off_Rconv_Mx = nan; end
    gmConvMod.rois(roi).Cit_on = Cit_on;
    gmConvMod.rois(roi).Cit_off = Cit_off;
    gmConvMod.rois(roi).Cst_on = Cst_on;
    gmConvMod.rois(roi).Cst_off = Cst_off;
    gmConvMod.rois(roi).nbts_on = nbts_on;
    gmConvMod.rois(roi).nbts_off = nbts_off;    
    gmConvMod.rois(roi).eyeDelta_on = eyeDelta_on;
    gmConvMod.rois(roi).eyeDelta_off = eyeDelta_off;
    gmConvMod.rois(roi).eyePrePost_on = eyePrePost_on;
    gmConvMod.rois(roi).eyePrePost_off = eyePrePost_off;    
end

%% Multiple comparisons correction (Benjamini & Hochberg/Yekutieli false discovery rate)

alpha = 0.05;

[~, ~, ~, Httpv_corr] = fdr_bh(Httpv(:),alpha,'pdep');
Httpv_corr = reshape(Httpv_corr,size(Httpv));
Htt_corr = Httpv_corr < alpha;

[~, ~, ~, Hmwpv_corr] = fdr_bh(Hmwpv(:),alpha,'pdep');
Hmwpv_corr = reshape(Hmwpv_corr,size(Hmwpv));
Hmw_corr = Hmwpv_corr < alpha;

%% Write to structure

% columns: Conv_ON_L; Conv_ON_R; Conv_OFF_L; Conv_OFF_R

gmConvMod.H_tt = Htt;
gmConvMod.H_tt_pv = Httpv;
gmConvMod.H_tt_ts = Htts;
gmConvMod.H_tt_corr = Htt_corr;
gmConvMod.H_tt_pv_corr = Httpv_corr;

gmConvMod.H_mw = Hmw;
gmConvMod.H_mw_pv = Hmwpv;
gmConvMod.H_mw_ts = Hmws;
gmConvMod.H_mw_corr = Hmw_corr;
gmConvMod.H_mw_pv_corr = Hmwpv_corr;

waitbar(1,h,'Saving')
save(fullfile(datadir,'gmConvMod'),'gmConvMod','-v7.3');
close(h)

end