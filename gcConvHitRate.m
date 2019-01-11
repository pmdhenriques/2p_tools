function [] = gcConvHitRate(datadir)
% ROI Ca2+ activity related to convergence events and calculates each ROIs
% Hit and Shuffle rates of activity for conv events
%
% Pedro Henriques, Jun 2018

%% Initial Variables

ffiltwidth = 5; % Fluorescence filter width
filter1 = ones(1,ffiltwidth)./ffiltwidth;   % Fluorescent filter
conv_offset = 4;    % Convergence offset for convergece window
conv_offset_high = 7;    % Iterations after/before convergence to consider pre/post window
nshuffle = 1000;    % Number of shuffled iteration to compute hr null dist for each ROI
active_thr = 2; % In zscored (std)
not_active_thr = 1;  % In zscored (std)
convfix = 55;   % Conv it in CvR F vectors

%% Begin

h = waitbar(0,'Loading structures');

load(fullfile(datadir,'gmranat'),'gmranat');
load(fullfile(datadir,'gmrxanat'),'gmrxanat');
load(fullfile(datadir,'gm'),'gm');
load(fullfile(datadir,'gmb'),'gmb');
load(fullfile(datadir,'gmbt'),'gmbt');
load(fullfile(datadir,'gmbf'),'gmbf');
load(fullfile(datadir,'gmv'),'gmv');
load(fullfile(datadir,'gmpc'),'gmpc');

% Experiment variables
nrois = size(gmrxanat.roi,2);
nfr = gm.nfr;
ne_t = size(gmv.visstim,1);

% Make CvR struct and add variables
CvR = struct;
CvR.ffiltwidth = ffiltwidth;
CvR.conv_offset = conv_offset;
CvR.nshuffle = nshuffle;
CvR.active_thr = active_thr;
CvR.not_active_thr = not_active_thr;
CvR.fishdir = datadir;

% Loop through all ROIs
k = 1;
for roi = 1:nrois
    waitbar(roi/nrois,h,'Processing rois');
    
    z = gmrxanat.roi(roi).z;    % z pos
    roiix = [gmranat.z(z).STATScrop.tix] == roi;
    centzbb = gmranat.z(z).STATScrop(roiix).Centroid_ZBB;   % ZBB reference centroid
    
    % Epochs to analyse (spots only)
    se = intersect(gm.zindices(z).e, ...
        find(gmv.visstim(:,1) > 77000 & gmv.visstim(:,1) < 82000));
    
    % ROI Ca2+ vector of entire experiment
    [v0,~,r0] = gcPresentationInfo(gm.zindices(z).e,gm,gmv);
    F0 = [];
    for i = 1:length(v0)
        F0 = cat(1,F0,gmrxanat.roi(roi).Vprofiles(v0(i)).zProfiles(r0(i),:)');
    end
    F0f = filtfilt(filter1,1,double(F0));   % zero-phase filter
    F0fz = zscore(F0f); % z-score
    
    % Get actual stimulus times
    vSTt = (gm.trfr+1)*gm.frtime./1000;
    vttime = gmv.visstim(:,6);
    vEDst = vSTt + vttime;
    
    % Bouts during convergences
    btps = [gmbf.b.p];
    btst = [gmbf.b.st];
    vbtix = ([gmbf.b.Vergbout] == 1 | [gmbf.b.Convbout] == 1) & ...
        btst < vEDst(btps)' & ...
        btst > vSTt;    
    vbtn = histcounts(btps(vbtix),0.5:ne_t+0.5);
    
    % Good convergences
    convix = gmb.convergences(:,3) == 1 & ismember(gmb.convergences(:,1),se);
    % Convergences in stim ON time
    convonix = gmb.convergences(:,2) >= vSTt & ...
        gmb.convergences(:,2) <= vEDst(gmb.convergences(:,1));
    
    % Stim ON convergence epochs
    pconv_on = intersect(se,gmb.convergences(convix & convonix,1));
    npconv_on = length(pconv_on);
    % Stim OFF (spontaneous) convergence epochs
    pconv_off = intersect(se,gmb.convergences(convix & ~convonix,1));
    npconv_off = length(pconv_off);
    
    % Epoch info
    [v_on,~,r_on] = gcPresentationInfo(pconv_on,gm,gmv);
    [v_off,~,r_off] = gcPresentationInfo(pconv_off,gm,gmv);
    nps = [length(v_on) length(v_off)];
    
    % Build conv ON and OFF tables
    Xon = table;
    Xoff = table;
    for ctype = 1:2
        for j = 1:nps(ctype)
            if ctype == 1
                vis = v_on(j);
                rep = r_on(j);
                p = pconv_on(j);
            else
                vis = v_off(j);
                rep = r_off(j);
                p = pconv_off(j);
            end
            
            nbts = vbtn(p); % Number of bouts (vergence on)            
            bix = find(vbtix & btps == p);  % gmbf bout index
            % Bout and eye features for conv bout #1
            if ~isempty(bix)
                bix = bix(1);
                tail_vig120 = gmbf.b(bix).vig120;
                tail_max_angl = gmbf.b(bix).max_angl;
                tail_max_TBF = gmbf.b(bix).max_TBF;
                tail_max_vel = gmbf.b(bix).max_vel;
                tail_pre_peak = gmbf.b(bix).pre_peak;
                tail_mean_TBF = gmbf.b(bix).mean_TBF;
                eye_Ldelta = gmbf.b(bix).Ldelta;
                eye_Rdelta = gmbf.b(bix).Rdelta;
            else
                tail_vig120 = nan;
                tail_max_angl = nan;
                tail_max_TBF = nan;
                tail_max_vel = nan;
                tail_pre_peak = nan;
                tail_mean_TBF = nan;
                eye_Ldelta = nan;
                eye_Rdelta = nan;
            end
            
            visp = find(v0 == vis & r0 == rep);
            
            gmpcix = [gmpc.con.pr] == p;
            convit = gmpc.con(gmpcix).fr;
            
            itst = nfr*(visp-1)+convit-nfr/2+1;
            ited = nfr*(visp-1)+convit+nfr/2;
            if itst >= 1 && ited <= length(F0fz)
                Fconv_cent = F0fz(itst:ited);  % Conv centered z F
                
                % Look for conv related activity
                ixx = find(abs(Fconv_cent) >= active_thr);
                ixx2 = find(abs(Fconv_cent) <= not_active_thr);
                if ~isempty(ixx)
                    % Is active around convergence?
                    if any(ismember(ixx, ...
                            convfix-conv_offset:convfix+conv_offset))
                        isconvactive = 1;
                        cixx = findnearest(convfix,ixx);
                        baseix = findnearest(ixx(cixx(1)),ixx2,-1);
                        ceilix = findnearest(ixx(cixx(1)),ixx2,1);
                        % Response start
                        if ~isempty(baseix)
                            rsp_st = ixx2(baseix);
                        else
                            rsp_st = ixx(cixx(1));
                        end
                        % Response end
                        if ~isempty(ceilix)
                            rsp_ed = ixx2(ceilix);
                        else
                            rsp_ed = ixx(cixx(1));
                        end
                    else
                        rsp_st = 0;
                        rsp_ed = 0;
                        isconvactive = 0;
                    end
                    
                    % Is active in pre convergence window?
                    if any(ismember(ixx, ...
                            convfix-conv_offset_high-(2*conv_offset)+1:convfix-conv_offset_high))
                        ispreactive = 1;
                    else
                        ispreactive = 0;
                    end
                    
                    % Is active in post convergence window?
                    if any(ismember(ixx, ...
                            convfix+conv_offset_high:convfix+conv_offset_high+(2*conv_offset)-1))
                        ispostactive = 1;
                    else
                        ispostactive = 0;
                    end
                else
                    rsp_st = nan;
                    rsp_ed = nan;
                    isconvactive = 0;
                    ispreactive = 0;
                    ispostactive = 0;
                end
                
                % Stimulus angle at convergence
                frix = find(gmb.p(p).vis_traj(:,2) == convit);
                if ~isempty(frix)
                    st_ang_px = gmb.p(p).vis_traj(frix(1),3);
                    st_ang_deg = spotangle_v2(st_ang_px, ...
                        gm.visgeom.sweepamp,gm.visgeom.amax, ...
                        gm.visgeom.R,gm.visgeom.fish2screen);
                else
                    st_ang_deg = nan;
                end
                
                % Concatenate tables
                if ctype == 1
                    Xon = cat(1,Xon, ...
                        table( ...
                        roi,p,vis,rep,convit,rsp_st,rsp_ed, ...
                        isconvactive,ispreactive,ispostactive, ...
                        st_ang_deg,nbts, ...
                        tail_vig120,tail_max_angl,tail_max_TBF,tail_max_vel, ...
                        tail_pre_peak,tail_mean_TBF,eye_Ldelta,eye_Rdelta, ...
                        {single(Fconv_cent)}));
                else
                    Xoff = cat(1,Xoff, ...
                        table( ...
                        roi,p,vis,rep,convit,rsp_st,rsp_ed, ...
                        isconvactive,ispreactive,ispostactive, ...
                        st_ang_deg,nbts, ...
                        tail_vig120,tail_max_angl,tail_max_TBF,tail_max_vel, ...
                        tail_pre_peak,tail_mean_TBF,eye_Ldelta,eye_Rdelta, ...
                        {single(Fconv_cent)}));
                end
            end
        end
    end
    Xon.Properties.VariableNames{end} = 'Fz';
    Xoff.Properties.VariableNames{end} = 'Fz';
    
    % Determine shuffled hit rates
    hr_shuff = NaN(nshuffle,2);
    npconvtypes = [npconv_on npconv_off];
    for ctype = 1:2
        for i = 1:nshuffle
            itst_s = randi(length(F0)-nfr,npconvtypes(ctype),1);
            ia = false(npconvtypes(ctype),1);
            for j = 1:npconvtypes(ctype)
                itst = itst_s(j);
                ited = itst+nfr-1;
                Fconv_cent_s = F0fz(itst:ited);  % Conv centered z F
                ixx = find(abs(Fconv_cent_s) >= active_thr);
                if any(ismember(ixx, ...
                        convfix-conv_offset:convfix+conv_offset))
                    ia(j) = 1;
                end
            end
            hr_shuff(i,ctype) = sum(ia)/npconvtypes(ctype);
        end
    end
    CvR.rois(k).roi = roi;
    CvR.rois(k).z = z;
    CvR.rois(k).cent_zbb = centzbb;
    CvR.rois(k).Xon = Xon;
    CvR.rois(k).Xoff = Xoff;
    
    %% Hit rates
    
    % Conv ON
    CvR.rois(k).hitrate_on_conv = sum(Xon.isconvactive)/size(Xon,1);
    CvR.rois(k).hitrate_off_conv = sum(Xoff.isconvactive)/size(Xoff,1);
    % Conv ON and no pre
    CvR.rois(k).hitrate_on_conv_only = sum(Xon.isconvactive & ...
        ~Xon.ispreactive)/size(Xon,1);
    CvR.rois(k).hitrate_off_conv_only = sum(Xoff.isconvactive & ...
        ~Xoff.ispreactive)/size(Xoff,1);
    % Pre conv
    CvR.rois(k).hitrate_on_pre = sum(Xon.ispreactive)/size(Xon,1);
    CvR.rois(k).hitrate_off_pre = sum(Xoff.ispreactive)/size(Xoff,1);
    % Post conv
    CvR.rois(k).hitrate_on_post = sum(Xon.ispostactive)/size(Xon,1);
    CvR.rois(k).hitrate_off_post = sum(Xoff.ispostactive)/size(Xoff,1);
    % Shuffle hr
    CvR.rois(k).hitrate_on_shuffle_mean = nanmean(hr_shuff(:,1));
    CvR.rois(k).hitrate_on_shuffle_std = nanstd(hr_shuff(:,1));
    CvR.rois(k).hitrate_off_shuffle_mean = nanmean(hr_shuff(:,2));
    CvR.rois(k).hitrate_off_shuffle_std = nanstd(hr_shuff(:,2));
    
    %% Z values
    
    CvR.rois(k).Z_on_conv = (CvR.rois(k).hitrate_on_conv - ...
        CvR.rois(k).hitrate_on_shuffle_mean) / ...
        CvR.rois(k).hitrate_on_shuffle_std;
    CvR.rois(k).Z_on_conv_only = (CvR.rois(k).hitrate_on_conv_only - ...
        CvR.rois(k).hitrate_on_shuffle_mean) / ...
        CvR.rois(k).hitrate_on_shuffle_std;
    CvR.rois(k).Z_on_pre = (CvR.rois(k).hitrate_on_pre - ...
        CvR.rois(k).hitrate_on_shuffle_mean) / ...
        CvR.rois(k).hitrate_on_shuffle_std;
    CvR.rois(k).Z_on_post = (CvR.rois(k).hitrate_on_post - ...
        CvR.rois(k).hitrate_on_shuffle_mean) / ...
        CvR.rois(k).hitrate_on_shuffle_std;
    
    CvR.rois(k).Z_off_conv = (CvR.rois(k).hitrate_off_conv - ...
        CvR.rois(k).hitrate_off_shuffle_mean) / ...
        CvR.rois(k).hitrate_off_shuffle_std;
    CvR.rois(k).Z_off_conv_only = (CvR.rois(k).hitrate_off_conv_only - ...
        CvR.rois(k).hitrate_off_shuffle_mean) / ...
        CvR.rois(k).hitrate_off_shuffle_std;
    CvR.rois(k).Z_off_pre = (CvR.rois(k).hitrate_off_pre - ...
        CvR.rois(k).hitrate_off_shuffle_mean) / ...
        CvR.rois(k).hitrate_off_shuffle_std;
    CvR.rois(k).Z_off_post = (CvR.rois(k).hitrate_off_post - ...
        CvR.rois(k).hitrate_off_shuffle_mean) / ...
        CvR.rois(k).hitrate_off_shuffle_std;
    
    k = k+1;
end

% Save CvR struct
waitbar(1,h,'Saving...');
save(fullfile(datadir,'CvR'),'CvR');
close(h)
end