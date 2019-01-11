function [] = gcXgc_allz(datadir)

load(fullfile(datadir,'gm'),'gm')
load(fullfile(datadir,'gmrxanat'),'gmrxanat');

if gm.frtime*gmrxanat.trfr>6000
    Flead = round(6000/gm.frtime); % 6 s baseline for estimating F
else
    Flead = gmrxanat.trfr-1;
end
lpfiltwidth = round(1650./gm.frtime); % 1.65 s box filter
if lpfiltwidth==6
    lpfiltwidth = 5;
end

Xgc = single(zeros(length(gmrxanat.roi), gmrxanat.nfr*length(gmrxanat.roi(1).Vprofiles))); % VRVs for all ROIs
for r = 1:length(gmrxanat.roi)
      
    zxx = gm.zindices(gmrxanat.roi(r).z).e;
    if isfield(gm, 'aligntwice')
        zgood = gm.zindices(gmrxanat.roi(r).z).aligngood2;
    else
        zgood = gm.zindices(gmrxanat.roi(r).z).aligngood;
    end
    
    Harrn = [];    
    for v = 1:length(gmrxanat.roi(r).Vprofiles)
        thisvis = gmrxanat.roi(r).Vprofiles(v).vistyp;
        z_v_ix = zxx(ismember(gm.visstim(zxx, 1:5), thisvis, 'rows'));
        good_ix = ismember(z_v_ix, zgood); % good epochs=1, bad epochs=0
        
        Zarr = gmrxanat.roi(r).Vprofiles(v).zProfiles(good_ix, :)'; %all good epochs as columns
        Farr = nanmean(Zarr(gmrxanat.trfr-Flead:gmrxanat.trfr,:), 1); %vector of Fs, for baseline period prior to trigger
        Sarr = nanstd(Zarr(gmrxanat.trfr-Flead:gmrxanat.trfr,:),[],1); %vector of Fs, for baseline period prior to trigger
        dFF = (Zarr-Farr)./Sarr; %dF/F
        dFF(~isfinite(dFF)) = nan;
        mdFF = nanmean(dFF,2); % mean response across presentations

        if ~isempty(mdFF) && ~all(isnan(mdFF))
            Harrn = [Harrn mdFF']; % this is the VRV           
        else
            Harrn = [Harrn zeros(1,gmrxanat.nfr)]; % if response missing, use zeros
        end       
    end
    
    % smoothing
    Harrn = single(filtfilt(ones(1,lpfiltwidth)./lpfiltwidth, 1, double(Harrn)));

    Xgc(r,:) = single(Harrn); % full response profile of each ROI is a row        
end

% rename Xgc
Xgc_allz = Xgc;
clear Xgc

save(fullfile(datadir,'Xgc_allz'),'Xgc_allz');
end
