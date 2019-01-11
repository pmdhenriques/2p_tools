function [] = gcROIfeature_vec(datadir,showplots)

if nargin < 2
    showplots = 0;
end

%% Load structures

h = waitbar(0, 'Init');
waitbar(0,h, 'Loading structures')

load(fullfile(datadir,'gm'),'gm');
load(fullfile(datadir,'gmv'),'gmv');
load(fullfile(datadir,'gmb'),'gmb');
load(fullfile(datadir,'gmbf'),'gmbf');
if ~isfield(gmbf.b,'Vergbout')
    gcbft_verg(datadir);
    load(fullfile(datadir,'gmbf'),'gmbf');
end
load(fullfile(datadir,'gmranat'),'gmranat');
if exist(fullfile(datadir,'gmrxanat4.mat'),'file')
    quadscan = 1;
    load(fullfile(datadir,'gmrxanat4'),'gmrxanat4');
    gmrxanat = gmrxanat4; clear gmrxanat4;
else
    quadscan = 0;
    load(fullfile(datadir,'gmrxanat'),'gmrxanat');
end
load(fullfile(datadir,'gmVCNV'),'gmVCNV');
load(fullfile(datadir,'Xgc_all'),'Xgc_all')
load(fullfile(datadir,'gmMRRV_kin'),'gmMRRV')

%%  Determine various convergence epochs

stims = [
    110,0,0,0;
    110,2,0,0;
    80006,0,5,30;
    80007,0,5,30;
    81006,0,5,30;
    81007,0,5,30;
    ];

vtypes = find(ismember(gmv.vistypz(:,1:4),stims,'rows'));

nvis = length(vtypes);
rois = 1:length(gmrxanat.roi);
nrois = length(rois);
nzs = size(gmranat.z,2);

se = find(ismember(gmv.visstim(:,1:4),stims,'rows'));

ne = length(se);
ne_t = size(gmv.visstim,1);

nfr_ori = gm.nfr;
vEDoffset = 15;
vEDoffsett = vEDoffset*gm.frtime./1000;
if quadscan
    frtime = gmrxanat.frtime/1000;
    nfr = gmrxanat.nfr;    
    vEDoffset = vEDoffset*4;
else
    frtime = gm.frtime/1000;
    nfr = gm.nfr;
end

% Get actual stimulus times
vST = repmat(gm.trfr + 2,nvis,1);
vSTt = (gm.trfr+1)*gm.frtime./1000;
vED = gm.trfr + 1 + ceil(1000.*gmv.vistypz(vtypes,5)./gm.frtime);
vEDt = vSTt + gmv.vistypz(vtypes,5);

vttime = gmv.visstim(:,6);
vEDs = gm.trfr + 1 + ceil(1000.*vttime./gm.frtime);
vEDst = vSTt + vttime;
if quadscan
    vST = vST*4;
    vED = vED*4;
    vEDs = vEDs*4;
end

% High vergence number of bouts during stim ON
vbtix = ([gmbf.b.Vergbout] == 1 | [gmbf.b.Convbout] == 1) & ...
    [gmbf.b.st] < vEDst([gmbf.b.p])' & ...
    [gmbf.b.st] > vSTt;

% High vergence bout count per epoch
vbtn = histcounts([gmbf.b(vbtix).p],0.5:ne_t+0.5);

% Good convergences
convix = gmb.convergences(:,3) == 1 & ismember(gmb.convergences(:,1),se);

% Epochs with convergences
pconv = intersect(se,gmb.convergences(convix,1));

% Epochs without convergences
npconv = setdiff(se,pconv);

% Convergences in stim ON time
convonix = gmb.convergences(:,2) >= vSTt & ...
    gmb.convergences(:,2) <= vEDst(gmb.convergences(:,1));

% Epochs with NO convergences in stim ON time
npnconvon = setdiff(se,gmb.convergences(convonix & convix,1));

% Convergences in PRE stim ON time
convpreix = gmb.convergences(:,2) < vSTt;

% Epochs with NO convergences PRE stimulus ON time
pconvpreon = setdiff(se,gmb.convergences(convpreix & convix,1));

% Convergences in stim ON time + offset
convonix_offs = gmb.convergences(:,2) >= vSTt & ...
    gmb.convergences(:,2) <= vEDst(gmb.convergences(:,1))+vEDoffsett;

% Epochs with NO convergences in stim ON time + offset
npnconvon_offs = setdiff(se,gmb.convergences(convonix_offs & convix,1));

% Epochs with convergences in stim ON time
pconvon = intersect(se,gmb.convergences(convix & convonix,1));

% Epochs with convergences outside stim ON time (spontaneous)
pnconvon = intersect(se,gmb.convergences(convix & ~convonix,1));

% Convergence epochs with multiple bout
msbtcp = pconvon(vbtn(pconvon) > 1);

% Convergence epochs with one bout
osbtcp = pconvon(vbtn(pconvon) == 1);

% Get v,z and rep indexes for presentations
[npconv_v,npconv_z,npconv_rep] = gcPresentationInfo(npconv,gm,gmv);
[pnconvon_v,pnconvon_z,pnconvon_rep] = gcPresentationInfo(pnconvon,gm,gmv);
[pconvon_v,pconvon_z,pconvon_rep] = gcPresentationInfo(pconvon,gm,gmv);
[osbtcp_v,osbtcp_z,osbtcp_rep] = gcPresentationInfo(osbtcp,gm,gmv);
[msbtcp_v,msbtcp_z,msbtcp_rep] = gcPresentationInfo(msbtcp,gm,gmv);
[npnconvon_v,npnconvon_z,npnconvon_rep] = gcPresentationInfo(npnconvon,gm,gmv);
[npnconvon_offs_v,npnconvon_offs_z,npnconvon_offs_rep] = gcPresentationInfo(npnconvon_offs,gm,gmv);
[pconvpreon_v,pconvpreon_z,pconvpreon_rep] = gcPresentationInfo(pconvpreon,gm,gmv);

%%  Compute visual response vectors and convergence modulation

filtwidth = 5;
filter1 = ones(1,filtwidth)./filtwidth;

repsratio = NaN(nzs,nvis);
Rmean = zeros(nrois,nvis,2);
MD = zeros(nrois,3,2);
for roi = 1:nrois
    waitbar(roi/nrois,h, 'Processing...')
    z = gmVCNV.ROIinfo(roi,2);
    
    X = cat(1,gmrxanat.roi(roi).Vprofiles(vtypes).zProfiles)';
    Xpre = X(1:vST(1)-1,:);
    baseline = nanmedian(Xpre(:));
    
    M = cell(5,1);
    for v = 1:nvis
        vtype = vtypes(v);
        nreps = size(gmrxanat.roi(roi).Vprofiles(vtype).zProfiles,1);
        
        F = gmrxanat.roi(roi).Vprofiles(vtype).zProfiles';
        F = reshape(interpolate_NaNs(F(:)),size(F));
        dFF = (F-baseline)./baseline;
        dFFf = reshape(filtfilt(filter1,1,dFF(:)), ...
            size(dFF));
        vreps = npnconvon_rep(npnconvon_z == z & npnconvon_v == vtype);
        Rmean(roi,v,1) = nanmean(nanmean(dFFf(vST(v):vED(v)+vEDoffset,vreps)));
        Rmean(roi,v,2) = nanstd(nanmean(dFFf(vST(v):vED(v)+vEDoffset,vreps)));
        
        convreps = pconvon_rep(pconvon_z == z & pconvon_v == vtype);
        obconvreps = osbtcp_rep(osbtcp_z == z & osbtcp_v == vtype);
        mbconvreps = msbtcp_rep(msbtcp_z == z & msbtcp_v == vtype);
        nconvreps = npnconvon_offs_rep(npnconvon_offs_z == z & npnconvon_offs_v == vtype);
        convprereps = pconvpreon_rep(pconvpreon_z == z & pconvpreon_v == vtype);
        
        M{1} = [M{1} ...
            nanmean(dFFf(vST(v):vED(v)+vEDoffset,convreps))];
        M{2} = [M{2} ...
            nanmean(dFFf(vST(v):vED(v)+vEDoffset,obconvreps))];
        M{3} = [M{3} ...
            nanmean(dFFf(vST(v):vED(v)+vEDoffset,mbconvreps))];
        M{4} = [M{4} ...
            nanmean(dFFf(vST(v):vED(v)+vEDoffset,nconvreps))];
        M{5} = [M{5} ...
            nanmean(dFFf(vST(v):vED(v)+vEDoffset,convprereps))];
        
        repsratio(z,v) = length(vreps)/nreps;
    end
    
    % Convergence modulation
    MD(roi,1,1) = nanmean(M{1})-nanmean(M{4});
    if ~isnan(MD(roi,1,1))
        MD(roi,1,2) = ranksum(M{1},M{4});
    else
        MD(roi,1,2) = nan;
    end
    MD(roi,1,3) = MD(roi,1,1)/nanstd(M{4});
    
    % Multiple response modulation
    MD(roi,2,1) = nanmean(M{3})-nanmean(M{2});
    if ~isnan(MD(roi,2,1))
        MD(roi,2,2) = ranksum(M{3},M{2});
    else
        MD(roi,2,2) = nan;
    end
    MD(roi,2,3) = MD(roi,2,1)/nanstd(M{2});
    
    % Pre spontaneous convergence modulation
    MD(roi,3,1) = nanmean(M{5})-nanmean(M{4});
    if ~isnan(MD(roi,3,1))
        MD(roi,3,2) = ranksum(M{5},M{4});
    else
        MD(roi,3,2) = nan;
    end
    MD(roi,3,3) = MD(roi,3,1)/nanstd(M{4});
end

%% Multiple comparisons correction (Benjamini & Hochberg/Yekutieli false discovery rate)

MDcorr = MD;
fdr = 0.05; % false discovery rate

for i = 1:size(MD,2)
    [~, ~, ~, MDcorr(:,i,2)] = fdr_bh(MD(:,i,2),fdr,'pdep');
    MDcorr(MDcorr(:,i,2) > 1,i,2) = 1;
end

%%  Plot modulation maps

clim = 0.5;
if showplots
    cmap = flip(brewermap(256,'RdBu'),1);
    imrange = [0 40];
    
    C = [];
    for z = 1:nzs
        C = cat(1,C,cat(1,gmranat.z(z).STATScrop.Centroid));
    end
    
    for z = 1:nzs
        zix = gmVCNV.ROIinfo(:,2) == z;                
        
        figure('name',sprintf('Z = %d',z), ...
            'position',[206,283,1492,625])
        subplot(1,2,1)
        imagesc(gm.zimg(z).aligngood2,imrange)
        axis('off','square')
        colormap(gca,'gray')
        hold on
        newcolorbar
        scatter(C(zix,1),C(zix,2),[],MDcorr(zix,1,1),'filled')
        axis('off','square')
        colormap(cmap)
        set(gca,'YDir','reverse')
        caxis([-clim clim])
        title('Convergence modulation (GO vs NO-GO trials)')
        
        subplot(1,2,2)
        imagesc(gm.zimg(z).aligngood2,imrange)
        axis('off','square')
        colormap(gca,'gray')
        hold on
        newcolorbar
        scatter(C(zix,1),C(zix,2),[],MDcorr(zix,2,1),'filled')
        axis('off','square')
        colormap(cmap)
        set(gca,'YDir','reverse')
        caxis([-clim clim])
        title('Multiple response modulation')
    end
end

%% Build ROI finger matrix

Xvis = [];
for v = 1:nvis
    vtype = vtypes(v);
    Xvis = cat(2,Xvis, ...
        Xgc_all(:,((vtype-1)*nfr_ori)+1:((vtype-1)*nfr_ori)+nfr_ori));
end

Xreg = gmMRRV.RRV(:,1:17);

Czbb = [];
Czbblab = [];
for z = 1:nzs
    Czbb = cat(1,Czbb,gmranat.z(z).STATScrop.Centroid_ZBB);
    Czbblab = cat(1,Czbblab,gmranat.z(z).STATScrop.Masks_ZBB);
end
Czbblab = logical(Czbblab);

gmROIfvec.Xv = Xvis;
gmROIfvec.Xr = Xreg;
gmROIfvec.reg_names = gmMRRV.Rnamez(1:17);
gmROIfvec.X_zpos = gmVCNV.ROIinfo(:,2);
gmROIfvec.X_Centroids_ZBB = Czbb;
gmROIfvec.X_ZBB_Labels = Czbblab;
gmROIfvec.Rmean = Rmean;
gmROIfvec.Rmax = gmVCNV.ROI_R(:,vtypes);
gmROIfvec.Modulation.MD = MD;
gmROIfvec.Modulation.MDcorr = MDcorr;


waitbar(1,h, 'Saving...')
save(fullfile(datadir,'gmROIfvec'),'gmROIfvec');
close(h)

end

