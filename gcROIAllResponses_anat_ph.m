function [gmrxanat, gmrxanat4] = gcROIAllResponses_anat_ph(datadir, quadoption)

% Computes anatomical ROI responses, grouped by visual stimuli
%
% requires `gmranat`
% quadoption :  true, for quadscan datasets (ScanType `Q1`) will output gmrxanat4
%               with timeseries data at 4x the native 2P frame rate

% Version notes
versionid = 180307;

% v180307:
%       i. outputs 4x upsampled data from quadscan datasets

% 170810:
%       i. response concat instead of stack concat, vectorized mask application


%%
if nargin<2
    quadoption = false;
end
if nargin<1
    if ismac
        datadir = uigetdir('/Volumes/DATA/EXPT/2P_functional');
    else
        datadir = uigetdir('C:\data');
    end
end

%%
disp(['This is gcROIallresponses_anat, processing ' datadir '...'])

h = waitbar(.3,'Loading `gm` structure...');
load([datadir filesep 'gm.mat']);
waitbar(.6,h,'Loading `gmv` structure...');
load([datadir filesep 'gmv.mat']);
waitbar(.9,h,'Loading `gmranat` structure...');
load([datadir filesep 'gmranat.mat']);

%%

% Check for temporary files
if exist(fullfile(datadir,'gmrxanat_tmp.mat'),'file')
    disp('Loading temporary files')
    tmp = 1;
    load(fullfile(datadir,'gmrxanat_tmp.mat'))
    if quadoption
        load(fullfile(datadir,'gmrxanat4_tmp.mat'))
    end
else
    tmp = 0;
end

%%

gmrxanat.name = gm.name;
gmrxanat.version = versionid;
gmrxanat.trfr = gm.trfr;
gmrxanat.nfr = gm.nfr;
gmrxanat.frtime = gm.frtime;

nvis = size(gmv.vistypz,1);

if isfield(gm, 'aligntwice')
   gmrxanat.maptype = 'ml_align_2';
else
    gmrxanat.maptype = 'ml_align';
end

if ~quadoption || ~isfield(gm, 'ScanType') || ~strcmp(gm.ScanType, 'Q1') || rem(gm.yDef, 8)~=0
    quadoption = false;
    fprintf(2, 'Not using quadoption... \n')
else
    fprintf(2, 'Using quadoption... \n')
end

if quadoption
    gmrxanat4.name = gm.name;
    gmrxanat4.version = versionid;
    gmrxanat4.upsamplefactor = 4; % Always 4 for quad scanning (Q1)!
    gmrxanat4.trfr = gm.trfr*4;
    gmrxanat4.nfr = gm.nfr*4;
    gmrxanat4.frtime = gm.frtime./4;
    
    [scan1, scan2, scan3, scan4] = getquadlines(gm);
    blankim = zeros(gm.yDef, gm.xDef);
    scanmaskz = [];
    for s = 1:4
        switch s
            case 1
                scanlines = scan1;
            case 2
                scanlines = scan2;
            case 3
                scanlines = scan3;
            case 4
                scanlines = scan4;
        end
        tim = blankim;
        tim(scanlines, :) = 1;
        scanmaskz(:,:,s) = tim;
    end
end

%%
 try
     parpool local
 end
%%

% total cell count and z initialization
if tmp
    tcc = size(gmrxanat.roi,2);
    zinit = gmrxanat.roi(end).z+1;
else
    tcc = 0; 
    zinit = 1;
end

for z = zinit:length(gm.zrange)

    % set-up rois
    goodcellz = setdiff(gmranat.z(z).allIdx, gmranat.z(z).excludeIdx);
    ncellz = length(goodcellz);
    maskz = false(size(gmranat.z(z).L,1), size(gmranat.z(z).L,2), ncellz);
    for x = 1:ncellz
        r = goodcellz(x);
        gmrxanat.roi(tcc+x).z = z;
        gmrxanat.roi(tcc+x).label = r;
        for v = 1:nvis
            gmrxanat.roi(tcc+x).Vprofiles(v).zProfiles = single([]);
            gmrxanat.roi(tcc+x).Vprofiles(v).meanprofile = single([]);
            gmrxanat.roi(tcc+x).Vprofiles(v).vistyp = gmv.vistypz(v,:);
        end
        if quadoption
            gmrxanat4.roi(tcc+x).z = z;
            gmrxanat4.roi(tcc+x).label = r;
            for v = 1:nvis
                gmrxanat4.roi(tcc+x).Vprofiles(v).zProfiles = single([]);
                gmrxanat4.roi(tcc+x).Vprofiles(v).meanprofile = single([]);
                gmrxanat4.roi(tcc+x).Vprofiles(v).vistyp = gmv.vistypz(v,:);
            end
        end
        maskz(:,:,x) = gmranat.z(z).L==r;
    end
    
    zxx = gm.zindices(z).e;
    vxx = gm.visstim(zxx, 1:5);
 
    for v = 1:nvis
        waitbar(v/nvis, h, ['Computing ROI responses for z = ' num2str(z) '/' num2str(length(gm.zrange)) '; v = ' num2str(v) '/' num2str(nvis) '...']);
        z_v_ix = zxx(ismember(vxx, gmv.vistypz(v,:), 'rows'));
        z_v_ix = z_v_ix(:);
        
        allcellzzProfiles = [];
        allcellzzProfiles4 = [];
        
        for p = z_v_ix'
            load([datadir filesep gmrxanat.maptype filesep num2str(p-1) '.mat']);
            thisstk = uint16(im_align);
            frms = size(thisstk, 3);
            tic
            stimProfiles = single(nan(ncellz, frms));
            parfor x = 1:ncellz  % parallel process improves speed ~2x (JNCS_B machine)
               mask = uint16(maskz(:,:,x));
               maskarea = single(sum(mask(:)));
               masked = bsxfun(@times, thisstk, mask);
               mask_seq = sum(sum(masked))./maskarea;
               pxts = squeeze(mask_seq);
               stimProfiles(x,:) = pxts;
            end
            allcellzzProfiles = [allcellzzProfiles stimProfiles];
            t1 = toc;
            
            tic
            if quadoption
                stimProfiles4 = single(nan(ncellz, frms*4));     
                for s = 1:4
                    scanmask = uint16(scanmaskz(:,:,s));
                    scanixx = [s:4:4*(frms-1)+s];
                    parfor x = 1:ncellz
                        mask = uint16(maskz(:,:,x));
                        mask = bsxfun(@times,mask,scanmask);
                        maskarea = single(sum(mask(:)));
                        masked = bsxfun(@times, thisstk, mask);
                        mask_seq = sum(sum(masked))./maskarea;
                        pxts = squeeze(mask_seq);
                        stimProfiles4(x, scanixx) = pxts;
                    end
                end
                allcellzzProfiles4 = [allcellzzProfiles4 stimProfiles4];
            end
            t2 = toc;
        end
        for x = 1:ncellz % both zProfiles and meanprofile will be single precision
            gmrxanat.roi(tcc+x).Vprofiles(v).zProfiles = reshape(allcellzzProfiles(x,:),[gmrxanat.nfr length(z_v_ix)])'; % note: '
            gmrxanat.roi(tcc+x).Vprofiles(v).meanprofile = nanmean(gmrxanat.roi(tcc+x).Vprofiles(v).zProfiles, 1);
            
            if quadoption
                gmrxanat4.roi(tcc+x).Vprofiles(v).zProfiles = reshape(allcellzzProfiles4(x,:),[gmrxanat4.nfr length(z_v_ix)])'; % note: '
                gmrxanat4.roi(tcc+x).Vprofiles(v).meanprofile = nanmean(gmrxanat4.roi(tcc+x).Vprofiles(v).zProfiles, 1);
            end
        end
    end
    
    tcc = tcc + ncellz;
    
    % Save temporary files
    waitbar(1,h,'Saving tmp data...') % Save at the end
    save([datadir filesep 'gmrxanat_tmp'],'gmrxanat','-v7.3');
    if quadoption
        waitbar(1,h,'Saving tmp quad data...')
        save([datadir filesep 'gmrxanat4_tmp'],'gmrxanat4','-v7.3');
    end

end

%%

% Rename temporary files at the end
movefile([datadir filesep 'gmrxanat_tmp.mat'],...
    [datadir filesep 'gmrxanat.mat']);
if quadoption
    movefile([datadir filesep 'gmrxanat4_tmp.mat'],...
    [datadir filesep 'gmrxanat4.mat']);
end

%%
try
    parpool close
end
disp('Done.')
close(h);

%% make figure
% % 
% for r=1%1:length(gmrx.roi)
%     
%     visperfig=5; %5 vistypes per figure
%     ymax=150;
%     
%     numv=length(gmrx.roi(r).Vprofiles);
%     numf=ceil(numv/visperfig); 
%     
%     v=0;
%     
%     for f=1:numf;
%         figure('Name',[gmrx.name ' ROI ' num2str(r) ' vProfiles' num2str(f) '/' num2str(numf)]);
%     
%         if (numv-v)>=visperfig;
%             thispanelnumber=visperfig;
%         else
%             thispanelnumber=numv-v;
%         end;
% 
%         sp=1;
%             for p=1:thispanelnumber;
%                 v=v+1;
%                 
%                 %Visstim ON/OFF etc
%                 fon=gmrx.trfr+2;
% 
%                 if gm.viscategory==2
%                     if any(floor(gmrx.roi(r).Vprofiles(v).vistyp)==[100 101]) % 3s + 0 padtime = 3s
%                         foff=gm.trfr+2+ceil(3000./gm.frtime)-1;
%                         von=fon;
%                         voff=foff;
%                     elseif any(floor(gmrx.roi(r).Vprofiles(v).vistyp)==[6 7 16 17]) % 6s + 2*2s padtime = 10s
%                         foff=gm.trfr+2+ceil(10000./gm.frtime)-1;
%                         von=gm.trfr+2+ceil(2000./gm.frtime)-1;
%                         voff=gm.trfr+2+ceil(8000./gm.frtime)-1;
%                     elseif any(floor(gmrx.roi(r).Vprofiles(v).vistyp)==[612 712 1612 1712]) % 12s + 2*2s padtime = 16s
%                         foff=gm.trfr+2+ceil(16000./gm.frtime)-1;
%                         von=gm.trfr+2+ceil(2000./gm.frtime)-1;
%                         voff=gm.trfr+2+ceil(14000./gm.frtime)-1;
%                     end;
%                 else
%                     if any(floor(gmrx.roi(r).Vprofiles(v).vistyp)==[4 5 100])
%                         foff=gm.trfr+2+ceil(3000./gm.frtime)-1;
%                     elseif any(floor(gmrx.roi(r).Vprofiles(v).vistyp)==[6 7])
%                         foff=gm.trfr+2+ceil(6000./gm.frtime)-1;
%                     end;
%                 end;
% 
%                 % normal subplot
%                 subplot(thispanelnumber,2,sp)
%                 hold on
%                 plot(gmrx.roi(r).Vprofiles(v).zProfiles','Color',[.5 .5 .5]);
%                 plot(gmrx.roi(r).Vprofiles(v).meanprofile,'Color','k','LineWidth',2);
%                 plot(nanmedian(gmrx.roi(r).Vprofiles(v).zProfiles',2),'Color','b','LineWidth',2);
%                 if gm.viscategory==2
%                     line([fon von voff foff;fon von voff foff],[0 0 0 0; ymax ymax ymax ymax],'Color',[.5 .5 1])
%                 else
%                     line([fon foff;fon foff],[0 0; ymax ymax],'Color',[.5 .5 1])
%                 end;
%                 ylim([0 ymax]);
%                 title(['Vis ' num2str(gmrx.roi(r).Vprofiles(v).vistyp)])
%                 % sem subplot
%                 subplot(thispanelnumber,2,sp+1)
%                 hold on
%                 s = std(gmrx.roi(r).Vprofiles(v).zProfiles,0,1)./sqrt(size(gmrx.roi(r).Vprofiles(v).zProfiles,1));
%                 plot(gmrx.roi(r).Vprofiles(v).meanprofile+s,'Color','c','LineWidth',1);
%                 plot(gmrx.roi(r).Vprofiles(v).meanprofile-s,'Color','c','LineWidth',1);
%                 plot(gmrx.roi(r).Vprofiles(v).meanprofile,'Color','k','LineWidth',2);
%                 if gm.viscategory==2
%                     line([fon von voff foff;fon von voff foff],[0 0 0 0; ymax ymax ymax ymax],'Color',[.5 .5 1])
%                 else
%                     line([fon foff;fon foff],[0 0; ymax ymax],'Color',[.5 .5 1])
%                 end;
%                 ylim([0 ymax]);
% 
%                 sp=sp+2;
%             end;
%         
%     end;
% end;

