function [] = gcROIwaveplot(datadir,rois,sep,option1,option2,imtype,gm,gmv,gmb,gmbt,gmrxanat)
% Plots the neuronal responses of input ROIs in either a wave plot or as a
% raster plot. If input has ROIs from multiple planes, a figure is
% generated for each plane.
%
% inputs:
%   sep: waveplot y separation between ROI responses
%   option1: 'TR' - for total responses (option2 must provide TR matrix)
%            'VIS' - for specific visual stimulus responses (option2 is
%            visstiim index)
%   option2: specific to option1
%   imtype: 1 for waveplot; 2 for raster
%
% Pedro Henriques, June 2017

if nargin < 7
    disp('Loading structures')
    load(fullfile(datadir,'gm'),'gm');
    load(fullfile(datadir,'gmv'),'gmv');
    load(fullfile(datadir,'gmb'),'gmb');
    load(fullfile(datadir,'gmbt'),'gmbt');
    load(fullfile(datadir,'gmrxanat'),'gmrxanat');
end

if size(rois,1) > size(rois,2)
    rois = rois';
end

[nzs,~,ic] = unique([gmrxanat.roi(rois).z]);

for i = 1:length(nzs)
    z = nzs(i);
    tix = find([gmrxanat.roi.z] == z)';
    ix = rois(ic == i)';
    zrois = tix(ismember(tix,ix,'rows'))';
    
    e = [gm.zindices(z).e(1) gm.zindices(z).e(end)];
    col = parula(size(gmv.vistypz,1));
    
    switch option1
        case 'TR'
            figure('name',sprintf('Z = %d',z),'Position',[200 200 1500 700])
            TR = option2;
            t = cat(1,gmbt.p(e(1):e(end)).resampled);
            b = cat(1,gmb.p(e(1):e(end)).resampled);
            
            if imtype == 2
                X = zeros(length(zrois),size(TR,2));
            end
            
            k = 1;
            for roi = zrois
                x = (TR(roi,:)-nanmean(TR(roi,:)))/nanstd(TR(roi,:));
                xf = filtfilt(ones(1,3)/3,1,x);
                
                ax1 = subplot(5,1,1:3);
                if imtype == 2
                    X(k,:) = x;
                else
                    plot(xf+(k-1)*sep,'k')
                    hold on
                end
                k = k+1;
            end
            
            if imtype == 2
                imagesc(X,[0 3]); colormap gray; colormap(flipud(colormap));
                set(ax1,'YDir','normal')
            end
            
            xlim([0 size(TR,2)]);
            ylims = ylim;
            hold off
            
            ax2 = subplot(5,1,4);
            plot(t(:,3));
            ylim([0 70]); xlim([0 size(TR,2)]);
            ax3 = subplot(5,1,5);
            plot(b(:,2)); hold on;
            plot(b(:,3)); hold off
            ylim([-40 40]); xlim([0 size(TR,2)]);
            
            for v = 1:size(gmv.visstim,1)
                vis = gmv.visstim(v,1:5);
                vtypes = ismember(gmv.vistypz(:,1:5),vis,'rows');
                [vST, vED] = gcStimInfo(vis, gm.trfr, gm.frtime, 0);
                xP = [vST vED vED vST]+(v-1)*gm.nfr;
                
                axes(ax1)
                patch(xP,[ylims(1) ylims(1) ylims(2) ylims(2)], ...
                    col(vtypes,:),'EdgeColor','none','FaceAlpha',.3);
                axes(ax2)
                patch(xP,[0 0 70 70], ...
                    col(vtypes,:),'EdgeColor','none','FaceAlpha',.3)
                axes(ax3)
                patch(xP,[-40 -40 40 40], ...
                    col(vtypes,:),'EdgeColor','none','FaceAlpha',.3)
            end
        case 'VIS'
            vtypes = option2;
            for v = vtypes
                figure('name',sprintf('Z = %d; V = %d',z,v),'Position',[200 200 1500 700])
                Lia = find(ismember(gmv.visstim(:,1:5),gmv.vistypz(v,:),'rows'));
                Lia = Lia(Lia >= e(1) & Lia <= e(2));
                
                t = cat(1,gmbt.p(Lia).resampled);
                b = cat(1,gmb.p(Lia).resampled);
                
                if imtype == 2
                    X = zeros(length(zrois),size(t,1));
                end
                
                k = 1;
                for roi = zrois
                    x = gmrxanat.roi(roi).Vprofiles(v).zProfiles;
                    if gmrxanat.nfr == 4*gm.nfr
                        nreps = size(x,1);
                        x = movsum(x(:),[0 3]);
                        x = x(1:4:end);
                        x = reshape(x,nreps,gm.nfr);
                    end
                    
                    x = double(reshape(x',1,numel(x)));
                    x = (x-nanmean(x))/nanstd(x);
                    xf = filtfilt(ones(1,3)/3,1,x);
                    
                    ax1 = subplot(5,1,1:3);
                    if imtype == 2
                        X(k,:) = x;
                    else
                        plot(xf+(k-1)*sep,'k')
                        hold on
                    end
                    k = k+1;
                end
                set(ax1,'Box','off','TickDir','out','XTick',[])
                if imtype == 2
                    imagesc(X,[0 3]); colormap gray; colormap(flipud(colormap));
                    set(ax1,'YDir','normal')
                end
                xlim([0 size(t,1)]);
                ylims = ylim;
                hold off
                
                ax2 = subplot(5,1,4);
                plot(t(:,3));
                ylim([0 70]); xlim([0 length(x)]);
                set(ax2,'Box','off','TickDir','out','XTick',[])
                ax3 = subplot(5,1,5);
                plot(b(:,2)); hold on;
                plot(b(:,3)); hold off
                set(ax3,'Box','off','TickDir','out')
                ylim([-40 40]); xlim([0 length(x)]);
                
                vis = gmv.vistypz(v,1:5);
                [vST, vED] = gcStimInfo(vis, gm.trfr, gm.frtime, 0);
                reps = size(gmrxanat.roi(zrois(1)).Vprofiles(v).zProfiles,1);
                for k = 1:reps
                    xP = [vST vED vED vST]+(k-1)*gm.nfr;
                    
                    axes(ax1)
                    patch(xP,[ylims(1) ylims(1) ylims(2) ylims(2)], ...
                        col(v,:),'EdgeColor','none','FaceAlpha',.3);
                    if k ~= 1
                        line([gm.nfr*(k-1) gm.nfr*(k-1)], [ylims(1) ylims(2)], ...
                            'Color','k','LineStyle','-')
                    end
                    axes(ax2)
                    patch(xP,[0 0 70 70], ...
                        col(v,:),'EdgeColor','none','FaceAlpha',.3)
                    if k ~= 1
                        line([gm.nfr*(k-1) gm.nfr*(k-1)], [0 70], ...
                            'Color','k','LineStyle','-')
                    end
                    axes(ax3)
                    patch(xP,[-40 -40 40 40], ...
                        col(v,:),'EdgeColor','none','FaceAlpha',.3)
                    if k ~= 1
                        line([gm.nfr*(k-1) gm.nfr*(k-1)], [-40 40], ...
                            'Color','k','LineStyle','-')
                    end
                end
            end
    end
end
end