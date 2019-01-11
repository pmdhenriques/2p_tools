function [] = gcZBtransformation(nt)
%%
% Computes the new registered ROIs by applying the transformation matrix
% spit out by CMTK registration onto the Z-Brain reference brain. Also
% retrieves the anatomical labels for each ROI.
%
% Requires regitration (ofc) and transformation matrix with name
% "tmat.csv".

theta = -90;
getMask = 0;

if nargin < 1
    nt = 1;
end

folders = uipickfiles('Prompt','Select Folders to analyse','FilterSpec','G:\2P\hunting');

h=waitbar(0.5,'Initializing...');
for f = 1:length(folders)
    datadir = folders{f};
    
    waitbar(0.5,h,'Loading `gm` structure...')
    load([datadir filesep 'gm.mat']);
    waitbar(0.5,h,'Loading `gmranat` structure...')
    load([datadir filesep 'gmranat.mat']);
    if getMask
        waitbar(0.5,h,'Loading MaskDatabase...')
        load('MaskDatabase');
    end
        
    tmatfile = dir(fullfile(datadir,'registration','*csv'));
    if ~isempty(tmatfile)
        tmat = importdata(fullfile(tmatfile.folder,tmatfile.name));
    end
    
    gmranat.tmat = tmat';
    zstep_aq = diff(gm.zrange(1:2));
    tform2 = invert(affine3d(gmranat.tmat));        % to referece 
    
    if nt == 2
        tmat2 = importdata([datadir filesep 'registration' filesep filesep 'tmat2.csv']);   % Load ZB tMAT
        gmranat.tmat2 = tmat2';
        tform1 = affine2d([cosd(theta) -sind(theta) 0; ...
            sind(theta) cosd(theta) 0; 0 0 1]);     % 90 deg rotation
        tform3 = invert(affine3d(gmranat.tmat2));       % GC6f to ZBrain
    end
    
    for i = 1:size(gmranat.z,2)
        waitbar(i/size(gmranat.z,2),h,['Analysing Z = ' num2str(i)])
        ref0xy = vertcat(gmranat.z(i).STATScrop.Centroid);     % Get original centroids
        
        xy0 = bsxfun(@times,ref0xy,[gm.xpx,gm.ypx]);
        xy0(:,3) = i*zstep_aq;      % z in um
        
        ref1xy = transformPointsForward(tform2,xy0);       % Affine transf to GC6f ref
        
        if nt == 2
            xy1 = ref1xy;
            xy1(:,2) = 537.91-xy1(:,2);
            xy1 = [transformPointsForward(tform1,xy1(:,1:2)) xy1(:,3)];    % -90deg rotation to ZB ref
            xy1(:,1) = xy1(:,1)+537.91;     % Translate to original axis by adding size of GC6f ref y axis
            xy1(:,2) = 768.44-xy1(:,2);
            xy1(:,3) = (153-xy1(:,3));        % Add z positions
            
            ref2xy = transformPointsForward(tform3,xy1);       % Affine transf to ZB ref
        end
        
        for j = 1:size(gmranat.z(i).STATScrop,2)
            gmranat.z(i).STATScrop(j).ref_Centroid = ref1xy(j,:);
            if nt == 2
                gmranat.z(i).STATScrop(j).ref2_Centroid = ref2xy(j,:);
            end
            if getMask && nt == 2
                progressbar('Zpos','Centroid','MaskID')
                if j == 1
                    Mask_ix = NaN(size(MaskDatabase,2),1);
                end
                for l = 1:size(MaskDatabase,2)
                    Mask = reshape(full(MaskDatabase(:,l)), [height, width, Zs]);
                    Mask_ix(l) = Mask(round(ref2xy(j,2)/0.7980), ...
                        round(ref2xy(j,1)/0.7980), ...
                        round(ref2xy(j,3)));
                    
                    frac3 = l/size(MaskDatabase,2);
                    frac2 = ((j-1) + frac3) / size(gmranat.z(i).STATScrop,1);
                    frac1 = ((i-1) + frac2) / size(gmranat.z,2);
                    progressbar(frac1, frac2, frac3)
                end
                gmranat.z(i).STATScrop(j).Masks = MaskDatabaseNames(logical(Mask_ix));
                progressbar(1)
            end
        end
    end
    save([datadir filesep 'gmranat'],'gmranat','-v7.3');
end
close(h)
end