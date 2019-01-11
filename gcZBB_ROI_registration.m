function [] = gcZBB_ROI_registration()
%%
% Computes the new registered ROIs by applying the transformation matrix
% spit out by CMTK registration onto the Z-Brain reference brain. Also
% retrieves the anatomical labels for each ROI.
%
% Pedro Henriques, June 2017

dim = [616 1030 420];

folders = uipickfiles('Prompt','Select Folders to analyse', ...
    'FilterSpec','\\128.40.155.187\data2\Bianco_lab\Pedro\2P\NI\NI_Optoma');

h = waitbar(0.5,'Initializing...');
for f = 1:length(folders)
    datadir = folders{f};
    
    waitbar(0.5,h,'Loading `gm` structure...')
    load([datadir filesep 'gm.mat']);
    waitbar(0.5,h,'Loading `gmranat` structure...')
    load([datadir filesep 'gmranat.mat']);
    waitbar(0.5,h,'Loading MaskDatabase...')
    load('ZBBMaskDatabase');
    
    % import affine transformation matrix
    tmatfile = dir(fullfile(datadir,'registration','*csv'));
    if ~isempty(tmatfile)
        tmat = importdata(fullfile(tmatfile.folder,tmatfile.name));
    end
    
    % transpose to meet MATLAB standard
    gmranat.tmat = tmat';
    
    % z steps from range (expect equal interval)
    zstep_aq = diff(gm.zrange(1:2));    
    
    % Inverse to make point transform
    tform = invert(affine3d(gmranat.tmat));
    
    % Loop through z positins
    for z = 1:size(gmranat.z,2)
        waitbar(z/size(gmranat.z,2),h,['Analysing Z = ' num2str(z)])
        
        % Get original centroids
        ref0xy = vertcat(gmranat.z(z).STATScrop.Centroid);
        
        % transform to position in um
        xy0 = bsxfun(@times,ref0xy,[gm.xpx,gm.ypx]);
        xy0(:,3) = z*zstep_aq;
        
        % Transform the centroids to ZBB space
        ref1xy = transformPointsForward(tform,xy0);       % Affine transf to GC6f ref
        
        % Extract mask ids for each centroid
        M = false(size(ref1xy,1),size(ZBBMaskDatabase,2));
        for m = 1:size(ZBBMaskDatabase,2)
            waitbar(m/size(ZBBMaskDatabase,2),h,sprintf('Z = %d; Extracting Mask %d ixs',z,m))
            Mask = reshape(full(ZBBMaskDatabase(:,m)),dim);
            xyz = round(ref1xy);
            M(:,m) = Mask(sub2ind(dim,xyz(:,2),xyz(:,1),xyz(:,3)));
        end
        
        % Write to structure
        for i = 1:size(gmranat.z(z).STATScrop,2)
            gmranat.z(z).STATScrop(i).ref_Centroid = ref1xy(i,:);
            gmranat.z(z).STATScrop(i).ref_Masks = M(i,:);
        end
    end
    
    % Save
    waitbar(1,h,'Saving...')
    save(fullfile(datadir,'gmranat'),'gmranat');
end
close(h)
end