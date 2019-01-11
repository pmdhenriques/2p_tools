
% Dir of mask files
datadir = '\\128.40.155.187\data2\Bianco_lab\Registration\ZBB_browser\anatomy';
files = dir(fullfile(datadir,'*tif'));

scalefactor = 2;    % Scale factor of the ZBB masks

dim1 = [308 515 210];
dim2 = dim1*scalefactor;

%
[y,x,z] = ndgrid( ...
    linspace(1,dim1(1),dim2(1)),...
    linspace(1,dim1(2),dim2(2)),...
    linspace(1,dim1(3),dim2(3)));

ZBBMaskDatabase = logical(sparse(prod(dim2),size(files,1)));

h = waitbar(0.5,'Initializing');
for f = 1:size(files,1)
    waitbar(f/size(files,1),h,'Processing...')
    info = imfinfo(fullfile(datadir,files(f).name));
    im = zeros(dim1,'single');
    
    for i = 1:dim1(3)
        im(:,:,i) = imread(fullfile(datadir,files(f).name),i);
    end
    
    imInt = interp3(im,x,y,z);
    imInt(imInt > 0) = 255;
    imInt = logical(imInt);
    imInt = reshape(imInt,numel(imInt),1);
    imInt = sparse(imInt);
    ZBBMaskDatabase(:,f) = imInt;
    clear im imInt
end
close(h)

% Create names array
ZBBMaskDatabaseNames = cell(size(files,1),1);
for f = 1:size(files,1)
    ZBBMaskDatabaseNames{f} = strtok(files(f).name,'.');
end