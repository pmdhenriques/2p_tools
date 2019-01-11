function [] = gcRotateIm(datadir)
% Rotates all images from imaging experiment so that they conform with the
% standards (anterior left)

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P');
end

h = waitbar(0.5,'Loading gm');
load(fullfile(datadir,'gm'));

for z = 1:size(gm.zimg,2)
    waitbar(z/size(gm.zimg,2),h,sprintf('Rotating Z = %d',z));
    gm.zimg(z).e = rot90(gm.zimg(z).e,2);
    gm.zimg(z).aligngood = rot90(gm.zimg(z).aligngood,2);
    gm.zimg(z).aligngood2 = rot90(gm.zimg(z).aligngood2,2);
end

waitbar(1,h,'Saving...')
save(fullfile(datadir,'gm'),'gm');
close(h)
end