function [im] = im_sub_empoffset(im)
% Return image with subtracted empirical offset

ts = double(im(:));
immode = mode(ts); % find mode: should be center of noise distribution
ts(ts>immode) = []; % remove all values greater than mode
ts = ts-immode; % subtract the mode, to give portion of noise distribution less then zero
ts = [ts;-ts]; % reflect that distribution about y axis
imstd = nanstd(ts); % determine its standard deviation
empoffset = floor(immode); % offset is mode

im = im-empoffset;
end