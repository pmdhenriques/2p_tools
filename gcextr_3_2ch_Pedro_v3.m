function gm = gcextr_3_2ch_Pedro_v3(datadir)

%extract data from 2p imaging expts

% v3: uses mode as zero value
% v3_2 changed average to mean instead of 2*mean (which makes little
% sense...) and added support to single channel processing
% v3 - changed to work with red channel alone. Now doesn't need channel num input;

%%

if nargin < 1
    if ismac
        datadir = uigetdir('/Volumes/DATA/EXPT/2P_functional');
    else
        datadir = uigetdir('C:\data');
    end;
end;

if ismac
    dirspp='/';
else
    dirspp='\';
end;

%%
h=waitbar(0.5,'Initializing...');

datadirparts = strsplit(datadir,dirspp);

gm.name=[datadirparts{end-1} '_' datadirparts{end}];

cd([datadir dirspp 'm'])
lst=dir('*.tif');

idz=[]; % the unique numerical id for each stack
zlvl=[]; % the z cocordinate for each stack
for n=1:length(lst);
    partz=strsplit(lst(n).name,'.');
    idz(n)=str2num(partz{1});
    zlvl(n)=str2num(partz{4});
    if n==1
        gm.frtime=str2num(partz{5});
        gm.trfr=str2num(partz{10});
    end;
end;
[~,ix]=sort(idz);
zlvl=zlvl(ix)';
gm.zrange=unique(zlvl);

for z=1:length(gm.zrange);
    gm.zindices(z).e=find(zlvl==gm.zrange(z));
end;

for n=ix;
    gm.pr{n,1}=lst(ix(n)).name;
end;

for z=1:length(gm.zrange)
    theseidx=gm.zindices(z).e;
    
    firstch0=1;
    firstch1=1;
    
    for i=1:length(theseidx);
        p=theseidx(i);
        
        waitbar(p/length(lst),h,['Processing presentation ' num2str(p) ' of ' num2str(length(lst))]);
        
        StackInfo=imfinfo(gm.pr{p},'tiff');
        im=zeros([StackInfo(1).Height,StackInfo(1).Width,length(StackInfo)],'uint16');
        for f=1:length(StackInfo)
            im(:,:,f)=uint16(imread(gm.pr{p}, 'tiff', f));
        end;
        
        % empirical determination of offset
        ts=double(im(:));
        gm.noise(p).mode=mode(ts); % find mode: should be center of noise distribution
        ts(ts>gm.noise(p).mode)=[]; % remove all values greater than mode
        ts=ts-gm.noise(p).mode; % subtract the mode, to give portion of noise distribution less then zero
        ts=[ts;-ts]; % reflect that distribution about y axis
        gm.noise(p).std=nanstd(ts); % determine its standard deviation
        gm.noise(p).empoffset=floor(gm.noise(p).mode+1*gm.noise(p).std); % offset is mode plus 2 standard deviations, rounded to lower integer
        [p gm.noise(p).empoffset]
        
        % subtract offset from im
        im=im-gm.noise(p).empoffset; % in MATLAB uint16 will stop at zero if you subtract from a number a number greater than it.
        
        % determine channel
        partz=strsplit(gm.pr{p},'.');
        
        if str2num(partz{end-1}) == 0 % green channel
            % Build average image
            if firstch0
                zimg=mean(im,3)./length(theseidx);      % was 2*mean for some reason
                firstch0=0;
            else
                zimg=zimg+(mean(im,3)./length(theseidx));   % was 2*mean
            end;
            
        elseif str2num(partz{end-1}) == 1 % red channel
            % Build average image
            if firstch1
                zimg1=mean(im,3)./length(theseidx);     % was 2*mean
                firstch1=0;
            else
                zimg1=zimg1+(mean(im,3)./length(theseidx));     % was 2*mean
            end;
        end;
        
    end;
    
    if exist('zimg','var')
        gm.zimg(z).e=zimg;
        imwrite(uint16(zimg),[datadir dirspp 'zproj_ml.tif'],'TIF','Compression','none','WriteMode','append');
    end
    if exist('zimg1','var')
        gm.zimg1(z).e=zimg1;
        imwrite(uint16(zimg1),[datadir dirspp 'zproj_ml_1.tif'],'TIF','Compression','none','WriteMode','append');
    end
end

gm.nfr=length(StackInfo);

save([datadir dirspp 'gm'],'gm')
%%

close(h);
