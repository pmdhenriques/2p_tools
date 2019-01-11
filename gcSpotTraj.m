function [gmv] = gcSpotTraj(datadir)
% Extracts spot trajectories from _traj.xls files
%
% struct is as follows:
% time(s) / 2P frame / target xpos(px) / target ypos(distance from ycenter)
% / target xpos(px) / target ypos (distrance form ycenter)
% -999 is not being presented
%
% Pedro Heriques, UCL
% Apr 18, 2017

disp(['This is gcSpotTraj, processing ' datadir '...'])
h = waitbar(.5,'Loading `gmv` structure...');
load(fullfile(datadir,'gmv.mat'))

if isdir(fullfile(datadir,'traj'))
    files = dir(fullfile(datadir,'traj','*_traj.xls'));
    trajdir = fullfile(datadir,'traj');
else
    files = dir(fullfile(datadir,'*_traj.xls'));
    trajdir = datadir;
end

if ~isempty(files)
    names = NaN(length(files),1);
    for i = 1:length(files)
        nameprt = strsplit(files(i).name,'_');
        names(i) = str2double(nameprt{1});
    end
    
    for f = 1:size(gmv.visstim,1)
        waitbar(f/size(gmv.visstim,1),h,'Processing spot trajectories...')
        ix = find(names == f-1);    % corrected for difference in indexing
        if any(ix)
            gmv.traj(f).flight = dlmread(fullfile(trajdir,files(ix).name));
        else
            gmv.traj(f).flight = [];
        end
    end
    waitbar(1,h,'Saving...')
    save(fullfile(datadir,'gmv.mat'),'gmv','-v7.3');
else
    fprintf('No traj files on %s\n',datadir)
end
close(h)
end