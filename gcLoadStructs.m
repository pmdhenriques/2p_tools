function [varargout] = gcLoadStructs(datadir,varargin)
% Loads gc structures added as multiple character input (eg.,
% 'gm','gmv','gmranat').
%
% Pedro Henriques, Aug 2017

if nargin >= 2
    k = 1;
    for i = 2:nargin
        arg = varargin{i-1};
        fprintf('Loading %s\n',arg)
        if exist(fullfile(datadir,[arg '.mat']),'file')
            load(fullfile(datadir,arg));
            varargout{k} = eval(arg);
            k = k+1;
        else
            fprintf('No %s struct in datadir...\n',arg);
        end
    end
end

end