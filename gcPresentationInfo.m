function [v,z,rep] = gcPresentationInfo(p,gm,gmv)
% Returns some handy indices for input stimulus presentation number (p)
%
% Pedro Henriques, Fev 2018

nzs = size(gmv.averageresp.z,2);
np = length(p); % Number of input presentations

z = NaN(np,1);
v = NaN(np,1);
rep = NaN(np,1);

for i = 1:np
    % Find which visstim type
    visstim = gmv.visstim(p(i),1:5);
    v(i) = find(ismember(gmv.vistypz(:,1:5),visstim,'rows'));
    
    % Find z
    for zz = 1:nzs
        ix = find(gm.zindices(zz).e == p(i),1);
        if ~isempty(ix)
            z(i) = zz;
            break
        end
    end
    
    % Find repetition
    e_range = gm.zindices(z(i)).e;
    visp = ismember(gmv.visstim(e_range,1:5),visstim,'rows');
    rep(i) = sum(e_range(visp) <= p(i));
end
end

