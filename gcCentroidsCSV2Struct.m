function [] = gcCentroidsCSV2Struct(datadir)
% Adds transformed centroids to gmranat data structure.
% Requires a CSV file of the transformed centroids ran using the
% AntsApplyTransformToPoints funtion.
%
% Pedro Henriques, Aug 2017

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
end

disp('Loading gmranat')
load(fullfile(datadir,'gmranat.mat'));
[csvname,csvdir] = uigetfile([datadir '\*.csv'], ...
    'Select CSV of transformed centroids');

disp('Reading CSV file')
Cref = readtable(fullfile(csvdir,csvname));

reffields = fieldnames(gmranat.z(1).STATScrop);
reffields = reffields(contains(reffields,'Centroid_'));
if ~isempty(reffields)
    nrefs = size(reffields,1);
    refs = cell(nrefs,1);
    for r = 1:nrefs
        refspp = strsplit(reffields{r},'Centroid_');
        refs{r} = refspp{2};
    end
    
    [Refidx,overwrite] = listdlg('PromptString','Overwite ref?:',...
        'SelectionMode','single',...
        'OKString','Yes','CancelString','No', ...
        'ListString',refs(:));
    if overwrite
        refname = refs{Refidx};
    else
        refname = inputdlg('Reference name');
        refname = refname{1};
    end
else
    refname = inputdlg('Reference name');
    refname = refname{1};
end

nzs = size(gmranat.z,2);

k = 1;
for z = 1:nzs
    for c = 1:size(gmranat.z(z).STATScrop,1)
        gmranat.z(z).STATScrop(c).(['Centroid_' refname]) = [Cref.x(k) Cref.y(k) Cref.z(k)];
        k = k+1;
    end
end

disp('Saving...')
save(fullfile(datadir,'gmranat.mat'),'gmranat');
end
