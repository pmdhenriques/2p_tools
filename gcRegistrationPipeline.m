% Pipeline for registration of Calcium imaging datasets onto reference brain

datadir = 'Y:\Pedro2\2P\NI_AF7_OT\HuC_H2B_GC6s\180306\f1_2';    % Data directory
rot180 = 1; % Rotate by 180 deg?

%%

zcorr = []; % Any z correction applied to points?

gcIm2Nrrd(datadir,rot180);  % Transform alinged stack to nrrd for registration

% ==> Register original stack to reference (ZBB)

gcCentroids2CSV(datadir,rot180,zcorr);    % Get centroid csv file for ANTs transformation

% ==> Transform csv file using antsApplyTransformToPoints function

gcCentroidsCSV2Struct(datadir);     % Add transformed csv to structure

gcROIaddZBBIdxMasks(datadir);   % Atribute mask ids to ROIs

%% Sanity check

% Get OT rois
[otrois,~] = gcROIextractZBBIdxMask(datadir, ...
    'optic tectum - stratum periventriculare', ...
    'Masks_ZBB');

% Get aNI rois
[anirois,~] = gcROIextractZBBIdxMask(datadir, ...
    'nucleus isthmus - anterior', ...
    'Masks_ZBB');

% Get pNI rois
[pnirois,~] = gcROIextractZBBIdxMask(datadir, ...
    'nucleus isthmus - posterior', ...
    'Masks_ZBB');

load(fullfile(datadir,'gm'),'gm');
load(fullfile(datadir,'gmranat'),'gmranat');

gcROIplotcentroids(catpad(2,otrois,anirois,pnirois),[],gm,gmranat);
