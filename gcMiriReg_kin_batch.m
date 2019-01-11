function [] = gcMiriReg_kin_batch(datadir)
% Requires:
% gm; gmb, gmbt, gmr, gmranat, gmranatdesig, ml_align_filt (dir)

set(0,'DefaultFigureVisible','off')
CIRFtau = 1.1;  % 0.420 (GC6f); 1.1 (GC6s)
use_filtered = 1;
build_type = 'kin';
RegTypz = {'Oculo', 'Tail'};
turbomode = 0;
ROItype = 'anat';

gcMiriBuild_kin(CIRFtau, RegTypz, turbomode, datadir);
gcMiriRun(use_filtered, build_type, datadir);
gcMiriROIRV_kin(ROItype, build_type, datadir);

end