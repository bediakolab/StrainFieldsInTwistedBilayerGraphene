% register_SPgeometry_DS4_Apr12.m

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/DS4');
load('DS4_prefactorfit_objectdata.mat');
m4.scan_stepsize = 0.5;

filteramp = 0.3;
filterrange = 7;
hotspot_suppression = 0;
AA_inflation_factor = 1.5;
edge_threshold = 5;
gaussian_filter_threshold = 0.1;
use_robust_multistart = 0;
AA_amp_threshold = 0.5; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold);

% m4.plotAACircleFitOnCustomColor(0.3,7,0)
% [pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)

m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor);
boundary_inflation_margin = 10;
boundary_trim_margin = 10;
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);


