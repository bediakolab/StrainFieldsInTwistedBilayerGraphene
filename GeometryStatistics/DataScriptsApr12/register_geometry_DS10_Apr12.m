% register_geometry_DS10_Apr12.m

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/DS10');
load('DS10_ringedgauss_prefactorfit_objectdata.mat');
m4.scan_stepsize = 0.5;
% m4.classifyDisplacementConvergence();

filteramp = 0.3;
filterrange = 7;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0);
hotspot_suppression = 0;
AA_inflation_factor = 2;
edge_threshold = 5;
gaussian_filter_threshold = 0.1;
use_robust_multistart = 1;
AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold);

m4.plotAACircleFitOnCustomColor(filteramp,filterrange,0)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)%,[1.2:0.2:2.4])

make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics()

m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor);
boundary_inflation_margin = 10;
boundary_trim_margin = 10;
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
