% register_geometry_DS26_trulyApr12.m
%
% Initial sensation upon playing around with this dataset is that the
% original version may be better behaved.
%
% Nathanael Kazmierczak 05/05/2020

load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/DS26/DS26_newfit_objectdata.mat');
m4.scan_stepsize = 0.5;


filteramp = 0.3;
filterrange = 7;
hotspot_suppression = 0;
AA_inflation_factor = 2;
edge_threshold = 10;
gaussian_filter_threshold = 0.15;
use_robust_multistart = 1;
AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold);
m4.plotAACircleFitOnCustomColor(filteramp,filterrange,0)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1,1:0.2:2.5)

make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics();

% m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,[]);
% m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
% m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
