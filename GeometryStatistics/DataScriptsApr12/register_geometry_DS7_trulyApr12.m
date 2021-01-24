% register_geometry_DS7_trulyApr12.m
%
% Nathanael Kazmierczak, 05/05/2020

load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/DS7/DS7_newfit_objectdata.mat');
m4.scan_stepsize = 0.5;
m4.reformAxes();
% m4.classifyDisplacementConvergence();

filteramp = 0.3;
filterrange = 5;
% filterrange = 9;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0);
hotspot_suppression = 0;
% AA_inflation_factor = 2;
AA_inflation_factor = 5;
edge_threshold = 5;
gaussian_filter_threshold = 0.1;
use_robust_multistart = 1;
AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold);

m4.plotAACircleFitOnCustomColor(filteramp,filterrange,0)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1,[1:0.2:2.5])%,[1.2:0.2:2.4])

make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics()

optional_SP_tolerance = 0.55;  % saturation has to be larger than this for saddle point registration. 0.3 default or leave blank
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance);
boundary_inflation_margin = 15;
boundary_trim_margin = 15;
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);

% AB domain information (04/28/2020)
disallow_edge_intersect = true;
[pixel_sizes,nm_sizes] = m4.getABareaStatistics(disallow_edge_intersect,boundary_trim_margin);
[mean_area,stdev_area,mean_sidelength,stdev_sidelength,...
    mean_smallest_angle,stdev_smallest_angle,mean_angle_difference,...
    stdev_angle_difference] = m4.plotABareaStatistics()
