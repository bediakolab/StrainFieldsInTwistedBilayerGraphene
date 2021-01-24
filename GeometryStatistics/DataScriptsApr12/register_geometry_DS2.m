% register_geometry_DS2.m

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset2');
load('DS2objectdata.mat');
m4.scan_stepsize = 0.5;
m4.reformAxes();

filteramp = 0.3;
filterrange = 5;
hotspot_suppression = 0;
% We will use the new filter here in a limited capacity. First try with old
% filter.


AA_inflation_factor = 1.2;
edge_threshold = 5;
gaussian_filter_threshold = 0.1;
use_robust_multistart = 1;


AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold);
m4.plotAACircleFitOnCustomColor(0.3,7,0)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)

% make_plots = true; % for debugging only
make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics(-90:10:90);

optional_SP_tolerance = 0.3;  % This has been the default.
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance);


% 04/28/2020: additionally get information about AB domains
% setup
boundary_inflation_margin = 0;  % For getting AB area, don't want this.
boundary_trim_margin = 15;
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
% actual statistics
disallow_edge_intersect = true;
[pixel_sizes,nm_sizes] = m4.getABareaStatistics(disallow_edge_intersect,boundary_trim_margin);
% figure; histogram(m4.ABareas_nm(:));
% title('AB skeleton domain areas');
% xlabel('Area (nm^2)');
% ylabel('Counts');
% figure; histogram(m4.ABsidelengths(:));
% title('AB skeleton domain sidelengths');
% xlabel('Nearest-neighbor AA straight line distance (nm)');
% ylabel('Counts');
% figure; histogram(m4.ABangles(:,1));
% title('AB skeleton domain smallest angle');
% xlabel('Smallest domain angle (degrees)');
% ylabel('Counts');
% figure; histogram(m4.ABangles(:,3)-m4.ABangles(:,1));
% title('AB skeleton domain triangular distortion');
% xlabel('Smallest domain angle minus largest domain angle (degrees)');
% ylabel('Counts');
[mean_area,stdev_area,mean_sidelength,stdev_sidelength,...
    mean_smallest_angle,stdev_smallest_angle,mean_angle_difference,...
    stdev_angle_difference] = m4.plotABareaStatistics()

%% Add a triangulation and prediction of approximate twist angles. 
m4.triangulateAA(2);  % Using Gaussian ellipse;
m4.setFacetTwistAngles();
[average_angle,std_angle,max_angle,min_angle] = m4.plotFacetedTwistAngles([0.9,inf]);std_angl

