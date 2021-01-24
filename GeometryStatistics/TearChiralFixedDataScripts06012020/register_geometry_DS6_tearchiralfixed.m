% register_geometry_DS6_tearchiralfixed.m
%
% The tear datasets also have the counterclockwise red/purple/orange on
% both sides of the tear, so they need to be flipped in the same way that
% the other datasets are flipped.
%
% Need to build a tear mask to mask off the tear region, and probably
% erosion filter immunity around the edges so we can get a complete
% geometry registration.
%
% Nathanael Kazmierczak, 06/01/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/Dataset6');
load('DS6objectdata.mat');
m4.scan_stepsize = 0.5;  % correct for tear datasets as well.
m4.reformAxes();
m4.DSC_fit_storage(:,:,1) = flipud(m4.DSC_fit_storage(:,:,1));
m4.DSC_fit_storage(:,:,2) = flipud(m4.DSC_fit_storage(:,:,2));
m4.multistart_displacement_convergence = flipud(m4.multistart_displacement_convergence);
m4.classified_displacement_convergence = flipud(m4.classified_displacement_convergence);

m4.makeTearRegionDisplacementMask();

filteramp = 0.3;
filterrange = 9;
hotspot_suppression = 0;
% We will use the new filter here in a limited capacity. First try with old
% filter.


AA_inflation_factor = 1.5;
% AA_inflation_factor = 2;
% edge_threshold = 5;
edge_threshold = -1;
% gaussian_filter_threshold = 0.1;
gaussian_filter_threshold = 0.065;
use_robust_multistart = 1;


AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold);
m4.plotAACircleFitOnCustomColor(filteramp,filterrange,0)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)

% make_plots = true; % for debugging only
make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics(-90:10:90);

% optional_SP_tolerance = 0.45;  % This has been the default.
% m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance);

optional_SP_tolerance = 0.4;  % This has been the default.
optional_AAcircle_min_radius = 0;  % in pixels.
include_untrimmed = false;
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance,optional_AAcircle_min_radius,include_untrimmed);

% 04/28/2020: additionally get information about AB domains
% setup
boundary_inflation_margin = 0;  % For getting AB area, don't want this.
% boundary_trim_margin = 15;
boundary_trim_margin = 6;
tear_boundary_trim_margin = 2;
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin,tear_boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
yn = input('Do you want to make manual adjustments to the saddle point registration?');
if yn
    m4.adjustSolitonWalls();
end

% actual statistics
% % % disallow_edge_intersect = true;
% % % [pixel_sizes,nm_sizes] = m4.getABareaStatistics(disallow_edge_intersect,boundary_trim_margin);


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


% Not sure that we need these anymore, at least while testing.
% [mean_area,stdev_area,mean_sidelength,stdev_sidelength,...
%     mean_smallest_angle,stdev_smallest_angle,mean_angle_difference,...
%     stdev_angle_difference] = m4.plotABareaStatistics()

%% Add a triangulation and prediction of approximate twist angles. 
m4.triangulateAA(2);  % Using Gaussian ellipse;
m4.setFacetTwistAngles();
[average_angle,std_angle,max_angle,min_angle] = m4.plotFacetedTwistAngles([0.5,inf]);

