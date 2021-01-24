% register_geometry_DS23S2.m
%
% Geometry registration script for the
% Session 2 data. 
%
% Nathanael Kazmierczak, 06/22/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/Session2/DS23S2');
load('DS23S2objectdata.mat');
m4.scan_stepsize = 0.5;  % correct for tear datasets as well.
m4.reformAxes();
m4.DSC_fit_storage(:,:,1) = flipud(m4.DSC_fit_storage(:,:,1));
m4.DSC_fit_storage(:,:,2) = flipud(m4.DSC_fit_storage(:,:,2));
m4.multistart_displacement_convergence = flipud(m4.multistart_displacement_convergence);
m4.classified_displacement_convergence = flipud(m4.classified_displacement_convergence);

filteramp = 0.3;
filterrange = 5;
hotspot_suppression = 0;
% We will use the new filter here in a limited capacity. First try with old
% filter.


AA_inflation_factor = 1.5;
% AA_inflation_factor = 2;
% edge_threshold = 5;
edge_threshold = -1;
gaussian_filter_threshold = 0.1;
% gaussian_filter_threshold = 0.065;
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
boundary_trim_margin = 3;
tear_boundary_trim_margin = 2;
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin,tear_boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
yn = input('Do you want to make manual adjustments to the saddle point registration?');
if yn
    m4.adjustSolitonWalls();
end

%% Add a triangulation and prediction of approximate twist angles. 
m4.triangulateAA(2);  % Using circle fit because of the poor boundary fits on DS8S2
m4.setFacetTwistAngles();
m4.repairTriangulation();
[average_angle,std_angle,max_angle,min_angle,ntriangles] = m4.plotFacetedTwistAngles([0,inf]);

