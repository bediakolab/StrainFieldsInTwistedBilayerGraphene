% simulated_data_driver_revisedVM.m
%
% This is not the script for beamwidth averaging.
% This script designed to make plots of the revised vonMises strain.
%
% Nathanael Kazmierczak, 05/29/2020, Bediako Lab, UC Berkeley.

%% Dataset construction
moire_angle = 1.1;
recon_angle = 1.0;
recon_distance = 45;
scan_dimensions = [200,200];
scan_step = 0.5;
scan_angle = 55;
add_noise_amount = 0;
plot_sample = true;
raster_start_pos = [5,7];


[ displacement_field ] = makeSimulatedDisplacementDataset(moire_angle, recon_angle, recon_distance, scan_dimensions,...
    scan_step, scan_angle, add_noise_amount, plot_sample, raster_start_pos);

%% Geometry registration and annealing.
m4 = FourDSTEM_Analysis_Engine;
m4.DSC_fit_storage = displacement_field;
m4.datacube_size = [0,0,scan_dimensions];
m4.scan_stepsize = scan_step;
m4.reformAxes();


filteramp = 10;
filterrange = 3;
hotspot_suppression = 0;
edge_threshold = 10;
gaussian_filter_threshold = 0.2;
AAblur_Gaussian_sigma = 7;
use_robust_multistart = 1;
AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold,AAblur_Gaussian_sigma);
m4.plotAACircleFitOnCustomColor([],[],0);

% [pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1);
make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics(-90:10:90);

optional_SP_tolerance = 0.5;  % This has been the default.
optional_AAcircle_min_radius = 4;  % in pixels.
AA_inflation_factor = 1.5;
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance,optional_AAcircle_min_radius);


% 04/28/2020: additionally get information about AB domains
% setup

boundary_inflation_margin = 10;  % For getting AB area, don't want this.
boundary_trim_margin = 20;  % For doing the small stepsizes
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
% actual statistics
% disallow_edge_intersect = true;
% [pixel_sizes,nm_sizes] = m4.getABareaStatistics(disallow_edge_intersect,boundary_trim_margin);

%% Anneal dataset 
m4.setEmitterOrientations();
m4.annealDisplacementField2(filteramp,filterrange,1);

method = 'outside input';
AA_type_string = [];
window_size = [];
window_increment = [];
crop_val_pixels = 40;
angle_input = 1.085;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

%% Map Strain (Can just run this part upon loading saved dataset
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/SimulatedDatasets';
m4.saved_plots_foldername = '05292020_1p085_recon1p0_45_200x200_0p5nmstep___VMrevisednofilter';

load('demo_colormap.mat');
filteramp = 10;
filterrange = 3;
saveplotflag = 1;
% m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
% filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
filterstruct = buildFilterStructPreset('blank filter');
% stringcell = {'median outlier','xdisp',[31,0.2];'median outlier','ydisp',[31,0.2];'median normal','xdisp',[15];'median normal','ydisp',[15];'moving average circle','xdisp',31;'moving average circle','ydisp',31};
% filterstruct = buildFilterStruct(stringcell);
sign_convention = [+1,-1];
trim_edge_pixel_num = [40,0,0];
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);
m4.saveObject();
%             