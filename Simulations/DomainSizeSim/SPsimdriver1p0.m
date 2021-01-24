% SPsimdriver1p0.m
%
% This is not the script for beamwidth averaging.
%
% Nathanael Kazmierczak, 05/29/2020, Bediako Lab, UC Berkeley.

%% Dataset construction
moire_angle = 1.0;
recon_struct.type = 'Tadmor';
recon_struct.recon_angle = 0;
recon_struct.recon_distance = 0;
heterostrain_struct = [];
% scan_dimensions = [100,100];
scan_dimensions = [50,50];
scan_step = 0.5;
% scan_angle = 0;
scan_angle = 10;
add_noise_amount = 0;
plot_sample = true;
% raster_start_pos = [5,7];
raster_start_pos = [4,2];
customCellDimXY = [];


[ displacement_field ] = makeSimulatedDisplacementDataset(moire_angle, recon_struct, heterostrain_struct, scan_dimensions,...
    scan_step, scan_angle, add_noise_amount, plot_sample, raster_start_pos, customCellDimXY);


%% Geometry registration and annealing.
m4 = FourDSTEM_Analysis_Engine;
m4.DSC_fit_storage = displacement_field;
m4.datacube_size = [0,0,scan_dimensions];
m4.scan_stepsize = scan_step;
m4.reformAxes();


filteramp = 10;
filterrange = 3;
hotspot_suppression = 0;
edge_threshold = 0;
% gaussian_filter_threshold = 0.2;
% AAblur_Gaussian_sigma = 7;
gaussian_filter_threshold = 0.2;
AAblur_Gaussian_sigma = 7;

% use_robust_multistart = 1;
AA_amp_threshold = 0.71; % Angstrom
% m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold,AAblur_Gaussian_sigma);
% m4.plotAACircleFitOnCustomColor([],[],0);

% [pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1);
make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma)
[radius_mean_nm,radius_std_nm,radius_stderr,radius_95confdelta] = m4.getGaussianCircleFitStatistics([0,inf])

% 
optional_SP_tolerance = 0.5;  % This has been the default.
optional_AAcircle_min_radius = 6;  % in pixels.
AA_inflation_factor = 1;
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance,optional_AAcircle_min_radius);

% 
% % 04/28/2020: additionally get information about AB domains
% % setup
% 
boundary_inflation_margin = 0;  % For getting AB area, don't want this.
boundary_trim_margin = 10;  % For doing the small stepsizes
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
% actual statistics
% disallow_edge_intersect = true;
% [pixel_sizes,nm_sizes] = m4.getABareaStatistics(disallow_edge_intersect,boundary_trim_margin);
% 
% %% Anneal dataset 
m4.setEmitterOrientations('v1',true);
m4.annealDisplacementField2(filteramp,filterrange,1);

%% Do the SP geometry registration, which was the original point of this.
filterstruct = [];
cut_number = 10;  % along the base_t direction
cut_spacing_nm = 0.1;
AA_center_type = 'Gaussian circle AA fit';
radius_value = 3;
image_boundary_buffer = 10;  % pixels
pixel_size_cutoff = 4;
window_dims_nm = [3,5];  % First value is for the width, second for the length
cutoff_score = 0.5;
[distcoords,av_scores,distance,endpoint1,endpoint2,rangemean,rangemeanstderr] = ...
    m4.getSPwidths(filterstruct,cut_number,cut_spacing_nm,AA_center_type,radius_value,image_boundary_buffer,pixel_size_cutoff,window_dims_nm,cutoff_score)

% 
% method = 'outside input';
% AA_type_string = [];
% window_size = [];
% window_increment = [];
% crop_val_pixels = 30;
% angle_input = 1.085;
% m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);
% 
% %% Map Strain (Can just run this part upon loading saved dataset
% m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/SimulatedDatasets';
% m4.saved_plots_foldername = '05292020_1p085_recon1p0_45_500x500_0p1nmstep___Customfilteredplots_heavy';
% 
% load('demo_colormap.mat');
% filteramp = 10;
% filterrange = 3;
% saveplotflag = 1;
% % m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
% % filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
% % filterstruct = buildFilterStructPreset('blank filter');
% stringcell = {'median outlier','xdisp',[31,0.2];'median outlier','ydisp',[31,0.2];'median normal','xdisp',[15];'median normal','ydisp',[15];'moving average circle','xdisp',31;'moving average circle','ydisp',31};
% filterstruct = buildFilterStruct(stringcell);
% sign_convention = [-1,+1];
% trim_edge_pixel_num = [80,0,0];
% [eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
%                 m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);
% m4.saveObject();
%             