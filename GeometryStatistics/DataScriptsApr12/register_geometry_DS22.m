% register_geometry_DS22.m

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/Dataset22');
    load('DS22objectdata.mat');
end
m4.scan_stepsize = 0.2;
m4.reformAxes();

filteramp = 0.3;
filterrange = 11;
hotspot_suppression = 0;
% We will use the new filter here in a limited capacity. First try with old
% filter.
% stringcell = {'soft multistart','amplitude',[0.3,5]};%,'hard multistart','amplitude',[0.3,5]};
% filterstruct = buildFilterStruct(stringcell);
% m4.makeCustomDisplacementColorPlot([],[],filterstruct,[],1,0,0);
useImmunity = false;
if useImmunity
    if isempty(m4.filter_immunity_mask)
        m4.drawFilterImmunityRegions();
    end
end
if hotspot_suppression
    if isempty(m4.hotspot_suppression_mask)
        m4.makeCustomDisplacementColorPlot();  % this can't be used for drawing the roipoly because of the cursor glitch, but keep it as a reference.
        m4.setManualHotspotSuppressionMask();
    end
end

% string_cell = {'hard multistart','amplitude',[0.3,11],'soft multistart','amplitude',[0.3,11],'soft multistart','amplitude',[0.3,7]};
% string_cell = {'hard multistart','amplitude',[0.3,11];'median outlier','amplitude',[3,0.3];'median outlier','amplitude',[3,0.3]};
% [ filterstruct ] = buildFilterStruct( string_cell );
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,hotspot_suppression);
% m4.makeCustomDisplacementColorPlot([],[],filterstruct,[],1,0,0);

AA_inflation_factor = 2.8;
edge_threshold = 0;
% gaussian_filter_threshold = 0.1;
gaussian_filter_threshold = 0.4;
AAblur_Gaussian_sigma = 7;
use_robust_multistart = 1;


AA_amp_threshold = 0.71; % Angstrom
m4.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,use_robust_multistart,AA_amp_threshold,AAblur_Gaussian_sigma);
m4.plotAACircleFitOnCustomColor(filteramp,filterrange,0)
[pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)

% make_plots = true; % for debugging only
make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well
m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma)
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics(-90:10:90);

optional_SP_tolerance = 0.7;  % This has been the default.
optional_AAcircle_min_radius = 7;  % in pixels.
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,optional_SP_tolerance,optional_AAcircle_min_radius);

%% Manual modification of SP wall here.
%%

% 04/28/2020: additionally get information about AB domains
% setup
boundary_inflation_margin = 0;  % For getting AB area, don't want this.
boundary_trim_margin = 10;
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
m4.triangulateAA(0);  % Using Gaussian ellipse;
m4.setFacetTwistAngles();
% [average_angle,std_angle,max_angle,min_angle,n_triangles_used] = m4.plotFacetedTwistAngles([0,inf])

% %% DS5 has an outlier AA region
% % Delete it in both the circlefit mask and Gaussian ellipse parameter lists
% % before computing reliable statistics.
% [val,idx] = min(m4.AA_circlefit_values(:,3));
% cvo = m4.AA_circlefit_values(idx,3);
% evo = m4.AAgaussianfit_FWHM_ellipses(idx,1:3);
% m4.AA_circlefit_values(idx,3) = nan;
% m4.AAgaussianfit_FWHM_ellipses(idx,1:3) = nan;
% % Now the statistics will be accurate
% [pixel_radii,nm_radii,nm_average,nm_std] = m4.getAACirclefitStatistics(1)
% [semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
%     m4.getAAellipticGaussianFitStatistics(-90:10:90)
% m4.AA_circlefit_values(idx,3) = cvo;
% m4.AAgaussianfit_FWHM_ellipses(idx,1:3) = evo;

%% See if we can get anything here
m4.setEmitterOrientations();
m4.annealDisplacementField2(filteramp,filterrange,1);
method = 'outside input';
AA_type_string = [];
window_size = [];
window_increment = [];
% window_size = 80;
% window_increment = 10;
crop_val_pixels = 30;
angle_input = 0.8453;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

load('demo_colormap.mat');
saveplotflag = 0;
% filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
% stringcell = {'median outlier','xdisp',[5,2]};
% stringcell = {'median outlier','xdisp',[5,0.2];'median outlier','ydisp',[5,0.2];'median normal','xdisp',3;'median normal','ydisp',3;'median normal','xdisp',3;'median normal','ydisp',3};

stringcell = {'median outlier','xdisp',[21,0.7];'median outlier','ydisp',[21,0.7];'median outlier','xdisp',[7,0.2];'median outlier','ydisp',[7,0.2];'median normal','xdisp',7;'median normal','ydisp',7;'median normal','xdisp',11;'median normal','ydisp',11;'moving average circle','xdisp',11;'moving average circle','ydisp',11};
% This is an interesting attempt without the moving average circle, but it
% is *so* messy. Still no sign of the mythic unstrained core, however.
% stringcell = {'median outlier','xdisp',[21,0.7];'median outlier','ydisp',[21,0.7];'median outlier','xdisp',[7,0.2];'median outlier','ydisp',[7,0.2];'median normal','xdisp',7;'median normal','ydisp',7;'median normal','xdisp',11;'median normal','ydisp',11};
filterstruct = buildFilterStruct(stringcell);
sign_convention = [-1,-1];
trim_edge_pixel_num = [0,0,0];
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);

figHandles = get(groot,'Children');
for i = 1:numel(figHandles)
    figure(figHandles(i));
    xlim([2,10]);
    ylim([0,8]);
    caxis([-5,5]);
end
