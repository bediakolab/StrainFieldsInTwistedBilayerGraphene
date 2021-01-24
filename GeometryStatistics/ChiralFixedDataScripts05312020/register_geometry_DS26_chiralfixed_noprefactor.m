% register_geometry_DS26_Apr12.m

load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/FirstTryFullData/04032020_ds26_fullannealing_m4only.mat')
m4.scan_stepsize = 0.5;
m4.reformAxes();

m4.DSC_fit_storage(:,:,1) = flipud(m4.DSC_fit_storage(:,:,1));
m4.DSC_fit_storage(:,:,2) = flipud(m4.DSC_fit_storage(:,:,2));
m4.multistart_displacement_convergence = flipud(m4.multistart_displacement_convergence);
m4.classified_displacement_convergence = flipud(m4.classified_displacement_convergence);

filteramp = 0.3;
filterrange = 9;
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
    m4.getAAellipticGaussianFitStatistics()

boundary_inflation_margin = 0;  % For getting AB area, don't want this.
boundary_trim_margin = 15;
m4.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor,[]);
m4.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
m4.plotGeometryFit(filteramp,filterrange,hotspot_suppression);

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

m4.triangulateAA(2,true);  % Using Gaussian ellipse, manual triangulation.;
m4.setFacetTwistAngles();
[average_angle,std_angle,max_angle,min_angle,n_triangles_used] = m4.plotFacetedTwistAngles([0,inf]);
