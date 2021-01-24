% DS4S2_GaussianCircleFits.m
%
% First in a series of new scripts for getting the AA geometry for ONR
% presentation and paper.
%
% Nathanael Kazmierczak, 07/05/2020

% addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/Session2');
load('06222020_DS4S2_prepared_for_annealing.mat');

filteramp = 0.3;
filterrange = 5;
hotspot_suppression = false;
edge_threshold = 5;  % Note: the original setting was -1 for geometry registration of this dataset.
% However, for AA radius quantification, we want to make sure we are far
% enough away from the edge to get accurate information.
% Setting edge_threshold = 5 controls this.
gauss_filter_threshold = 0.1;
AA_amp_threshold = 0.71; % Angstrom
make_plots = false;  % should validate these at some point
mask_radius = 10;
AAblur_Gaussian_sigma = 5;  % the default


m4.fitAAtoGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma);
% statsrange = [1,inf];  % to get rid of that one bad outlier circle.
[radius_mean_nm,radius_std_nm,radius_stderr_nm,radius_95confdelta_nm] = m4.getGaussianCircleFitStatistics()

AA_type = 3;  % this will give us Gaussian circle fit for the AA domain centers.
m4.triangle_sidelengths = [];
m4.triangle_moire_angle = [];
m4.triangulateAA(AA_type,false);
m4.setFacetTwistAngles();
m4.repairTriangulation();
[average_angle,std_angle,max_angle,min_angle,n_triangles_used,angles_used] = ...
    m4.plotFacetedTwistAngles([0,inf],pink,false,false,true);

