% DS8_AAGaussianCircle_testing_driver.m
%
% Driver script for bootstrapping the AA Gaussian ellipse fit and performing
% hypothesis tests.
%
% Nathanael Kazmierczak, 06/26/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
    load('05312020_chiralfixed_DS8_postannealing.mat');
end

gauss_filter_threshold = 0.1;
use_robust_multistart = 1;
AA_amp_threshold = 0.71;
AAblur_Gaussian_sigma = 5;
filteramp = 0.3;
filterrange = 5;
bidx = [];

hotspot_suppression = false;
edge_threshold = 5;
make_plots = false;
mask_radius = 10;  % 25 is way too large and includes all kinds of other AA domains as well

m4.fitAAtoGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma);
% Delete ellipse 12, which really went off the rails

[radius_mean_nm,radius_std_nm] = m4.getGaussianCircleFitStatistics();


