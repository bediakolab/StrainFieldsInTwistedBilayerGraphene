% DS8_AAGaussianCircle_bootstrapping_driver.m
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

bootstrap_idx = 61;
nboot = 1000;

% m4.fitAAtoGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma);
% Delete ellipse 12, which really went off the rails
[stderrs_of_fit,CI_lb,CI_ub,raw_bootstrap_data] = m4.bootstrapAAGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma,bootstrap_idx,nboot);
boot_stderrs = stderrs_of_fit*m4.scan_stepsize
boot_CI_lb = (CI_lb-[1,1,0])*m4.scan_stepsize
boot_CI_ub = (CI_ub-[1,1,0])*m4.scan_stepsize
bootdata = (raw_bootstrap_data-[1,1,0])*m4.scan_stepsize;

figure
histogram(bootdata(:,3));
xlabel('Radius (nm)');
ylabel('Counts');
title('Gaussian circle fit radius bootstrap distribution');

figure
scatter(bootdata(:,1),bootdata(:,2));
xlabel('x position (nm)');
ylabel('y position (nm)');
title('Gaussian circle fit position distribution');


