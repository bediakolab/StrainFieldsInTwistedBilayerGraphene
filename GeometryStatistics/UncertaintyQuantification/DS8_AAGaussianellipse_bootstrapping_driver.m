% DS8_AAGaussianellipse_bootstrapping_driver.m
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

% converge_storage = m4.fitAAtoGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma);
% Delete ellipse 12, which really went off the rails
% m4.AAgaussianfit_FWHM_ellipses(12,:) = nan;
% m4.AAgaussianfit_raw_parameters(12,:) = nan;
m4.plotAAGaussianEllipseFitOnCustomColor(filteramp,filterrange,hotspot_suppression);
[semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
    m4.getAAellipticGaussianFitStatistics(-90:10:90);

bootstrap_idx = 61;
nboot = 1000;
[stderrs_of_fit,CI_lb,CI_ub,raw_bootstrap_data] = m4.bootstrapAAGaussianEllipse(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma,bootstrap_idx,nboot);
% convert Std errors and confidence intervals out here for consistency with
% the last function.
timevec = [1,0.5,0.5,0.5,0.5];
minusvec = [0,0,0,1,1];
boot_stderrs = stderrs_of_fit .* timevec
boot_CIlb = (CI_lb-minusvec) .* timevec
boot_CIub = (CI_ub-minusvec) .* timevec
scaled_bootdata = raw_bootstrap_data .* timevec;

% Make histograms
figure
ellipseangles = rad2deg(scaled_bootdata(:,1));
ellipseangles(ellipseangles<0) = ellipseangles(ellipseangles<0) + 180;
histogram(ellipseangles);
title('Ellipse angle distribution');
xlabel('Angle of rotation (deg)');
ylabel('Counts');

figure
histogram(scaled_bootdata(:,2));
title('Semimajor bootstrap distribution');
xlabel('Axis length (nm)');
ylabel('Counts');

figure
histogram(scaled_bootdata(:,3));
title('Semiminor bootstrap distribution');
xlabel('Axis length (nm)');
ylabel('Counts');

figure
scatter(scaled_bootdata(:,2),scaled_bootdata(:,3));
title('Semimajor and semiminor scatterplot');
xlabel('Semimajor axis length (nm)');
ylabel('Semiminor axis length (nm)');

figure
scatter(scaled_bootdata(:,4),scaled_bootdata(:,5));
title('Ellipse center');
xlabel('x location (nm)');
ylabel('y location (nm)');

% Actual vartest here
% Get AA test on semimajor
% Perform Chi-square test
varvec = m4.AAgaussianfit_FWHM_ellipses(:,2)*m4.scan_stepsize;
fitvar = 0.00649636; % Averaged radius variance over three bootstraps -- see 06262020AAbootstrappingworkbook.xlsx.
[h,p,ci,stats] = vartest(varvec,fitvar);



