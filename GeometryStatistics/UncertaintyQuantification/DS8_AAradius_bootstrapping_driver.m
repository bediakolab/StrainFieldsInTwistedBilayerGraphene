% DS8_AAradius_bootstrapping_driver.m
%
% Driver script for bootstrapping the AA circle fit and performing
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
bidx = 53;
nboot = 100;
[boot_95CIs,boot_stderrs,bootstrap_bestcoords,bootstrap_bestresiduals] = m4.getAAFitUncertainty(gauss_filter_threshold,use_robust_multistart,AA_amp_threshold,AAblur_Gaussian_sigma,filteramp,filterrange,bidx,nboot);

% convert to nm for analysis
bootstrap_bestcoords = bootstrap_bestcoords*m4.scan_stepsize;
boot_95CIs = boot_95CIs*m4.scan_stepsize
boot_stderrs = boot_stderrs*m4.scan_stepsize

figure; histogram(bootstrap_bestcoords(:,3));
xlabel('Radius (nm)');
ylabel('Counts');
title('Bootstrap distribution of AA circle radii');
figure; scatter(bootstrap_bestcoords(:,1),bootstrap_bestcoords(:,2))
xlabel('x position (nm)');
ylabel('y position (nm)');
title('Bootstrap distribution of AA circle centers');

% Get AA circlefit radii and variance
nmradii = m4.AA_circlefit_values(:,3)*m4.scan_stepsize;
varnmradii1 = std(nmradii)^2;
varnmradii2 = var(nmradii);

% Perform Chi-square test
varvec = m4.AA_circlefit_values(:,3)*m4.scan_stepsize;
fitvar = 0.00976; % Averaged radius variance over three bootstraps -- see 06262020AAbootstrappingworkbook.xlsx.
[h,p,ci,stats] = vartest(varvec,fitvar);


