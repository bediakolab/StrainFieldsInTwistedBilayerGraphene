% DS1S2_GaussianCircleFits.m
%
% First in a series of new scripts for getting the AA geometry for ONR
% presentation and paper.
%
% Nathanael Kazmierczak, 07/05/2020


addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/Session2/DS1S2');
load('DS1S2objectdata.mat');

m4.scan_stepsize = 0.5;   % appears to have been left out in original script.
filteramp = 0.3;
filterrange = 3;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,0);


hotspot_suppression = false;
edge_threshold = 5;
gauss_filter_threshold = 0.1;
AA_amp_threshold = 0.71; % Angstrom
make_plots = false;  % should validate these at some point
mask_radius = 10;
AAblur_Gaussian_sigma = 5;  % the default


m4.fitAAtoGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,AA_amp_threshold,AAblur_Gaussian_sigma);
% statsrange = [1,inf];  % to get rid of that one bad outlier circle.
[radius_mean_nm,radius_std_nm,radius_stderr_nm,radius_95confdelta_nm] = m4.getGaussianCircleFitStatistics()


% Perform triangulation for the first time here.
m4.triangulateAA(3);  % Using Gaussian circle fit
m4.setFacetTwistAngles();
m4.repairTriangulation();
[average_angle,std_angle,max_angle,min_angle,ntriangles] = m4.plotFacetedTwistAngles([0,inf]);


