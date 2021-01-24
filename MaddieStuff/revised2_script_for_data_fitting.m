% revised2_script_for_data_fitting.m
%
% Nathanael Kazmierczak, 04/11/2020

% For Maddie: It is best if you rename the files to different filenames so
% there is no danger of accidentally loading the wrong file. 
filename = 'MV_02262020_10.h5';
probe_filename = []; % Replace this once Maddie uploads it. '60kV_1mrad_spot9_0p5s_130CL_defocus = 0.dm3';
currentdir = pwd;
cd ..
addpath(genpath(pwd));
cd(currentdir);
addpath('D:\Maddie Local 2\4D STEM\20200226\10-200x197-0p5nm large unit cell');
results_path = 'D:\Maddie Local 2\4D STEM\20200226 Results\10-200x197-0p5nm large unit cell Results';
% This is the manual way of setting where the pictures will be saved to. If
% you don't want automated saving of pictures, it really doesn't matter!
disks_name = [];
stepsize = 1; 
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);
m4.saved_plots_folderpath = 'D:\Maddie Local 2\4D STEM\20200226 Results';
m4.saved_plots_foldername = '10-200x197-0p5nm large unit cell Results';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m4.setAveragedDP('all');
m4.setAveragedDP('group1');
m4.plotAveragedDP(1);
m4.makeMasksForElasticFit();
m4.plotElasticMask(0);
m4.plotElasticMask(1);
m4.setBlinkingDisks();
m4.fitElasticScattering('ringed gaussian');
% Note for Maddie: it's good to examine the plot with the four panels to
% ensure that the ringed gaussian function is working well. If it is not,
% then you can use 'lorentzian' instead and proceed with the rest of the
% script normally.
m4.getCOMbeamCenter();
m4.setLineAverageElasticResiduals('com');
m4.interpolateThroughPeakResiduals(1);
% m4.getLineAverageElasticResiduals();

subtract_elastic_flag = 1;
use_average_flag = 1;
exponent = 1;
removeHBN = 1;
removeGraphene = 1;
removeBeamstop = 1;
m4.setAveragedDP('group2');
subtractScattering = 1;
subtractRadialAverage = 1;
subtractCartesianInterpolation = 0;
subtractPolarInterpolation = 1;
% Note for Maddie: This plot (in grayscale) should show very little
% systematic signal, probably just a bit of noise around max amplitude of
% 0.05. If it is much larger, something is wrong.
m4.plotBaselinedData(exponent,removeHBN,removeGraphene,removeBeamstop,...
     subtractScattering,subtractRadialAverage,subtractCartesianInterpolation,subtractPolarInterpolation);
m4.saveObject();

saveplotflag = 1;
shading_type = 'flat';
normalize = 0;

% m4.integrateDisks(subtract_elastic_flag,use_average_flag,radial_elastic_correction_flag,interpolation_correction_flag);  % omit loadnumcroprange argument for default, loading all.
% loadnum_crop_range = 1:8;
m4.integrateDisks(use_average_flag,subtractScattering,subtractRadialAverage,subtractCartesianInterpolation,...
                subtractPolarInterpolation);
m4.makeBlinkingPlots([],saveplotflag,shading_type,normalize);
m4.saveObject();

m4.setTrigHeight();
fit_boundaries = [];
useParallel = true;
useAnalyticJacobian = true;
useFittedPrefactors = false;
m4.fitBlinking3([],useParallel,fit_boundaries,useAnalyticJacobian,useFittedPrefactors);
m4.saveObject();
m4.makeCustomDisplacementColorPlot([],[],[],[],1,0);
title('Unfiltered displacement after first fit.');

m4.fitBlinking4([],useParallel);
m4.makeCustomDisplacementColorPlot([],[],[],[],1,1);
title('Unfiltered displacement after second fit.');
useFittedPrefactors = true;
m4.fitBlinking3([],useParallel,fit_boundaries,useAnalyticJacobian,useFittedPrefactors);
m4.makeCustomDisplacementColorPlot([],[],[],[],1,0);
title('Unfiltered displacement after third fit.');
m4.classifyDisplacementConvergence();
m4.makeCustomDisplacementColorPlot([],[],0.4,7,1,0);
title('Filtered displacement after third fit.');

individual_residual_flag = 1;
median_filter_threshold = [];
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residual_flag,shading_type);
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.saveObject();

m4.makeInteractiveDisplacementDarkFieldPlot();
m4.makeInteractiveInverseDarkFieldPlot();
