% revised_script_for_data_fitting.m
%
% Script for using some of the revised fitting functions developed by NPK.
%
% Goal #1: use this to fit the day1, dataset9 (dataset which has the tear
% and two different domains flanking on either side).
%
% Nathanael Kazmierczak

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
%m4.setAveragedDP(1);
m4.plotAveragedDP(1);
m4.makeMasksForElasticFit();
m4.fitElasticScattering('lorentzian');
m4.interpolateThroughPeakResiduals();
% m4.getLineAverageElasticResiduals();
m4.plotElasticMask(0);
m4.plotElasticMask(1);
subtract_elastic_flag = 1;
use_average_flag = 1;
radial_elastic_correction_flag = 0;
interpolation_correction_flag = 1;
exponent = 1;
removeHBN = 1;
removeGraphene = 1;
removeBeamstop = 1;
m4.setAveragedDP(2);
m4.plotBaselinedData(exponent,removeHBN,removeGraphene,removeBeamstop)
m4.saveObject();

saveplotflag = 1;
shading_type = 'flat';
normalize = 0;
m4.setBlinkingDisks();
m4.integrateDisks(subtract_elastic_flag,use_average_flag,radial_elastic_correction_flag,interpolation_correction_flag);  % omit loadnumcroprange argument for default, loading all.
m4.makeBlinkingPlots([],saveplotflag,shading_type,normalize);
m4.saveObject();

m4.setTrigHeight();
% fit_boundaries = {[1,80],[1,80]};
useParallel = 1;
m4.fitBlinking3([],useParallel,[]);
m4.saveObject();

individual_residual_flag = 1;
median_filter_threshold = [];
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residual_flag,shading_type);
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.saveObject();

m4.makeInteractiveDisplacementDarkFieldPlot();
m4.makeInteractiveInverseDarkFieldPlot();
