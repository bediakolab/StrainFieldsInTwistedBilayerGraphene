% dataset10_polar_interpolation.m
%
% Nathanael Kazmierczak, 04/11/2020

filename = 'TEAM1Day1Dataset10.h5';
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset10');
probe_filename = [];
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/BackgroundFitTesting/Dataset10PolarInterpolation';
% This is the manual way of setting where the pictures will be saved to. If
% you don't want automated saving of pictures, it really doesn't matter!
disks_name = [];
stepsize = 0.5; 
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m4.setAveragedDP('all');
m4.setAveragedDP('group1');
m4.plotAveragedDP(1);
m4.makeMasksForElasticFit();
m4.plotElasticMask(0);
m4.plotElasticMask(1);
m4.setBlinkingDisks();
m4.fitElasticScattering('ringed gaussian');
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
m4.plotBaselinedData(exponent,removeHBN,removeGraphene,removeBeamstop,...
     subtractScattering,subtractRadialAverage,subtractCartesianInterpolation,subtractPolarInterpolation);
m4.saveObject();

saveplotflag = 1;
shading_type = 'flat';
normalize = 0;

% m4.integrateDisks(subtract_elastic_flag,use_average_flag,radial_elastic_correction_flag,interpolation_correction_flag);  % omit loadnumcroprange argument for default, loading all.
loadnum_crop_range = 1:8;
m4.integrateDisks(use_average_flag,subtractScattering,subtractRadialAverage,subtractCartesianInterpolation,...
                subtractPolarInterpolation,loadnum_crop_range);
m4.makeBlinkingPlots([],saveplotflag,shading_type,normalize);
m4.saveObject();

m4.setTrigHeight();
fit_boundaries = {[1,35],[1,35]};
useParallel = false;
useAnalyticJacobian = true;
useFittedPrefactors = false;
m4.fitBlinking3([],useParallel,fit_boundaries,useAnalyticJacobian,useFittedPrefactors);
m4.saveObject();
m4.makeCustomDisplacementColorPlot([],[],[],[],1,0);

m4.fitBlinking4([],useParallel,[35,35]);
m4.makeCustomDisplacementColorPlot([],[],[],[],1,1);
useFittedPrefactors = true;
m4.fitBlinking3([],useParallel,fit_boundaries,useAnalyticJacobian,useFittedPrefactors);
m4.makeCustomDisplacementColorPlot([],[],[],[],1,0);
m4.classifyDisplacementConvergence();
m4.makeCustomDisplacementColorPlot([],[],0.3,9,1,0);

individual_residual_flag = 1;
median_filter_threshold = [];
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residual_flag,shading_type);
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.saveObject();

m4.makeInteractiveDisplacementDarkFieldPlot();
m4.makeInteractiveInverseDarkFieldPlot();
