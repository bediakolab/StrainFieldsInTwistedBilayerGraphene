% ringed_Gaussian_testing_script.m
%
% Nathanael Kazmierczak, 04/09/2020


filename = 'TEAM1Day1Dataset4DiffractionSI.h5';
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset4');
probe_filename = [];
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/BackgroundFitTesting/Dataset4RingedGaussian';
% This is the manual way of setting where the pictures will be saved to. If
% you don't want automated saving of pictures, it really doesn't matter!
disks_name = [];
stepsize = 0.5; 
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m4.setAveragedDP('all');
m4.setAveragedDP(1);
m4.plotAveragedDP(1);
m4.makeMasksForElasticFit();
m4.setBlinkingDisks();
m4.fitElasticScattering('ringed gaussian');
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

m4.integrateDisks(subtract_elastic_flag,use_average_flag,radial_elastic_correction_flag,interpolation_correction_flag);  % omit loadnumcroprange argument for default, loading all.
m4.makeBlinkingPlots([],saveplotflag,shading_type,normalize);
m4.saveObject();

m4.setTrigHeight();
fit_boundaries = {[1,25],[1,25]};
useParallel = false;
useAnalyticJacobian = true;
useFittedPrefactors = false;
m4.fitBlinking3([],useParallel,fit_boundaries,useAnalyticJacobian,useFittedPrefactors);
m4.saveObject();
m4.makeCustomDisplacementColorPlot([],[],[],[],1,0);

m4.fitBlinking4([],useParallel,[25,25]);
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
