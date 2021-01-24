% team1_ds4_radialresidualsinterp_script.m

% addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset7');
% addpath('E:\BediakoLab\Maddie4DSTEM\20200226\Dataset7');
% addpath(genpath('H:\matlab4DSTEM\'));
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset4');
filename = 'TEAM1Day1Dataset7DiffractionSIcropped.h5';
disk_filename = [];
resultsfolderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset4RadialResidualsInterp_cropped';
scan_stepsize = 1;
probe_filename = [];
m4 = FourDSTEM_Analysis_Engine(filename,disk_filename,resultsfolderpath,scan_stepsize,probe_filename);
% m4.saved_plots_folderpath = 'E:\BediakoLab\Maddie4DSTEM\NPKResults';
% m4.saved_plots_foldername = 'Dataset7';

m4.setAveragedDP('all');
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
m4.setAveragedDP(1);
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
fit_boundaries = {[1,80],[1,80]};
m4.fitBlinking3([],0,fit_boundaries);
m4.saveObject();

individual_residual_flag = 1;
median_filter_threshold = [];
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residual_flag,shading_type);
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.saveObject();


%%% Old stuff from this file, left here in case we need it.
% 
% saveplotflag = 1;
% shading_type = 'flat';
% normalize = 0;
% subtract_elastic_flag = 1;
% use_average_flag = 0;
% radial_elastic_correction_flag = 1;
% m4.integrateDisks(subtract_elastic_flag,use_average_flag,radial_elastic_correction_flag);  % omit loadnumcroprange argument for default, loading all.
% m4.setTrigHeight();
% m4.fitBlinking3([],0,[]);
% 
% 
% individual_residual_flag = 1;
% median_filter_threshold = [];
% m4.makeBlinkingPlots([],saveplotflag,shading_type,normalize);
% m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residual_flag,shading_type);
% m4.makeQuiverPlotFromBlinking(saveplotflag);
