% team1_ds7_radialresidualsinterp_script.m

% addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset7');
% addpath('E:\BediakoLab\Maddie4DSTEM\20200226\Dataset7');
% addpath(genpath('H:\matlab4DSTEM\'));
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset7');
filename = 'TEAM1DS7DiffractionSI.h5';
disk_filename = [];
resultsfolderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset7RadialResidualsInterp_Fulldataset';
scan_stepsize = 1;
probe_filename = [];
m4 = FourDSTEM_Analysis_Engine(filename,disk_filename,resultsfolderpath,scan_stepsize,probe_filename);
% m4.saved_plots_folderpath = 'E:\BediakoLab\Maddie4DSTEM\NPKResults';
% m4.saved_plots_foldername = 'Dataset7';


m4.setAveragedDP();
m4.plotAveragedDP(1);
m4.makeMasksForElasticFit();
m4.fitElasticScattering('lorentzian');
m4.plotElasticMask(0);
m4.plotElasticMask(1);
m4.saveObject();

saveplotflag = 1;
shading_type = 'flat';
normalize = 0;
subtract_elastic_flag = 1;
use_average_flag = 1;
radial_elastic_correction_flag = 1;
m4.setBlinkingDisks();
m4.integrateDisks(subtract_elastic_flag,use_average_flag,radial_elastic_correction_flag);  % omit loadnumcroprange argument for default, loading all.
m4.saveObject();

m4.setTrigHeight();
m4.fitBlinking3([],0,[]);
m4.saveObject();

individual_residual_flag = 1;
median_filter_threshold = [];
m4.makeBlinkingPlots([],saveplotflag,shading_type,normalize);
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residual_flag,shading_type);
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.saveObject();