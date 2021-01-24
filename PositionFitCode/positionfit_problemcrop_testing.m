% positionfit_problemcrop_testing.m
%
% Trying both the blinking and the position strain analysis methods.

filename = 'Diffraction SI crop 20x20.h5';
probe_filename = []; % Replace this once Maddie uploads it. '60kV_1mrad_spot9_0p5s_130CL_defocus = 0.dm3';
currentdir = pwd;
cd ..
addpath(genpath(pwd));
cd(currentdir);
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/20200226';

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
% disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
disks_name = [];
stepsize = 1; % given the RS binning of 2
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);

% explicit_chunks = [15,16]
% m4.fitElasticPowerLaw([]);
% m4.plotPowerLawMask();
% m4.eliminateHBN();
m4.setBlinkingDisks([]);
% subtract_power_law_flag = 1;
subtract_power_law_flag = 0;
m4.integrateDisks(subtract_power_law_flag);
saveplotflag = 1;
trig_prefactor = 1;
individual_residuals = 1;
shadingtype = 'flat';
m4.makeBlinkingPlots([],saveplotflag,shadingtype);
% weight_vector_override = [1,1,0,1,1,0,1,1,1,1,1,1];
m4.fitBlinking(trig_prefactor); % no override
median_filter_threshold = [];
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residuals,'flat');
saveplotflag = 0;
colormaptype = 'hsv';
trimwidth = 0;

% m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth);
m4.saveObject();
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.makeInteractiveDisplacementDarkFieldPlot();
