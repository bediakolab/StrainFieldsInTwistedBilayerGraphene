% positionfit_1023_DS1_testing.m
%
% Using the dataset where a real space binning of 2 was applied in
% py4DSTEM. Note that because no additionally binning has taken place in
% reciprocal space, we can load in the disks directly from the previous
% analysis of the full-resolution dataset 7.
%
% Nathanael Kazmierczak, 01/29/2020

filename = '1_MV_1p5_9_40x40_ss0p5nm_1s_spot9_alpha=1_bin2_cl=130_60kV_normalized.h5';
probe_filename = '60kV_1mrad_spot9_0p5s_130CL_defocus = 0.dm3';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/New2019_1023/Dataset1');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/PositionFit/Results/1023_dataset1_positionfit_try1';

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
% disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
disks_name = [];
stepsize = 0.5; % given the RS binning of 2
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);
% explicit_chunks = [15,16];
m4.setBlinkingDisks();
m4.integrateDisks();
saveplotflag = 1;
trig_prefactor = 0.9;
individual_residuals = 1;
shadingtype = 'flat';
m4.makeBlinkingPlots([],saveplotflag,shadingtype);
m4.fitBlinking(saveplotflag,trig_prefactor,individual_residuals);
saveplotflag = 1;
colormaptype = 'hsv';
trimwidth = 2;
m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth);
m4.saveObject();
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.makeInteractiveDisplacementDarkFieldPlot();
