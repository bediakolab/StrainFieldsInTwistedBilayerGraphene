% blinkingfit_0123_dataset7_RSbin2_allouter.m
%
% Using the dataset where a real space binning of 2 was applied in
% py4DSTEM. Note that because no additionally binning has taken place in
% reciprocal space, we can load in the disks directly from the previous
% analysis of the full-resolution dataset 7.
%
% Nathanael Kazmierczak, 01/29/2020

filename = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized_RSbin2.h5';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200123');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/0123_dataset7_RSbin2_try4allouter';

thisd = pwd;
cd('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/0123_dataset7_RSbin2_try3newnegativitytreatment');
load('objectdata.mat');
m4old = m4;
clear m4
cd(thisd);

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
stepsize = 4; % given the RS binning of 2
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize);
% explicit_chunks = [15,16];

m4.disk_centers = m4old.disk_centers;
m4.disk_radius = m4old.disk_radius;
m4.skipped_disks = m4old.skipped_disks;
m4.power_law_fit_mask = m4old.power_law_fit_mask;
m4.hBN_mask = m4old.hBN_mask;
m4.beamstop_mask = m4old.beamstop_mask;
m4.graphene_mask = m4old.graphene_mask;
m4.detector_mask = m4old.detector_mask;
m4.power_law_fit_parameters = m4old.power_law_fit_parameters;
m4.optimal_power_law_fit = m4old.optimal_power_law_fit;
% fitElasticPowerLaw(obj,power_law_mask);

subtract_power_law_flag = 1;
m4.setBlinkingDisks();
m4.integrateDisks(subtract_power_law_flag);
saveplotflag = 1;
trig_prefactor = 0.9;
individual_residuals = 1;
shadingtype = 'flat';
m4.makeBlinkingPlots([],saveplotflag,shadingtype);
% The key reason for this script
weight_vector_override = [0,0,0,0,0,0,1,1,1,0,1,1];
m4.fitBlinking(saveplotflag,trig_prefactor,individual_residuals,weight_vector_override);
saveplotflag = 1;
colormaptype = 'hsv';
trimwidth = 2;
m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth);
m4.saveObject();
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.makeInteractiveDisplacementDarkFieldPlot();
