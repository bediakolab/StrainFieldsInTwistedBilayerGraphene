% blinkingfit_0123_dataset7.m
%
% Nathanael Kazmierczak, 01/29/2020

filename = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized.h5';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200123');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/0123_dataset7';

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
stepsize = 2; % nm, see the filename
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize);
explicit_chunks = [15,16];
m4.setBlinkingDisks(explicit_chunks);
m4.integrateDisks();
saveplotflag = 0;
trig_prefactor = 0.9;
individual_residuals = 1;
m4.makeBlinkingPlots([],saveplotflag);
m4.fitBlinking(saveplotflag,trig_prefactor,individual_residuals);
saveplotflag = 1;
shadingtype = 'flat';
colormaptype = 'hsv';
trimwidth = 3;
m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth);
m4.saveObject();
