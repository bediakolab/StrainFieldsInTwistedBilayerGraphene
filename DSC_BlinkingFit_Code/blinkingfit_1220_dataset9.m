% blinkingfit_1220_dataset9.m

filename = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130 bin=2_scan35_60kV_NPK_Normalized.h5';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/2019_1220');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/1220_dataset9_prefactor0p8';

% disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
disks_name = 'ds9disks.mat';
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path);
m4.setBlinkingDisks();
m4.integrateDisks();
saveplotflag = 1;
trig_prefactor = 0.8;
individual_residuals = 1;
m4.makeBlinkingPlots([],saveplotflag);
m4.fitBlinking(saveplotflag,trig_prefactor,individual_residuals);


