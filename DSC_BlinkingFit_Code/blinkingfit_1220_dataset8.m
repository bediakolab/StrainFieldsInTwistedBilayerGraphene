% blinkingfit_1220_dataset8.m

filename = '8_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130 bin=2_scan0_60kV_NPK_Normalized_and_Cropped.h5';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/2019_1220');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/1220_dataset7';

disks_name = [];
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path);
m4.setBlinkingDisks();
m4.integrateDisks();
saveplotflag = 1;
m4.makeBlinkingPlots([],saveplotflag);
m4.fitBlinking(saveplotflag);


