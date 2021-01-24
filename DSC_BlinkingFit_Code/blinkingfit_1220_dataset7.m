% blinkingfit_1220_dataset7.m

filename = '7_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130 bin=2_scan90_60kV_NPKNormalized_andcropped.h5';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/2019_1220');

m4 = FourDSTEM_Analysis_Engine(filename);
m4.setBlinkingDisks();
m4.integrateDisks();
m4.makeBlinkingPlots();
m4.fitBlinking();


