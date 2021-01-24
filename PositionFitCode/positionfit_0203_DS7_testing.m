% positionfit_0203_DS7_testing.m
%
% Trying both the blinking and the position strain analysis methods.

filename = '7_80x80_ss=1nm_C2=40um_alpha=1mrad_spot7_100ms_CL=130 bin=4_60kV.h5';
probe_filename = []; % Replace this once Maddie uploads it. '60kV_1mrad_spot9_0p5s_130CL_defocus = 0.dm3';
currentdir = pwd;
cd ..
addpath(genpath(pwd));
cd(currentdir);
addpath('G:\Maddie 4D STEM\20200203');
results_path = 'G:\Maddie 4D STEM\20200203\Results';

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
% disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
disks_name = [];
stepsize = 1; % given the RS binning of 2
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);

% explicit_chunks = [15,16]
m4.fitElasticPowerLaw([]);
m4.plotPowerLawMask();
m4.setBlinkingDisks();
subtract_power_law_flag = 1;
m4.integrateDisks(subtract_power_law_flag);
saveplotflag = 1;
trig_prefactor = 0.95;
individual_residuals = 1;
shadingtype = 'flat';
m4.makeBlinkingPlots([],saveplotflag,shadingtype);
m4.fitBlinking(trig_prefactor); % no override
median_filter_threshold = 1.0;
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residuals);
saveplotflag = 0;
colormaptype = 'hsv';
trimwidth = 0;

m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth);
m4.saveObject();
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.makeInteractiveDisplacementDarkFieldPlot();