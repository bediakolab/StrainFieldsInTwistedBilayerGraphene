% Script_for_data_fitting.m
%
% NPK and MVW, 02/09/2020

%% All sorts of setup
% filename = '9_40x40_ss=1nm_C2=40um_alpha=1mrad_spot7_100ms_CL=130 bin=4_rot=90_60kV.h5';
filename = 'Diffraction SI.h5';
% filename = 'Diffraction SI.h5';
probe_filename = []; % Replace this once Maddie uploads it. '60kV_1mrad_spot9_0p5s_130CL_defocus = 0.dm3';
currentdir = pwd;
cd ..
addpath(genpath(pwd));
cd(currentdir);
% addpath('G:\Maddie 4D STEM\20200203');
addpath('D:\Maddie Local 2\4D STEM\20200226\4-200x200-0.5nm small unit cell with filter');
results_path = 'C:\Users\Bediako Lab\Desktop\Maddie local\Maddie 4D STEM\20200226\Results\new';
% results_path = 'D:\Maddie Local 2\4D STEM\20200226\7-300x300-0.5nm large unit cell_NPKResults';
% This is the manual way of setting where the pictures will be saved to. If
% you don't want automated saving of pictures, it really doesn't matter!

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
% disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
disks_name = [];
stepsize = 1; % given the RS binning of 2
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);
m4.saved_plots_folderpath = 'C:\Users\Bediako Lab\Desktop\Maddie local\Maddie 4D STEM\20200226\Results';
m4.saved_plots_foldername = 'new';


%% Define disks for interferometry 
m4.setBlinkingDisks();
subtract_power_law_flag = 0;
m4.integrateDisks(subtract_power_law_flag);

%% Make blinking plot
saveplotflag = 1;
shadingtype = 'flat';
m4.makeBlinkingPlots([],saveplotflag,shadingtype);

%% Fit displacement field to blinking
trig_prefactor = 0.95;
weight_vector_override = [];
useParallel = 1;
m4.fitBlinking(trig_prefactor,weight_vector_override,useParallel); 

%% Make displacement pictures:
individual_residuals = 1;
% median_filter_threshold = [0,0];
% median_filter_threshold = 1.0;
median_filter_threshold = [1,1];
shading_type = 'flat';
m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_threshold,individual_residuals,shading_type);

%% Strain, quiver, interactiveStuff
colormaptype = 'gray';
trimwidth = 2;
m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth,median_filter_threshold);
m4.saveObject();
m4.makeQuiverPlotFromBlinking(saveplotflag);
m4.makeInteractiveDisplacementDarkFieldPlot();


