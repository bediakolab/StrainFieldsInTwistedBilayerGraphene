% team1_ds7_script.m

% addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset7');
addpath('E:\BediakoLab\Maddie4DSTEM\20200226\Dataset7');
addpath(genpath('H:\matlab4DSTEM\'));
filename = 'TEAM1DS7DiffractionSI.h5';
disk_filename = [];
resultsfolderpath = 'a/b';
scan_stepsize = 1;
probe_filename = [];
m4 = FourDSTEM_Analysis_Engine(filename,disk_filename,resultsfolderpath,scan_stepsize,probe_filename);
m4.saved_plots_folderpath = 'E:\BediakoLab\Maddie4DSTEM\NPKResults';
m4.saved_plots_foldername = 'Dataset7';
m4.setBlinkingDisks();

m4.fitElasticPowerLaw();

subtract_power_law_flag = 1;
use_average_flag = 0;
m4.integrateDisks(subtract_power_law_flag,use_average_flag);

