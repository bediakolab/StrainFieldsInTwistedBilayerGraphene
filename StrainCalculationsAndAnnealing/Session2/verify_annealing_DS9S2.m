% verify_annealing_DS9S2.m
%
% Phase unwrapping script for the session 2 data.
% These scripts no longer contain strain mapping functionality; that has
% been moved to /Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Filters/TGVstrainmappingscripts
%
% Nathanael Kazmierczak, 06/22/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/Session2');
% addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/ChiralFixed');
% load('06022020_DS6temp_preparedforannealing.mat');
load('06222020_DS9S2_prepared_for_annealing.mat');
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/ChiralFixed';
m4.saved_plots_foldername = 'DS9S2';

%% This works decently well we know
filteramp = 0.3;
filterrange = 5;
m4.setEmitterOrientations('v1',true);
m4.annealDisplacementField2(filteramp,filterrange,1);

method = 'outside input';
% AA_type_string = 'AA Gaussian ellipse';
% window_size = 80;
% window_increment = 10;
AA_type_string = [];
window_size = [];
window_increment = [];
crop_val_pixels = [];
angle_input = 0.6632;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

