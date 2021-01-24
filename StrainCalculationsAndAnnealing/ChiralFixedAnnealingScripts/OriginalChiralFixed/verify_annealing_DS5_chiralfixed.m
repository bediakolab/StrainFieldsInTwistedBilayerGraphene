% verify_annealing_DS5.m
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/ChiralFixed');
load('05312020_chiralfixed_DS5_forannealing.mat');
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/ChiralFixed';
m4.saved_plots_foldername = 'DS5_05312020';

% Filteramp and filterrange should be populated by the loaded .mat file.
filteramp = 0.3;
filterrange = 5;
direction = 'v1';
m4.setEmitterOrientations(direction);
m4.annealDisplacementField2(filteramp,filterrange,1);

method = 'outside input';
AA_type_string = [];
window_size = [];
window_increment = [];
% window_size = 80;
% window_increment = 10;
crop_val_pixels = 30;
angle_input = 1.2343;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

%% This works decently well we know

% As a good first approximation
saveplotflag = 1;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
sign_convention = [+1,-1];
trim_edge_pixel_num = [35,0,0];
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);
