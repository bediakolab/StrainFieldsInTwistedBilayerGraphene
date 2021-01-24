% verify_annealing_DS10.m
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets');
load('05232020_DS10_GeometryFit_PreparedForAnnealing.mat');
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults';
m4.saved_plots_foldername = 'DS10_05232020';

% Filteramp and filterrange should be populated by the loaded .mat file.
m4.setEmitterOrientations();
m4.annealDisplacementField2(filteramp,filterrange,1);

method = 'outside input';
AA_type_string = [];
window_size = [];
window_increment = [];
% window_size = 80;
% window_increment = 10;
crop_val_pixels = 30;
angle_input = 0.6272;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

%% This works decently well we know

% As a good first approximation
saveplotflag = 1;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
sign_convention = [+1,+1];
trim_edge_pixel_num = [30,0,0];
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);
