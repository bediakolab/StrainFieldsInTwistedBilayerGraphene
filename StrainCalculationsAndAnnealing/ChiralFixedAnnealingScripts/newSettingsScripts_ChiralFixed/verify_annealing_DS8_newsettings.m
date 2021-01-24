% verify_annealing_DS8_newsettings.m
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
load('05312020_chiralfixed_DS8_postannealing.mat');
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/ChiralFixed';
m4.saved_plots_foldername = 'DS8_06032020_newsettings';

% Filteramp and filterrange should be populated by the loaded .mat file.
filteramp = 0.3;
filterrange = 5;
% m4.setEmitterOrientations();
% m4.annealDisplacementField2(filteramp,filterrange,1);

method = 'outside input';
AA_type_string = [];
window_size = [];
window_increment = [];
% window_size = 80;
% window_increment = 10;
crop_val_pixels = 25;
angle_input = 1.3672;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

%% This works decently well we know

% As a good first approximation
saveplotflag = 1;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
sign_convention = [+1,-1];
trim_edge_pixel_num = [30,0,0];
trim_tear_pixel_num = [];
overlay_registration = false;
divide_by_two_for_intralayer = true;
SP_num_for_xaxis = 1;
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention,trim_tear_pixel_num,...
                overlay_registration,divide_by_two_for_intralayer,SP_num_for_xaxis);
            
overlay_registration = true;          
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention,trim_tear_pixel_num,...
                overlay_registration,divide_by_two_for_intralayer,SP_num_for_xaxis);
            