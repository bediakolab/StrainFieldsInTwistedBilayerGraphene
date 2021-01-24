% verify_annealing_DS3.m
%
% Using a smaller-range filter for the two 1 nm stepsize datasets to see if
% this brings them into consistency with the other 0.5 nm datasets.
%
% Annealing proper is not needed here, just re-computation of the strain
% maps.
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
load('06012020_chiralfixed_tear_DS3_postannealing.mat');
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/ChiralFixed';
m4.saved_plots_foldername = 'DS3_06012020';

%% This works decently well we know
filteramp = 0.3;
filterrange = 9;
% m4.setEmitterOrientations();
% m4.annealDisplacementField2(filteramp,filterrange,1);

% % % method = 'outside input';
% % % % AA_type_string = 'AA Gaussian ellipse';
% % % % window_size = 80;
% % % % window_increment = 10;
% % % AA_type_string = [];
% % % window_size = [];
% % % window_increment = [];
% % % crop_val_pixels = [];
% % % angle_input = 0.1203;
% % % m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);


% As a good first approximation
saveplotflag = 0;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
sign_convention = [0,0];
trim_edge_pixel_num = [10,0,0];
trim_tear_boundary_pixelnum = 2;
overlay = true;
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention,trim_tear_boundary_pixelnum,overlay);

% % The untrimmed may be closer to the center of the AA domain, as I recall
% [averages_mat,stds_mat,n_values,names] = m4.getStrainInAARingedBufferMasks('untrimmed AA','circular',1,0)

