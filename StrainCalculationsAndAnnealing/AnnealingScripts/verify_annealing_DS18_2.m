% verify_annealing_DS18_2.m
%
% Using a smaller-range filter for the two 1 nm stepsize datasets to see if
% this brings them into consistency with the other 0.5 nm datasets.
%
% Annealing proper is not needed here, just re-computation of the strain
% maps.
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload');
load('dataset18_deterministicanneal_05242020.mat');
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults';
m4.saved_plots_foldername = 'DS18_2_05242020';

%% This works decently well we know

% As a good first approximation
saveplotflag = 1;
filteramp = 0.3;
filterrange = 9;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
filterstruct = buildFilterStructPreset('DS15 annealed displacement #1');
sign_convention = [-1,-1];
trim_edge_pixel_num = [50,0,0];
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);

% The untrimmed may be closer to the center of the AA domain, as I recall
[averages_mat,stds_mat,n_values,names] = m4.getStrainInAARingedBufferMasks('untrimmed AA','circular',1,0)

