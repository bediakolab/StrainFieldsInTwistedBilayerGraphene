% verify_annealing_DS26.m
%
% Using the new deterministic annealing algorithm (fast) to construct an
% extended-zone vector field for dataset 26 that will be comparable to the
% ones we constructed for datasets 4 and 7.
%
% Nathanael Kazmierczak, 04/27/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealingStudies');
load('04062020AnnealingStartingData.mat');
load('demo_colormap.mat');
% m4.scan_stepsize  % Already correctly set.
m4.annealed_dfield = [];
m4.annealDisplacementField2(0.3,9,1);


trim_edge_pixel_num = [30,0,10];
median_filter_threshold = [];
TVpre_fraction = [0.01,0.01];
TVpost_fraction = 0.90;
% TVpre_fraction = [0.1,0.1];
% TVpost_fraction = 0.85;
m4.makeStrainMapsFromBlinking(0,'flat',cmap,trim_edge_pixel_num,median_filter_threshold,TVpre_fraction,TVpost_fraction);
