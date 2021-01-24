% verify_annealing_DS7.m
%
% Nathanael Kazmierczak, 04/26/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets');
load('04062020DS7_Day1_GeometryEmitterPreparedforAnnealing.mat');
load('demo_colormap.mat');
% m4.scan_stepsize  % Already correctly set.
m4.annealDisplacementField2(0.3,7,1);


trim_edge_pixel_num = [30,0,10];
median_filter_threshold = [];
TVpre_fraction = [0.01,0.01];
TVpost_fraction = 0.90;
% TVpre_fraction = [0.1,0.1];
% TVpost_fraction = 0.85;
m4.makeStrainMapsFromBlinking(0,'flat',cmap,trim_edge_pixel_num,median_filter_threshold,TVpre_fraction,TVpost_fraction);
