% verify_annealing_DS4.m
%
% Nathanael Kazmierczak, 04/26/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets');
load('04252020DS4_Day1_GeometryEmitterPreparedforAnnealing.mat');
load('demo_colormap.mat');
% m4.scan_stepsize  % Already correctly set.
m4.annealDisplacementField2(0.3,5,1);
%
%
% trim_edge_pixel_num = [30,5,50];
% median_filter_threshold = [];
% TVpre_fraction = [0.01,0.01];
% % TVpost_fraction = 0.95;
% TVpost_fraction = 0.9;
% m4.makeStrainMapsFromBlinking(0,'flat',cmap,trim_edge_pixel_num,median_filter_threshold,TVpre_fraction,TVpost_fraction);

%   changeOfBasis = 0;
%             METHOD = 'old';
%             MEAN = true;
%             FILTER = true;
%             MEDFILTTHRES = true;
%             SMOOTH = false;
%             scale_color_map_flag = true;
%             vonMises = true;
%             PAD = true;
%             ZEROFILL = true;


%% This works decently well we know

trim_edge_pixel_num = [30,0,10];
median_filter_threshold = [];
TVpre_fraction = [0.01,0.01];
% TVpost_fraction = 0.95;
TVpost_fraction = 0.90;
m4.makeStrainMapsFromBlinking(0,'flat',cmap,trim_edge_pixel_num,median_filter_threshold,TVpre_fraction,TVpost_fraction);

%  
%             changeOfBasis = 0;
%             METHOD = 'old';
%             MEAN = true;
%             FILTER = true;
% %             MEDFILTTHRES = true;
%             SMOOTH = false;
%             scale_color_map_flag = true;
%             vonMises = true;
