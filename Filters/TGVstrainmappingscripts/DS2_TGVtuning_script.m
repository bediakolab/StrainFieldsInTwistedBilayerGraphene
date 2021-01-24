% DS2_TGVtuning_script.m
%
% Script trying to find a new set of filter settings for dataset 2 so that
% we can make principal strain plots for a series of datasets.
%
% Nathanael Kazmierczak, 06/19/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
%     load('05312020_chiralfixed_DS8_postannealing.mat');
    load('05312020_chiralfixed_DS2_postannealing.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS2_tgv_tuned';
end

% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.15]};
pre_filterstruct = buildFilterStruct(stringcell);
% Testing with the heavy mean.
% [ pre_filterstruct ] = buildFilterStructPreset( 'DS26 annealed displacement #2' );
% pre_filterstruct = [];

TGVid = 'DS2testing';  % versus Colin's original parameter settings.
SPdirection = 1;
colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = true;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 0;
trim_value = [0,30];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

saveplotflag = false;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);

f1 = figure(13);
f2 = figure(32);
f1.HandleVisibility = 'off';
f2.HandleVisibility = 'off';
close all
f1.HandleVisibility = 'on';
f2.HandleVisibility = 'on';
figure(f2);
figure(f1);



