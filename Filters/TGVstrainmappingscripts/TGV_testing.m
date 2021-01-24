% TGV_testing.m
%
% Function testing the total generalized variation filters Colin found and
% trying to get them tuned correctly.
%
% Nathanael Kazmierczak, 06/13/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
    load('05312020_chiralfixed_noprefactor_DS26_postannealing.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS26_tgv';
end

% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct = buildFilterStruct(stringcell);
useNPKsettings = false;  % versus Colin's original parameter settings.
SPdirection = 1;
colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = true;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 0;
m4.computeStrainMaps3(30,pre_filterstruct,useNPKsettings,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

saveplotflag = false;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);         

%% Comparison to my original function.
if false
    saveplotflag = false;
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    colormaptype = cmap;
    trim_edge_pixel_num = [30,0,0];
    filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
    sign_convention = [0,0];
    trim_tear_pixel_num = 0;
    overlay_registration = false;
    divide_by_two_for_intralayer = false;
    % SP_num_for_xaxis = 0;  % implying unrotated strain mape.
    [eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
        m4.makeStrainMapsFromBlinking2(saveplotflag,colormaptype,trim_edge_pixel_num,...
        filterstruct,sign_convention,trim_tear_pixel_num,overlay_registration,divide_by_two_for_intralayer,SPdirection);
    
    
end
