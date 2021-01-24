% TGV_testing_DS8.m
%
% Function testing the total generalized variation filters Colin found and
% trying to get them tuned correctly. 
% 
% Note TGV_testing function operates on DS26.
%
% Nathanael Kazmierczak, 06/13/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
%     load('05312020_chiralfixed_DS8_postannealing.mat');
    load('06012020_chiralfixed_tear_DS3_postannealing.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS3_tgv_dividebytwo';
end

% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct = buildFilterStruct(stringcell);
% Testing with the heavy mean.
% [ pre_filterstruct ] = buildFilterStructPreset( 'DS26 annealed displacement #2' );
% pre_filterstruct = [];
rot_calibration_id = 'SP';
tensor_rotate_SPid = 1;
TGVid = 'DS3';  % versus Colin's original parameter settings.
colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [0,0];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 5;
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rot_calibration_id,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);

saveplotflag = false;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);         


m4.saved_plots_foldername = 'DS3_tgv_nodividebytwo';
TGVid = 'DS3';  % versus Colin's original parameter settings.
SPdirection = 1;
colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [0,0];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 5;
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

saveplotflag = true;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);   
