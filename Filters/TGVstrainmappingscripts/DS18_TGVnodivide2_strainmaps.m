% DS18_TGVnodivide2_strainmaps.m
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
    load('05312020_chiralfixed_DS18_2_postannealing.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS18_tgv_nodivideby2';
end

% very light outlier filtering to start off
% Note that the median filtering is a little bit lightened here.
stringcell = {'median outlier','both displacements',[3,1]};
pre_filterstruct = buildFilterStruct(stringcell);
% Testing with the heavy mean.
% [ pre_filterstruct ] = buildFilterStructPreset( 'DS26 annealed displacement #2' );
% pre_filterstruct = [];

TGVid = 'DS18';  % versus Colin's original parameter settings.
tensor_rotate_SPid = 1;
rotcalibration_ID = 'SP';
colormaptype = {cmap,inferno};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 0;
trim_value = [0,30];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotcalibration_ID,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);

saveplotflag = true;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);         

