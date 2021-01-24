% DS8S2_TGV_strainmaps.m
%
% Strain map function built for Session 2 data. Strain
% mapping and annealing are now found in separate functions.
%
% Nathanael Kazmierczak, 06/22/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2');
%     load('05312020_chiralfixed_DS8_postannealing.mat');
    load('06222020_DS8S2_annealed.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS8S2';
end

% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct = buildFilterStruct(stringcell);
% Testing with the heavy mean.
% [ pre_filterstruct ] = buildFilterStructPreset( 'DS26 annealed displacement #2' );
% pre_filterstruct = [];

TGVid = 'DS8S2';  % versus Colin's original parameter settings.
SPdirection = 1;
colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = true;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = [];
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

saveplotflag = true;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);         


m4.saved_plots_foldername = 'DS30S2_nodividebytwo';
TGVid = 'DS8S2';  % versus Colin's original parameter settings.
SPdirection = 1;
colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = [];
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

saveplotflag = true;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);   
