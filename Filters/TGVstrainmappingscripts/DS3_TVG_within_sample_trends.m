% DS3_TVG_within_sample_trends.m
%
% This script makes AA, AB, and SP trends within a single dataset.
% Each value is assigned to the surrounding 
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

%% FBR trends without divide by two
TGVid = 'DS3';  % versus Colin's original parameter settings.
SPdirection = 1;
% colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [0,0];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 5;
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

[AAaverages,AAstds,AAmeantwistangles,AAstdtwistangles] = m4.getTriangulationStrainTrends(1);
% Make FBR plot for AA regions
figure;
errorbar(AAmeantwistangles,AAaverages(:,6),AAstds(:,6),AAstds(:,6),AAstdtwistangles,AAstdtwistangles,'o');
hold on
errorbar(AAmeantwistangles,AAaverages(:,6)-AAmeantwistangles,AAstds(:,6),AAstds(:,6),AAstdtwistangles,AAstdtwistangles,'o');
plot(AAmeantwistangles,AAmeantwistangles,'-');
title('AA fixed body rotations in DS3 triangulation');
xlabel('Triangulation Moire twist angle (deg)');
ylabel('Rotation (deg)');
legend('Total AA fixed-body rotation','Net reconstruction AA fixed-body rotation','Moire twist angle');


% Build combined plot also using the AA FBR from DS6
% Different colorization
figure;
errorbar(AAmeantwistangles,AAaverages(:,6),AAstds(:,6),AAstds(:,6),AAstdtwistangles,AAstdtwistangles,'o');
hold on
errorbar(AAmeantwistangles,AAaverages(:,6)-AAmeantwistangles,AAstds(:,6),AAstds(:,6),AAstdtwistangles,AAstdtwistangles,'o');
errorbar(DS6_AAmeantwistangles,DS6_AAaverages(:,6),DS6_AAstds(:,6),DS6_AAstds(:,6),DS6_AAstdtwistangles,DS6_AAstdtwistangles,'o');
errorbar(DS6_AAmeantwistangles,DS6_AAaverages(:,6)-DS6_AAmeantwistangles,DS6_AAstds(:,6),DS6_AAstds(:,6),DS6_AAstdtwistangles,DS6_AAstdtwistangles,'o');
plot([AAmeantwistangles;DS6_AAmeantwistangles],[AAmeantwistangles;DS6_AAmeantwistangles],'-');
title('AA fixed body rotations in DS3 and DS6 triangulation');
xlabel('Triangulation Moire twist angle (deg)');
ylabel('Rotation (deg)');
legend('DS3 Total AA fixed-body rotation','DS3 Net reconstruction AA fixed-body rotation','DS6 Total AA fixed-body rotation','DS6 Net reconstruction AA fixed-body rotation','Moire twist angle');

% Same colorization

% % % %% Strain trends with divide by two
% % % m4.saved_plots_foldername = 'DS6_tgv_nodividebytwo';
% % % TGVid = 'DS3';  % versus Colin's original parameter settings.
% % % SPdirection = 1;
% % % % colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
% % % divide_by_two_for_intralayer = false;
% % % sign_convention = [0,0];
% % % % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
% % % trim_tear_value = 5;
% % % trim_value = [0,10];
% % % m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);
% % % 
