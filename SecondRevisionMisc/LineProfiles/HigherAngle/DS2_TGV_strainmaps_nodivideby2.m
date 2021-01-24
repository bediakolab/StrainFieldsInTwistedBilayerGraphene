% DS2_TGV_strainmaps_nodivideby2.m
%
% Function testing the total generalized variation filters Colin found and
% trying to get them tuned correctly.
%
% Nathanael Kazmierczak, 06/13/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
    load('05312020_chiralfixed_DS2_postannealing.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS2_tgv';
end

% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct = buildFilterStruct(stringcell);
TGVID = 'DS2';  % versus Colin's original parameter settings.
rot_calibration_id = 'SP';
tensor_rotate_SPid = 1;
colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 0;
m4.computeStrainMaps3(30,pre_filterstruct,TGVID,rot_calibration_id,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);

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


figure;
imagesc(m4.xbase_for_strainmaps,m4.ybase_for_strainmaps,m4.dfield_filtered_for_strain(:,:,1)); set(gca,'yDir','normal'); axis equal; title('xdisp');
colormap pink
figure
imagesc(m4.xbase_for_strainmaps,m4.ybase_for_strainmaps,m4.dfield_filtered_for_strain(:,:,2)); set(gca,'yDir','normal'); axis equal; title('ydisp');
colormap pink
