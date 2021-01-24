% DS11S2_rotcalibrated_strainmaps.m
%
% Strain map function built for Session 2 data. Strain
% mapping and annealing are now found in separate functions.
%
% The new series of rotationally calibrated strain maps performs several
% new functions. (1) Computes the purple-soliton wall aligned strain map
% (same as before, unneeded for DS8S2 and DS9S2, (2) computes the
% rotationally calibrated strain maps, (3) returns the various rotational
% alignment parameters for storage. Eventually, this should also include
% conversion of uniaxial analysis back into a realspace strain axis
% description.
%
% Nathanael Kazmierczak, 06/22/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2');
    %     load('05312020_chiralfixed_DS8_postannealing.mat');
    load('06222020_DS11S2_annealed.mat');
    load('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
    m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/PostColin';
    m4.saved_plots_foldername = 'DS11S2';
end

lightTGVid = 'DS11S2';
heavyTGVid = 'DS26';

% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct = buildFilterStruct(stringcell);

%% Get both SP and nanoparticle rotation values
m4.scan_direction_theta_deg = 180;
m4.getDiffractionPatternRotation('integration disks');
[nanoparticle_rotation_deg,SProtation_deg] = m4.getRealSpaceRotationAngle()
nanoparticle_rotation_deg_reduced = mod(nanoparticle_rotation_deg,180)
% Check manually and make sure these values are sensible. If not, check the
% scan direction against the spreadsheet, for starters.

%% Perform strain mapping with the new nanoparticle rotational calibration
m4.saved_plots_foldername = 'DS11S2_nanoparticle_calibration_SP1';
TGVid = lightTGVid;  % versus Colin's original parameter settings.
rotational_calib_method = 'nanoparticle mod';
tensor_rotate_SPid = 1;  % after setting the x-axis via rotational calibration, rotate tensor to SP1.
colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = true;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = [];
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);

saveplotflag = true;
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);



%% Can do this with a heavier filter as well
m4.saved_plots_foldername = 'DS11S2_nanoparticle_calibration_DS26heavyfilter_SP1';
TGVid = heavyTGVid;  % Just to see what it will look like
rotational_calib_method = 'nanoparticle mod';
tensor_rotate_SPid = 1;  % after setting the x-axis via rotational calibration, rotate to tensor.
colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = true;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = [];
trim_value = [0,10];
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);

saveplotflag = true;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
overlay_registration = true;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);

%% Get strain maps for a different tensor rotation
if true
    close all
    m4.saved_plots_foldername = 'DS11S2_nanoparticle_calibration_DS26heavyfilter_SP2';
    TGVid = heavyTGVid;  % Just to see what it will look like
    rotational_calib_method = 'nanoparticle mod';
    tensor_rotate_SPid = 2;  % after setting the x-axis via rotational calibration, rotate to tensor.
    colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = [];
    trim_value = [0,10];
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);
    
    saveplotflag = true;
    overlay_registration = false;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    overlay_registration = true;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    
    
    m4.saved_plots_foldername = 'DS11S2_nanoparticle_calibration_DS26heavyfilter_SP3';
    TGVid = heavyTGVid;  % Just to see what it will look like
    rotational_calib_method = 'nanoparticle mod';
    tensor_rotate_SPid = 3;  % after setting the x-axis via rotational calibration, rotate to tensor.
    colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = [];
    trim_value = [0,10];
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);
    
    saveplotflag = true;
    overlay_registration = false;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    overlay_registration = true;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    
    
    
    % More tensor rotated strain maps, normal filter
    close all
    m4.saved_plots_foldername = 'DS11S2_nanoparticle_calibration_SP2';
    TGVid = lightTGVid;  
    rotational_calib_method = 'nanoparticle mod';
    tensor_rotate_SPid = 2;  % after setting the x-axis via rotational calibration, rotate to tensor.
    colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = [];
    trim_value = [0,10];
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);
    
    saveplotflag = true;
    overlay_registration = false;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    overlay_registration = true;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    
    
    m4.saved_plots_foldername = 'DS11S2_nanoparticle_calibration_SP3';
    TGVid = lightTGVid;
    rotational_calib_method = 'nanoparticle mod';
    tensor_rotate_SPid = 3;  % after setting the x-axis via rotational calibration, rotate to tensor.
    colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = [];
    trim_value = [0,10];
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);
    
    saveplotflag = true;
    overlay_registration = false;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    overlay_registration = true;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
end


%% Perform strain mapping as before with the soliton rotational calibration
if false   % will return to this for other datasets, but 8S2 is already computed.
    m4.saved_plots_foldername = 'DS8S2';
    TGVid = 'DS8S2';  % versus Colin's original parameter settings.
    rotational_calib_method = 'SP';
    colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = [];
    trim_value = [0,10];
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value);
    
    saveplotflag = true;
    overlay_registration = false;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    overlay_registration = true;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    
    
    m4.saved_plots_foldername = 'DS8S2_nodividebytwo';
    TGVid = 'DS8S2';  % versus Colin's original parameter settings.
    colormaptype = {cmap,cmap};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = false;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = [];
    trim_value = [0,10];
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,rotational_calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value);
    
    saveplotflag = true;
    overlay_registration = false;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
    overlay_registration = true;
    m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);
end