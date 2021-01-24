% DS1_cropped_strain.m
%
% Nathanael Kazmierczak, 11/09/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo');
load('06022020_chiralfixed_DS1_postannealing.mat');
load('demo_colormap.mat');

calculate_strain_maps = false;
plot_strain_maps = false;
AA_strain_statistics = false;
pseudstacking_statistics = false;
SP_strain_statistics = false;
AB_strain_statistics = false;
geometry_statistics = true;

SPdirection = 1;

% Get the twist angle for the cropped region
crop_region = {[10,100],[100,190]};
strain_crop_range = {91:180,1:90}; % because the edges have already been taken off by this point by the strain filter.
% crop_region = [];
% strain_crop_range = [];
AA_type = 2;  % Gaussian ellipse
use_manual_triangulation = false;
m4.triangulateAA(AA_type,use_manual_triangulation,crop_region);
m4.repairTriangulation(false);
m4.setFacetTwistAngles();
[average_angle,std_angle,max_angle,min_angle,n_triangles_used,angles_used] = ...
    m4.plotFacetedTwistAngles([-inf,inf],inferno,true,false,true,[])

% obj.moire_angle_estimates
method = 'outside input';
AA_type_string = [];
window_size = [];
window_increment = [];
% window_size = 80;
% window_increment = 10;
crop_val_pixels = 0;
angle_input = 0.3991;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);


if calculate_strain_maps
% The lighter prefilterstruct because of the larger step size.
stringcell = {'median outlier','both displacements',[3,1]};
pre_filterstruct = buildFilterStruct(stringcell);


tensor_calib_setting = 'SP';
TGVid = 'DS15quant';
%     colormaptype = {cmap,fire};  % option to have one colormap for principal strain.

%%% This needs to be set to "false" for calculating the rotation angle, but
%%% "true" for calculating the strain itself.
divide_by_two_for_intralayer = false;  % Changing this to false is key for the revised calculation and plotting.
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 0;
trim_value = [0,10];  % [0,10]
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,tensor_calib_setting,divide_by_two_for_intralayer,sign_convention,trim_tear_value,SPdirection);
end
if plot_strain_maps
saveplotflag = false;
overlay_registration = true;
colormaptype = {cmap,fire};

m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype,strain_crop_range);
end

if AA_strain_statistics
id = 'Gaussian ellipse AA fit';
AA_radius = 1;
[averages_mat,stds_mat,n_values,names,stderrs,nconn] = m4.getStrainInAARingedBufferMasks(id,'circular',AA_radius,false,crop_region,strain_crop_range);
end




% Perform the pseudostacking analysis
if pseudstacking_statistics
%     ps_crop_range = {101:200,1:100};
ps_crop_range = strain_crop_range;
use_annealed = 2; % New setting, makes other things unused because references operations by computeStrainMaps3.
croprange = [];  % unused
filterstruct = [];
%     r = 0.3;
r = 'modified wigner-seitz';
make_plots = true;
[AA_perc,AB_perc,SP_perc,pseudostacking_mat] = m4.getPseudostackingAssignments(make_plots,r,filterstruct,use_annealed,croprange,[],ps_crop_range);
end



if SP_strain_statistics
    AAremovalradius = 6;
%     SPdirection = 2;
    SPbufferdist = 1; % nm
%     useStrainFilterMatchedSP = 0;
    useStrainFilterMatchedSP = false;
    id = 'Gaussian ellipse AA fit';
    [values,stds,npixels,names,stderrs,nconn] = m4.getStrainInSPBufferMask(SPdirection,SPbufferdist,...
        AAremovalradius,id,useStrainFilterMatchedSP,strain_crop_range);
end


if AB_strain_statistics
    SPbufferdist = 5;  % DS26S1, which had similar dimensions, had 11, but this appears to be a pixel count and not nm and here the 1nm step size
    AAremovalradius = 6;
    useStrainFilterMatchedSP = false;
    id = 'Gaussian ellipse AA fit';
    [values,stds,npixels,names,stderrs,nconn] = m4.getStrainInABBufferMask(SPbufferdist,AAremovalradius,id,useStrainFilterMatchedSP,strain_crop_range);
end

if geometry_statistics
    filteramp = 0.3;
    filterrange = 3;
    hotspot_suppression = false;
    edge_threshold = -1;
    gaussian_filter_threshold = 0.05;
    AAblur_Gaussian_sigma = gaussian_filter_threshold;
    AA_amp_threshold = 0.71; % Angstrom
    make_plots = false;
    mask_radius = 10;
    m4.fitAAtoGaussianCircle(filteramp,filterrange,hotspot_suppression,edge_threshold,gaussian_filter_threshold,AA_amp_threshold,make_plots,mask_radius,[],AAblur_Gaussian_sigma);
    
    m4.AA_Gaussian_Circle_Fit_Params(m4.AA_Gaussian_Circle_Fit_Params(:,1)>100 | m4.AA_Gaussian_Circle_Fit_Params(:,1)<0 | m4.AA_Gaussian_Circle_Fit_Params(:,2)>200 | m4.AA_Gaussian_Circle_Fit_Params(:,2)<100,:) = [];
    
    [radius_mean_nm,radius_std_nm,radius_stderr,radius_95confdelta] = m4.getGaussianCircleFitStatistics([-inf,inf]);
end