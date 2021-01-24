% SPstrain_SP1_newestdata_23fixed.m
%
% This series of scripts will be for the final paper, ideally. 
%
% Nathanael Kazmierczak, 07/27/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2');
dataset_cell = {'05312020_chiralfixed_DS2_postannealing.mat';
                '05312020_chiralfixed_DS4_postannealing.mat';
                '05312020_chiralfixed_DS5_postannealing.mat';
                '05312020_chiralfixed_DS8_postannealing.mat';
                '05312020_chiralfixed_DS10_postannealing.mat';
                '05312020_chiralfixed_DS15_2_postannealing.mat';
                '05312020_chiralfixed_DS18_2_postannealing.mat';
                '05312020_chiralfixed_noprefactor_DS26_postannealing.mat';
                '06222020_DS30S2_annealed.mat';  % New datasets from here on below
                '06222020_DS8S2_annealed.mat';
                '06222020_DS9S2_annealed.mat';
                '06222020_DS10S2_annealed.mat';
                '06222020_DS11S2_annealed.mat';
                '06222020_DS12S2_annealed.mat';
                '06222020_DS23S2_annealed.mat';
                '06222020_DS24S2_annealed.mat';
                '06222020_DS25S2_annealed.mat';
                '06222020_DS3S2_annealed.mat';
                '06222020_DS15S2_annealed.mat';
                '06222020_DS1S2_annealed.mat'};
twist_angle_data = [1.0325,1.2276,1.2343,1.3672,0.6272,0.1203,0.1355,0.2609,0.8442,0.633,0.6632,0.6385,0.6543,0.7534,1.3067,1.3178,1.3468,1.1669,0.16,0.2894];
twist_angle_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248,0.0547,0.0239,0.0229,0.0053,0.0096,0.0299,0.0317,0.0472,0.0344,0.0444,0.03,0.0066];
dataset_id = [2,4,5,8,10,15,18,26,30.2,8.2,9.2,10.2,11.2,12.2,23.2,24.2,25.2,3.2,15.2,1.2];
TGVids = {'DS2','DS4','DS5','DS8','DS10','DS15quant','DS18quant','DS26','DS30S2','DS8S2','DS9S2','DS10S2','DS11S2','DS12S2','DS8','DS8','DS8','DS3S2','DS15S2','DS1S2'};
use_smaller_trim = [8.2,9.2,10.2,11.2,12.2,23.2,24.2,25.2,15.2,1.2];

            
        
flip_signs = [];
use_circles = [4,15,18,26];
use_smaller_radius = [4,5,8,23.2,24.2,25.2,3.2];
use_smaller_filter = [15,18];
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct_heavier = buildFilterStruct(stringcell);
stringcell = {'median outlier','both displacements',[3,1]};
pre_filterstruct_lighter = buildFilterStruct(stringcell);

n_datasets = numel(dataset_cell);
larger_radius = 6;
smaller_radius = 5;
useStrainFilterMatchedSP = false;

eyx_storage = zeros(n_datasets,2);
exy_storage = zeros(n_datasets,2);

nMats = 8;  % number of different matrices for which we are computing values. This could change.
vals_storage = zeros(n_datasets,nMats);
stds_storage = zeros(n_datasets,nMats);
stderrs_storage = zeros(n_datasets,nMats);
nconn_storage = zeros(n_datasets,1);

for i = 1:n_datasets
    load(dataset_cell{i});
    dataset_number = dataset_id(i); 
    
    if any(use_smaller_trim == dataset_number)
        trim_value = [0,10];
    else
        trim_value = [0,35];  
    end
    if nnz(dataset_number == use_circles) > 0
        id = 'Mask circle AA fit';
    else
        id = 'Gaussian ellipse AA fit';
    end
    if nnz(dataset_number == use_smaller_radius) > 0
        AAremovalradius = smaller_radius;
    else
        AAremovalradius = larger_radius;
    end
    if nnz(dataset_number == use_smaller_filter) > 0
        pre_filterstruct = pre_filterstruct_lighter;
    else
        pre_filterstruct = pre_filterstruct_heavier;
    end
    
    % Run the strain calculation
    TGVid = TGVids{i};  % versus Colin's original parameter settings.
%     SPdirection = 1;
    calibration_method = 'SP';
    tensor_rotation_direction = 1;
%     colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = 0;
%     trim_value = [0,35];  % Txhis is needed for DS5, so do for all.
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,calibration_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotation_direction);

    
    
%     AA_buffer_radius = 2;
%     optional_SP_tolerance = 0.6;
%     m4.setStrainFilterMatchedSPRegistration(filterstruct,id,AA_buffer_radius,optional_SP_tolerance);
%     The above settings were used to assign all of the SP registrations,
%     by hand to verify correctness.
    
%     SPdirection = 'all';
    SPdirection = 1;
    SPbufferdist = 1; % nm
%     useStrainFilterMatchedSP = 0;
    
    [values,stds,npixels,names,stderrs,nconn] = m4.getStrainInSPBufferMask(SPdirection,SPbufferdist,...
        AAremovalradius,id,useStrainFilterMatchedSP);
    
    vals_storage(i,:) = values';
    stds_storage(i,:) = stds';
    stderrs_storage(i,:) = stderrs';
    nconn_storage(i) = nconn;
    
    if nnz(dataset_number == flip_signs) > 0
        exy_storage(i,:) = [-values(2),stds(2)];
        eyx_storage(i,:) = [-values(3),stds(3)];
    else
        exy_storage(i,:) = [values(2),stds(2)];
        eyx_storage(i,:) = [values(3),stds(3)];
    end
    close all
    
end

% figure;
% x = 1:10:100;
% y = [20 30 45 40 60 65 80 75 95 90];
% yneg = eyx_storage(:,2)';
% ypos = eyx_storage(:,2)';
% xneg = twist_angle_stds;
% xpos = twist_angle_stds;
% errorbar(twist_angle_data,eyx_storage(:,1)',yneg,ypos,xneg,xpos,'o');
% title('eyx simple shear perpendicular to soliton wall 1');
% xlabel('Twist angle (degrees)');
% ylabel('Strain %');
% current_eyx_storage = eyx_storage;

% hold on
% load('SPstrain_05262020_orgSPs.mat');
% % eyx_storage is now the previous calculation values
% yneg = eyx_storage(:,2)';
% ypos = eyx_storage(:,2)';
% errorbar(twist_angle_data,eyx_storage(:,1)',yneg,ypos,xneg,xpos,'o');
% legend('New SP registration','Original SP registration');

%% Make a plot of fixed body rotation in SP domains!
% Sort before extraction so that we can make lines
vals_aug = horzcat(twist_angle_data',vals_storage);
stderrs_aug = horzcat(twist_angle_data',twist_angle_stds',stderrs_storage);
vals_sort = sortrows(vals_aug);
stderrs_sort = sortrows(stderrs_aug);
vals_use = vals_sort(:,2:end);
stderrs_use = stderrs_sort(:,3:end);
twist_angle_data_use = vals_sort(:,1)';
twist_angle_stderrs_use = stderrs_sort(:,2)';

fixed_body_rotation_avs = vals_use(:,6);
fixed_body_rotation_stderrs = stderrs_use(:,6);
sxy_avs = vals_use(:,2);
sxy_stderrs = stderrs_use(:,2);
syx_avs = vals_use(:,3);
syx_stderrs = stderrs_use(:,3);
principal_shear_avs = vals_use(:,8);
principal_shear_stderrs = stderrs_use(:,8);

%%% Original without drawing the lines.
% fixed_body_rotation_avs = vals_storage(:,6);
% fixed_body_rotation_stderrs = stderrs_storage(:,6);
% sxy_avs = vals_storage(:,2);
% sxy_stderrs = stderrs_storage(:,2);
% syx_avs = vals_storage(:,3);
% syx_stderrs = stderrs_storage(:,3);
% principal_shear_avs = vals_storage(:,8);
% principal_shear_stderrs = stderrs_storage(:,8);

% Multiplication arises on account of the "divide by two for intralayer"
% flag.
figure
errorbar(twist_angle_data_use,2*fixed_body_rotation_avs',2*2*fixed_body_rotation_stderrs',2*2*fixed_body_rotation_stderrs',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
title('Fixed-body rotation in SP domains, Stderr from connected components');
xlabel('Twist angle (degrees)');
ylabel('Rotation angle (degrees)');
set(gca,'FontSize',14);
legend('SP rotation (2 layers), 95% CI error bars');



figure
errorbar(twist_angle_data_use,sxy_avs',2*sxy_stderrs',2*sxy_stderrs',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
hold on
errorbar(twist_angle_data_use,syx_avs',2*syx_stderrs',2*syx_stderrs',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
hold on
errorbar(twist_angle_data_use,principal_shear_avs',2*principal_shear_stderrs',2*principal_shear_stderrs',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
title('Shear strain in SP domains, Stderr from connected components');
xlabel('Twist angle (degrees)');
ylabel('Strain %')
set(gca,'FontSize',14);
legend('sxy','syx','Gamma max');

disp('Number of connected components for each dataset id (left column) displayed below (right column)');
[dataset_id',nconn_storage]


