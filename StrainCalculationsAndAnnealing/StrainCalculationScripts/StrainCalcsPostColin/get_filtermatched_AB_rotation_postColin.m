% get_filtermatched_AA_rotation.m
%
% Goals: first, demonstrate that we can reproduce the existing excel graphs
% using this script. Second, generate new graphs with the re-registered SPs 
% and see if they are any better. Will need a plotting method to validate
% the re-registration.
%
% Nathanael Kazmierczak, 05/25/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload');
dataset_cell = {'05312020_chiralfixed_DS2_postannealing.mat';
                '05312020_chiralfixed_DS4_postannealing.mat';
                '05312020_chiralfixed_DS5_postannealing.mat';
                '05312020_chiralfixed_DS8_postannealing.mat';
                '05312020_chiralfixed_DS10_postannealing.mat';
                '05312020_chiralfixed_DS15_2_postannealing.mat';
                '05312020_chiralfixed_DS18_2_postannealing.mat';
                '05312020_chiralfixed_noprefactor_DS26_postannealing.mat'};
twist_angle_data = [1.0325,1.2276,1.2343,1.3672,0.6272,0.1203,0.1355,0.2609];
twist_angle_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248];
TGVids = {'DS2','DS4','DS5','DS8','DS10','DS15','DS18','DS26'};
dataset_id = [2,4,5,8,10,15,18,26];
use_circles = [4,15,18,26];
use_smaller_filter = [15,18];
n_datasets = numel(dataset_cell);
use_smaller_radius = [4,5,8];
larger_AAradius = 6;
smaller_AAradius = 5;
useStrainFilterMatchedSP = false;
% first column gives the dataset number, second column gives the SP width.
% Do this in accordance with the previous studies in the 2D strain
% calculation workbook.
% SP_bufferwidth_dict = [8,3;4,5;5,5;2,6;10,8;15,11;18,11;26,11];
SP_bufferwidth_dict = [8,3;4,5;5,5;2,6;10,9;15,11;18,11;26,11];
% The SP widths for the last three datasets have been increased by 3 nm.
% There should be room to spare and it should be a bit more accurate.
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct_heavier = buildFilterStruct(stringcell);
stringcell = {'median outlier','both displacements',[3,1]};
pre_filterstruct_lighter = buildFilterStruct(stringcell);

nMat = 7;
all_storage_avs = zeros(n_datasets,nMat);  % six quantitites to store for now
all_storage_stds = zeros(n_datasets,nMat);
all_storage_stderrs = zeros(n_datasets,nMat);
n_pixels_storage = zeros(n_datasets,1);
for i = 1:n_datasets
    load(dataset_cell{i});
    dataset_number = dataset_id(i);
    if nnz(dataset_number == use_circles) > 0
        id = 'Mask circle AA fit';
    else
        id = 'Gaussian ellipse AA fit';
    end
    if nnz(dataset_number == use_smaller_radius) > 0
        AAremovalradius = smaller_AAradius;
    else
        AAremovalradius = larger_AAradius;
    end
    if nnz(dataset_number == use_smaller_filter) > 0
        pre_filterstruct = pre_filterstruct_lighter;
    else
        pre_filterstruct = pre_filterstruct_heavier;
    end
    idx = find(SP_bufferwidth_dict(:,1) == dataset_number);
    SPbufferdist = SP_bufferwidth_dict(idx,2);
    
    
    % Run the strain calculation
    TGVid = TGVids{i};  % versus Colin's original parameter settings.
    SPdirection = 1;
%     colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = 0;
    trim_value = [0,35];  % This is needed for DS5, so do for all.
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,SPdirection,divide_by_two_for_intralayer,sign_convention,trim_tear_value);

    
    
    [values,stds,npixels,names,stderrs,nconn] = m4.getStrainInABBufferMask(SPbufferdist,AAremovalradius,id,useStrainFilterMatchedSP);
%     [averages_mat,stds_mat,n_values,names] = m4.getStrainInAARingedBufferMasks(id,'circular',AA_radius,false);
    
    
    all_storage_avs(i,:) = values';
    all_storage_stds(i,:) = stds';
    all_storage_stderrs(i,:) = stderrs';
    n_pixels_storage(i) = npixels;
    close all
    
end

% Remember, this is now fixed-body rotation of the AB domains, not the AA
% domains.
fixed_body_rotation_avs = all_storage_avs(:,6);
fixed_body_rotation_stds = all_storage_stds(:,6);
fixed_body_rotation_stderrs = all_storage_stderrs(:,6);

figure;
yneg = fixed_body_rotation_stderrs';
ypos = fixed_body_rotation_stderrs';
xneg = twist_angle_stds;
xpos = twist_angle_stds;
errorbar(twist_angle_data,fixed_body_rotation_avs',yneg,ypos,xneg,xpos,'o');
hold on
title('Fixed body rotation in AB domains');
xlabel('Moire twist angle (degrees)');
ylabel('Twist angle (degrees)');
legend('AB reconstruction rotation (1 layer)');

set(gca,'FontSize',14)

