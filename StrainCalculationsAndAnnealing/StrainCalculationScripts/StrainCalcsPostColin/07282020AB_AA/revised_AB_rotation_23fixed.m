% revised_AB_rotation_23fixed.m
%
% Goals: first, demonstrate that we can reproduce the existing excel graphs
% using this script. Second, generate new graphs with the re-registered SPs 
% and see if they are any better. Will need a plotting method to validate
% the re-registration.
%
% Nathanael Kazmierczak, 05/25/2020
% Revised for final analysis 07/28/2020

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

% addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload');
% dataset_cell = {'05312020_chiralfixed_DS2_postannealing.mat';
%                 '05312020_chiralfixed_DS4_postannealing.mat';
%                 '05312020_chiralfixed_DS5_postannealing.mat';
%                 '05312020_chiralfixed_DS8_postannealing.mat';
%                 '05312020_chiralfixed_DS10_postannealing.mat';
%                 '05312020_chiralfixed_DS15_2_postannealing.mat';
%                 '05312020_chiralfixed_DS18_2_postannealing.mat';
%                 '05312020_chiralfixed_noprefactor_DS26_postannealing.mat'};
% twist_angle_data = [1.0325,1.2276,1.2343,1.3672,0.6272,0.1203,0.1355,0.2609];
% twist_angle_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248];
% TGVids = {'DS2','DS4','DS5','DS8','DS10','DS15','DS18','DS26'};
% dataset_id = [2,4,5,8,10,15,18,26];


use_circles = [4,15,18,26];
use_smaller_filter = [15,18];
n_datasets = numel(dataset_cell);
use_smaller_radius = [4,5,8,23.2,24.2,25.2,3.2];
larger_AAradius = 6;
smaller_AAradius = 5;
useStrainFilterMatchedSP = false;
% first column gives the dataset number, second column gives the SP width.
% Do this in accordance with the previous studies in the 2D strain
% calculation workbook.
% SP_bufferwidth_dict = [8,3;4,5;5,5;2,6;10,8;15,11;18,11;26,11];
SP_bufferwidth_dict = [8,3;
                       4,5;
                       5,5;
                       2,6;
                       10,9;
                       15,11;
                       18,11;
                       26,11;
                       30.2,9;
                       8.2,9;
                       9.2,9;
                       10.2,9;
                       11.2,9;
                       12.2,9;
                       23.2,3;
                       24.2,3;
                       25.2,3;
                       3.2,5;
                       15.2,11;
                       1.2,11];
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
    SPtensordirection = 1;
%     colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = false;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = 0;
     if any(use_smaller_trim == dataset_number)
        trim_value = [0,10];
    else
        trim_value = [0,35];  
     end
    calib_method = 'SP';
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,calib_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,SPtensordirection);

    
    
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
yneg = 2*fixed_body_rotation_stderrs';
ypos = 2*fixed_body_rotation_stderrs';
xneg = twist_angle_stds;
xpos = twist_angle_stds;
errorbar(twist_angle_data,fixed_body_rotation_avs',yneg,ypos,xneg,xpos,'o');
hold on
title('Fixed body rotation in AB domains');
xlabel('Moire twist angle (degrees)');
ylabel('Twist angle (degrees)');
legend('AB reconstruction rotation (2 layers)');

set(gca,'FontSize',14)

sortmat = [twist_angle_data',fixed_body_rotation_avs,yneg'];
sorted = sortrows(sortmat);
