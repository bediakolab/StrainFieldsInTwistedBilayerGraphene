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
dataset_id = [2,4,5,8,10,15,18,26];
use_circles = [4,15,18,26];
use_smaller_filter = [15,18];
n_datasets = numel(dataset_cell);
AA_radius = 1;
useStrainFilterMatchedSP = false;

all_storage_avs = zeros(n_datasets,6);  % six quantitites to store for now
all_storage_stds = zeros(n_datasets,6);  
n_pixels_storage = zeros(n_datasets,1);
for i = 1:n_datasets
    load(dataset_cell{i});
    dataset_number = dataset_id(i);
    if nnz(dataset_number == use_circles) > 0
        id = 'Mask circle AA fit';
    else
        id = 'Gaussian ellipse AA fit';
    end
    if nnz(dataset_number == use_smaller_filter) > 0
        filterstruct = buildFilterStructPreset( 'DS15 annealed displacement #1' );
    else
        filterstruct = buildFilterStructPreset( 'DS26 annealed displacement #2' );
    end
    
    [averages_mat,stds_mat,n_values,names] = m4.getStrainInAARingedBufferMasks(id,'circular',AA_radius,false);
    all_storage_avs(i,:) = averages_mat';
    all_storage_stds(i,:) = stds_mat';
    n_pixels_storage(i) = n_values;
    close all
    
end

fixed_body_rotation_avs = all_storage_avs(:,6);
fixed_body_rotation_stds = all_storage_stds(:,6);
total_twist_angle = twist_angle_data' + fixed_body_rotation_avs;

figure;
yneg = fixed_body_rotation_stds';
ypos = fixed_body_rotation_stds';
xneg = twist_angle_stds;
xpos = twist_angle_stds;
errorbar(twist_angle_data,fixed_body_rotation_avs',yneg,ypos,xneg,xpos,'o');
hold on
plot(twist_angle_data,twist_angle_data,'-');
errorbar(twist_angle_data,total_twist_angle',yneg,ypos,xneg,xpos,'o');
title('Fixed body rotation in AA domains');
xlabel('Moire twist angle (degrees)');
ylabel('Twist angle (degrees)');
legend('Reconstruction rotation (net, over two layers)','Moire rotation','Total rotation');


% title('Fixed body rotation in AA domains');
% xlabel('Twist angle (degrees)');
% ylabel('Intralayer reconstruction fixed-body rotation');
% 
% 
% 
% 
% figure;
% yneg = fixed_body_rotation_stds;
% ypos = fixed_body_rotation_stds;
% xneg = theta_stds;
% xpos = theta_stds;
% errorbar(moire_angles,fixed_body_rotation_means,yneg,ypos,xneg,xpos,'o');
% hold on
% plot(moire_angles,moire_angles,'-');
% errorbar(moire_angles,total_twist_angle,yneg,ypos,xneg,xpos,'o');
% title('Fixed body rotation in AA domains');
% xlabel('Moire twist angle (degrees)');
% ylabel('Twist angle (degrees)');


