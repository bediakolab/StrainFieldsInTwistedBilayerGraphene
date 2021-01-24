% get_filtermatched_SP_strain.m
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
flip_signs = [];
use_circles = [4,15,18,26];
use_smaller_radius = [4,5,8];
use_smaller_filter = [15,18];
n_datasets = numel(dataset_cell);
larger_radius = 6;
smaller_radius = 5;
useStrainFilterMatchedSP = false;

eyx_storage = zeros(n_datasets,2);
exy_storage = zeros(n_datasets,2);

vals_storage = zeros(n_datasets,6);
stds_storage = zeros(n_datasets,6);
for i = 1:n_datasets
    load(dataset_cell{i});
    dataset_number = dataset_id(i); 
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
        filterstruct = buildFilterStructPreset( 'DS15 annealed displacement #1' );
    else
        filterstruct = buildFilterStructPreset( 'DS26 annealed displacement #2' );
    end
    
%     AA_buffer_radius = 2;
%     optional_SP_tolerance = 0.6;
%     m4.setStrainFilterMatchedSPRegistration(filterstruct,id,AA_buffer_radius,optional_SP_tolerance);
%     The above settings were used to assign all of the SP registrations,
%     by hand to verify correctness.
    
    SPdirection = 1;
    SPbufferdist = 1; % nm
%     useStrainFilterMatchedSP = 0;
    
    [values,stds,npixels,names] = m4.getStrainInSPBufferMask(SPdirection,SPbufferdist,...
        AAremovalradius,id,useStrainFilterMatchedSP);
    
    vals_storage(i,:) = values';
    stds_storage(i,:) = stds';
    
    if nnz(dataset_number == flip_signs) > 0
        exy_storage(i,:) = [-values(2),stds(2)];
        eyx_storage(i,:) = [-values(3),stds(3)];
    else
        exy_storage(i,:) = [values(2),stds(2)];
        eyx_storage(i,:) = [values(3),stds(3)];
    end
    close all
    
end

figure;
x = 1:10:100;
y = [20 30 45 40 60 65 80 75 95 90];
yneg = eyx_storage(:,2)';
ypos = eyx_storage(:,2)';
xneg = twist_angle_stds;
xpos = twist_angle_stds;
errorbar(twist_angle_data,eyx_storage(:,1)',yneg,ypos,xneg,xpos,'o');
title('eyx simple shear perpendicular to soliton wall 1');
xlabel('Twist angle (degrees)');
ylabel('Strain %');
current_eyx_storage = eyx_storage;

% hold on
% load('SPstrain_05262020_orgSPs.mat');
% % eyx_storage is now the previous calculation values
% yneg = eyx_storage(:,2)';
% ypos = eyx_storage(:,2)';
% errorbar(twist_angle_data,eyx_storage(:,1)',yneg,ypos,xneg,xpos,'o');
% legend('New SP registration','Original SP registration');

%% Make a plot of fixed body rotation in SP domains!
fixed_body_rotation_avs = vals_storage(:,6);
fixed_body_rotation_stds = stds_storage(:,6);
pure_shear_strain_avs = vals_storage(:,5);
pure_shear_strain_stds = stds_storage(:,5);
simple_shear_strain_avs = vals_storage(:,3);
simple_shear_strain_stds = stds_storage(:,3);

figure
errorbar(twist_angle_data,fixed_body_rotation_avs',fixed_body_rotation_stds',fixed_body_rotation_stds',xneg,xpos,'o');
title('Fixed-body rotation in SP domains');
xlabel('Twist angle (degrees)');
ylabel('Rotation angle (degrees)');
set(gca,'FontSize',14);
legend('Net SP rotation (2 layers)');

% If this becomes a narrative, eventually combine these onto one plot with 
% the above, using two y-axes.
figure
errorbar(twist_angle_data,pure_shear_strain_avs',pure_shear_strain_stds',pure_shear_strain_stds',xneg,xpos,'o');
hold on
errorbar(twist_angle_data,simple_shear_strain_avs',simple_shear_strain_stds',simple_shear_strain_stds',xneg,xpos,'o');
title('Shear strain in SP domains');
xlabel('Twist angle (degrees)');
ylabel('Strain %')
set(gca,'FontSize',14);
legend('Pure shear strain (net, 2 layers)','Simple shear strain (net, 2 layers)');



