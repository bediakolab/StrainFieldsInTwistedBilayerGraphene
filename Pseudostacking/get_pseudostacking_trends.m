% get_pseudostacking_trends.m
%
% This function iteratively loads all datasets in a list, applies the TGV
% filter used for strain mapping, calculates the pseudostacking area % for
% each, and makes a plot.
%
% Nathanael Kazmierczak, 06/20/2020

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

ps_storage = zeros(n_datasets,3);  % for the three pseudostacked regions.
% all_storage_stds = zeros(n_datasets,nMat);
% all_storage_stderrs = zeros(n_datasets,nMat);
% n_pixels_storage = zeros(n_datasets,1);
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
    calibration_method = 'SP';
%     colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = 0;
    trim_value = [0,35];  % This is needed for DS5, so do for all.
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,calibration_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,SPdirection);

    use_annealed = 2; % New setting, makes other things unused because references operations by computeStrainMaps3.
    croprange = [];  % unused
    filterstruct = [];
    r = 0.3;
    make_plots = true;
    [AA_perc,AB_perc,SP_perc,pseudostacking_mat] = m4.getPseudostackingAssignments(make_plots,r,filterstruct,use_annealed,croprange);
%     [values,stds,npixels,names,stderrs,nconn] = m4.getStrainInABBufferMask(SPbufferdist,AAremovalradius,id,useStrainFilterMatchedSP);
%     [averages_mat,stds_mat,n_values,names] = m4.getStrainInAARingedBufferMasks(id,'circular',AA_radius,false);
    
    
    ps_storage(i,:) = [AA_perc,AB_perc,SP_perc];
%     all_storage_stds(i,:) = stds';
%     all_storage_stderrs(i,:) = stderrs';
%     n_pixels_storage(i) = npixels;

    % Store the pretty figures
    save_figures = true;
    if save_figures
        h = figure(3);
        h.HandleVisibility = 'off';
        close all
        h.HandleVisibility = 'on';
        filestring = sprintf('Pseudostacking_DS_%d',dataset_id(i));
        filestring2 = sprintf('Pseudostacking_DS_%d.png',dataset_id(i));
        savefig(h,filestring);
        saveas(h,filestring2);
        excel_name = sprintf('Pseudostacking_DS_%d.csv',dataset_id(i));
        xlswrite(excel_name,pseudostacking_mat);
    end
    close all
end


tosort = [twist_angle_data',ps_storage];
sorted = sortrows(tosort);

figure;
hold on
% lay down the regime shading first
alpha = 0.1;
ylimvals = [0,80];
xpatch = [0,0,0.6,0.6];
ypatch = [ylimvals(1),ylimvals(2),ylimvals(2),ylimvals(1)];
p = patch(xpatch,ypatch,'m');
p.FaceAlpha = alpha;
xpatch = [0.6,0.6,1.4,1.4];
ypatch = [ylimvals(1),ylimvals(2),ylimvals(2),ylimvals(1)];
p = patch(xpatch,ypatch,'c');
p.FaceAlpha = alpha;

% now plot data
plot(sorted(:,1),sorted(:,2),'-o');
hold on
plot(sorted(:,1),sorted(:,3),'-o');
plot(sorted(:,1),sorted(:,4),'-o');
title('r = 0.3 Angstrom pseudostacking');
xlabel('Moire twist angle (degrees)');
ylabel('Pseudostacking %');
legend('AB reconstruction regime','AA reconstruction regime','AA','AB','SP');
set(gca,'FontSize',14)

