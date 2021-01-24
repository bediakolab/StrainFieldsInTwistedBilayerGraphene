% get_pseudostacking_trends_newdatacrit_minimalTGV.m
%
% This function iteratively loads all datasets in a list, applies the TGV
% filter used for strain mapping, calculates the pseudostacking area % for
% each, and makes a plot.
%
% This script has included the data from session 2, as well as removing
% some pointless variables. 
%
% Nathanael Kazmierczak, 07/12/2020

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
                '06222020_DS3S2_annealed.mat'};
% twist_angle_data = [1.0325,1.2276,1.2343,1.3672,0.6272,0.1203,0.1355,0.2609];
% twist_angle_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248];
% TGVids = {'DS2','DS4','DS5','DS8','DS10','DS15','DS18','DS26'};
twist_angle_data = [1.0325,1.2276,1.2343,1.3672,0.6272,0.1203,0.1355,0.2609,0.8442,0.633,0.6632,0.6385,0.6543,0.7534,1.3067,1.3178,1.3468,1.1669];
twist_angle_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248,0.0547,0.0239,0.0229,0.0053,0.0096,0.0299,0.0317,0.0472,0.0344,0.0444];
dataset_id = [2,4,5,8,10,15,18,26,30.2,8.2,9.2,10.2,11.2,12.2,23.2,24.2,25.2,3.2];
TGVid = 'minimal';
% TGVids = {'DS2','DS4','DS5','DS8','DS10','DS15','DS18','DS26','DS30S2','DS8S2','DS9S2','DS10S2','DS11S2','DS12S2','DS8S2','DS8S2','DS8S2','DS3S2'};
use_smaller_trim = [8.2,9.2,10.2,11.2,12.2,23.2,24.2,25.2];

% dataset_id = [2,4,5,8,10,15,18,26];
% % % use_circles = [4,15,18,26];
use_smaller_filter = [15,18];
n_datasets = numel(dataset_cell);
% % % use_smaller_radius = [4,5,8];
% % % larger_AAradius = 6;
% % % smaller_AAradius = 5;
% % % useStrainFilterMatchedSP = false;
% first column gives the dataset number, second column gives the SP width.
% Do this in accordance with the previous studies in the 2D strain
% calculation workbook.
% SP_bufferwidth_dict = [8,3;4,5;5,5;2,6;10,8;15,11;18,11;26,11];
% % % SP_bufferwidth_dict = [8,3;4,5;5,5;2,6;10,9;15,11;18,11;26,11];
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
    fprintf('Analyzing dataset %d of %d.\n',i,n_datasets);
    load(dataset_cell{i});
    dataset_number = dataset_id(i);
% % %     if nnz(dataset_number == use_circles) > 0
% % %         id = 'Mask circle AA fit';
% % %     else
% % %         id = 'Gaussian ellipse AA fit';
% % %     end
%     if nnz(dataset_number == use_smaller_radius) > 0
%         AAremovalradius = smaller_AAradius;
%     else
%         AAremovalradius = larger_AAradius;
%     end
    if nnz(dataset_number == use_smaller_filter) > 0
        pre_filterstruct = pre_filterstruct_lighter;
    else
        pre_filterstruct = pre_filterstruct_heavier;
    end
% % %     idx = find(SP_bufferwidth_dict(:,1) == dataset_number);
% % %     SPbufferdist = SP_bufferwidth_dict(idx,2);
    
    
    % Run the strain calculation
%     TGVid = TGVids{i};  % versus Colin's original parameter settings.
    SPdirection = 1;
    calibration_method = 'SP';
%     colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
    divide_by_two_for_intralayer = true;
    sign_convention = [+1,-1];
    % sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
    trim_tear_value = 0;
    if any(use_smaller_trim == dataset_number)
        trim_value = [0,10];
    else
        trim_value = [0,35];  % This is needed for DS5, so do for all.
    end
    m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,calibration_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,SPdirection);

    use_annealed = 2; % New setting, makes other things unused because references operations by computeStrainMaps3.
    croprange = [];  % unused
    filterstruct = [];
%     r = 0.3;
    r = 'modified wigner-seitz';
    make_plots = true;
    [AA_perc,AB_perc,SP_perc,pseudostacking_mat] = m4.getPseudostackingAssignments(make_plots,r,filterstruct,use_annealed,croprange);

    ps_storage(i,:) = [AA_perc,AB_perc,SP_perc];

    % Store the pretty figures
    save_figures = true;
    if save_figures
        h = figure(3);
        h.HandleVisibility = 'off';
        close all
        h.HandleVisibility = 'on';
        if mod(dataset_id(i),1) ~= 0
            dsid_to_use = round(dataset_id(i));
            filestring = sprintf('Pseudostacking_DS_%dS2',dsid_to_use);
            filestring2 = sprintf('Pseudostacking_DS_%dS2.png',dsid_to_use);
        else
            filestring = sprintf('Pseudostacking_DS_%d',dataset_id(i));
            filestring2 = sprintf('Pseudostacking_DS_%d.png',dataset_id(i));
        end
        savefig(h,filestring);
        saveas(h,filestring2);
        excel_name = sprintf('Pseudostacking_DS_%d.csv',dataset_id(i));
        xlswrite(excel_name,pseudostacking_mat);
    end
    close all
end


tosort = [twist_angle_data',ps_storage];
sorted = sortrows(tosort);

tosortstds = [twist_angle_data',twist_angle_stds'];
sorted_stds = sortrows(tosortstds);

figure;
hold on
% lay down the regime shading first
alpha = 0.1;
ylimvals = [0,90];
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

