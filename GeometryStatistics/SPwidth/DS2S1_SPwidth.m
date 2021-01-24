% DS2S1_SPwidth.m
%
% Nathanael Kazmierczak, 07/09/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
    load('05312020_chiralfixed_DS2_postannealing.mat');
end

filterstring = {'soft multistart','amplitude',[0.3,5]};
filterstruct = buildFilterStruct(filterstring);
cut_number = 10;  % along the base_t direction
cut_spacing_nm = 0.1;
AA_center_type = 'Gaussian ellipse AA fit';
radius_value = 5;
image_boundary_buffer = 25;  % pixels
pixel_size_cutoff = 5;
window_dims_nm = [3,5];  % First value is for the width, second for the length
cutoff_score = 0.5;
[distcoords,av_scores,distance,endpoint1,endpoint2,rangemean,rangemeanstderr] = ...
    m4.getSPwidths(filterstruct,cut_number,cut_spacing_nm,AA_center_type,radius_value,image_boundary_buffer,pixel_size_cutoff,window_dims_nm,cutoff_score)


