% track_interference_motion.m
%
% Script for tracking the how the interference patterns move as the real
% space probe moves in MV_1.5_9 data collection, particularly dataset 2.
%
% Nathanael Kazmierczak, 09/03/2019
% Bediako Lab, UC Berkeley.


addpath('/Volumes/SANDISK/BediakoLabLexar');
addpath('/Volumes/SANDISK/BediakoLabLexar/4DSTEMdata_updating');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM');
format long

GENERATE_NEW_DATA = 0;   % set to zero if we want to use previously defined 
% roi cuts and data saves, to speed debugging.

filename = '2_MV_1.5_9_20x20_ss3nm_1s_spot 9_alpha=1_bin1_cl=130_60kV__NPK_nobin.h5';
if GENERATE_NEW_DATA
%% (0) Load the top 1/10th horizontal slice of the diffraction pattern.
% This should be managable in RAM.


% filename = '3_MV_1.5_9_20x20_ss3nm_2s_spot 7_alpha=1_bin1_cl=130_60kV_NPK_nobin.h5';
info = h5info(filename);
info.Groups.Groups(1).Groups(1).Groups.Datasets(1).Dataspace  % looks like it is 70 x 8

start = [1,1,1,1];
count = [2048,2048,2,20];  % chunks are in the fourth dimension
this_data = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
firstDP = this_data(:,:,1,1);
figure; image(uint8(double(firstDP).^0.5),'CDataMapping','Scaled'); colormap(hsv); colorbar; pbaspect([1 1 1]);

%% (1) Cut down the image to a single disk of interest
disp('Select the location on the graph for the roi center.');
[central_position] = fliplr(ginput(1));
disp('Please select a horizontal location for the perimeter of the roi');
[pos2] = fliplr(ginput(1));
% Take the horizontal distance as what matters.
window_width = int16(2*abs(central_position(2) - pos2(2)));

[ roicut, IRangeCut, JRangeCut ] = cutROI_fromPosition( firstDP, central_position, window_width );
else
    loadfile = '09032019_intereference_pattern_tracking_ROIcut.mat';
    load(loadfile,'roicut','IRangeCut','JRangeCut');
    load('09032019_first_DP_save_interference_tracking.mat','firstDP');
end

%% (2) Draw a polygonal ROI around the approximate interference rod region
% % 
% % disp('Please draw the polygonal ROI region that will be used as the preliminary probe kernel.');
% % BW = roipoly();
% % % Verify the ROI has been constructed correctly
% % figure; 
% % masked = roicut;
% % masked(~BW) = 0;
% % image(uint8(masked),'CDataMapping','Scaled'); colormap(hsv); colorbar; pbaspect([1 1 1]);
flat_kernel = 0;
boundary_kernel = 0;
[ kernel ] = getROIKernel( roicut, flat_kernel, boundary_kernel );

% [ roicut ] = cutROI_fromRanges( IRangeCut, JRangeCut )

% Refactored into register_interference_rods so that we can be consistent
% with the testing on 1 DP and what goes down on all of them.
%% (3) Register polygonal rods via cross correlation
% Only taking the arguments which pertain to pixel disk registration.
sigma = 0;   % all of these setting are applied to all DPs.
subpixel = 'False';
% % % [subpixel_mask,conv_map,~,~,~,~] = braggDiskRegistrationRecipe_innards(roicut,kernel,sigma,subpixel);
% going into the thresholding, the typical format is [intensity, Icoord?,
% Jcoord?, window id]
% % % [ intensities, xinds, yinds, maxima_storage ] = convMax( subpixel_mask, conv_map );

%% (4) Filter for reliable registrations
% New function built for thresholding 
rel_intensity_threshold = 0.7;  % will need to tune this parameter.
min_disk_distance = 20;

%%%% New code after refactoring (same as below for all DPs).
[ sortedmaxcombo ] = registerInterferenceRods( roicut,kernel,sigma,subpixel, rel_intensity_threshold, min_disk_distance );
% % % [ sortedmaxcombo ] = relIntensityThreshold( maxima_storage, rel_intensity_threshold );
% % % [sortedmaxcombo] = removeTooNear(sortedmaxcombo,min_disk_distance);
% visualize:
plotDPandDisks( roicut, sortedmaxcombo(:,2:3) );
% For now, do we even need a reliable registration near the edges of the
% Bragg disks? I don't see why we would -- should be able to get the
% appropriate pattern in the middle.

%% (5) Stack on top of each other and average through to get averaged kernel
% roi cuts for the stacking are thickened versions of the probe kernel
thickening_iterations = 9;
[ thickened_mask ] = thickenMask( kernel, thickening_iterations );
figure; image(thickened_mask,'CDataMapping','Scaled'); colormap(pink); colorbar; pbaspect([1 1 1]);

%% (3-5) Register all disks
sizex = 7;
sizey = 7;
roicut_storage = zeros(size(roicut,1),size(roicut,2),sizey,sizex);
registration_storage = cell(sizex,sizey);
roi_masked_cut_storage = zeros(size(roicut,1),size(roicut,2),1);  % This is a manual collapse of the fourth dimension in effect.
% Comes out to about 90 MB -- not too bad for storage in RAM.
first_one = 1;
center = [round(size(roicut,1)/2), round(size(roicut,2)/2)];
for i = 1:sizey
    for j = 1:sizex
        start = [1,1,i,j];
        count = [2048,2048,1,1];
        this_DP = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
        [ this_roicut ] = cutROI_fromRanges( this_DP, IRangeCut, JRangeCut ); % where the ranges had better have been defined before.
        roicut_storage(:,:,i,j) = this_roicut; 
        this_registration_max_storage = registerInterferenceRods( this_roicut,kernel,sigma,subpixel, rel_intensity_threshold, min_disk_distance );
        registration_storage{i,j} = this_registration_max_storage;
        num_registered = size(this_registration_max_storage,1);  
        
        % now we are going to loop through here and eliminate extraneous
        % peaks for all of our registration locations -- then we will get
        % accurate self-convolutions.
        for k = 1:num_registered
            fprintf('%d, %d, %d, (i,j,k)\n',i,j,k);
            [ this_roi_mask_cut ] = cutROI_from_shape_and_position( this_roicut, this_registration_max_storage(k,2:3), thickened_mask );
            % We have performed the roi_mask_cut, but we now want to
            % realign for accurate averaging.
            %
            % Aligning the com of the mask rather than the interference rod
            % is super hacky, but might work since we already have a
            % convolution step in play.
            [ this_com_coords ] = getCOMofMASK( logical(this_roi_mask_cut) );
            shift = center - round(this_com_coords);
            this_roi_mask_cut = shiftImage(this_roi_mask_cut,shift(2),shift(1));
            if first_one
                roi_masked_cut_storage(:,:,1) = this_roi_mask_cut;
                first_one = 0;
                plotDPandDisks( this_roicut, this_registration_max_storage(:,2:3) );
                figure; image(uint8(this_roi_mask_cut),'CDataMapping','Scaled'); colormap(pink); colorbar; pbaspect([1 1 1]); title 'ROI Mask Cut for Convolution Alignment';
            else
                roi_masked_cut_storage(:,:,end+1) = this_roi_mask_cut;  % This is now a 3D array with everything we need to make the accurate probe.
            end
        end
    end
end


[ average_probe_convolution, average_probe_stacking ] = getAveragedDiskShape_noIndices( roi_masked_cut_storage,subpixel,max_convolve );
figure; image(uint8(average_probe_convolution),'CDataMapping','Scaled'); colormap(pink); colorbar; pbaspect([1 1 1]); title 'Average probe convolution';
figure; image(uint8(average_probe_stacking),'CDataMapping','Scaled'); colormap(pink); colorbar; pbaspect([1 1 1]); title 'Average probe stacking';

%% (6) Register with the averaged kernel
% Now that we have the probes from the first approximate cross correlation,
% we will go through again in much the same way
sizeytot = 20;
sizextot = 20;
registration_storage_final = cell(sizeytot,sizextot);
DProi_storage_final = cell(sizeytot,sizextot);
storage_final_array = [];
count = [2048,2048,1,1];
for i = 1:sizeytot
    for j = 1:sizextot
        fprintf('%d %d (i,j)\n',i,j);
        start = [1,1,i,j];
        this_DP = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
        [ this_roicut ] = cutROI_fromRanges( this_DP, IRangeCut, JRangeCut ); % where the ranges had better have been defined before.
        
%         this_roicut = roicut_storage(:,:,i,j);   % use what we have stored about the first ROI cuts origianll
        this_registration_max_storage = registerInterferenceRods( this_roicut,average_probe_stacking,sigma,subpixel, rel_intensity_threshold, min_disk_distance );
        registration_storage_final{i,j} = this_registration_max_storage;
        num_registered = size(this_registration_max_storage,1);
        
        DProi_storage_final{i,j} = this_roicut;
        storage_final_array = vertcat(storage_final_array,[this_registration_max_storage,repmat(i,size(this_registration_max_storage,1),1),repmat(j,size(this_registration_max_storage,1),1)]);
    end
end


%% (7) Trace the distances 

% As we move across i = 2, how does the x coord (i.e. third column) change
% with j?
to_sort = storage_final_array;
to_sort(to_sort(:,4) ~= 2,:) = [];
figure;
scatter(to_sort(:,5),to_sort(:,3),'filled');
title 'I = 2 pinned, allowing J (the column in real space) to move'
xlabel('Column # in real space');
ylabel('X-Pixels of registered interference rods');

% As we move across j = 5, how does the y coord (i.e., second column)
% chnage with i?
to_sort_2 = storage_final_array;
to_sort_2(to_sort_2(:,5) ~= 5,:) = [];
figure;
scatter(to_sort_2(:,4),to_sort_2(:,2),'filled');
title 'J = 5 pinned, allowing I (the row in real space) to move'
xlabel('Row # in real space');
ylabel('Y-Pixels of registered interference rods');


%% (8) Validate the disk registration


plotAllRods( DProi_storage_final,storage_final_array,2,6 )
plotAllRods( DProi_storage_final,storage_final_array,2,17 )
plotAllRods( DProi_storage_final,storage_final_array,6,5 )
plotAllRods( DProi_storage_final,storage_final_array,13,5 )
