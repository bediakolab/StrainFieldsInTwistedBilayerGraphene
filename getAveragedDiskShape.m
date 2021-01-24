function [average_probe_stacking, average_probe_convolution] = getAveragedDiskShape(DPs,disk_positions_cell,disk_indices_cell,index_choice,window_width,subpixel)
% getAveragedDiskShape.m
%
% Function for obtaining the average probe from auto-detected disks in
% diffraction pattern. 
%
% Externally, disks are registered through correlation to generate
% bragg_disk_locations. This should take place for maybe 50-100 disks at
% least.
%
% Then, we will convolve the disks against each other. Convolving the disks
% one at a time, align via the maximum inner product which must fall on a
% pixel, essentially using cross-correlation. But we could also use phase
% correlation, e.g.
% 
% DPs should be preprocessed so that they are a 3D array, with all real
% space information compressed into the third dimension.
%
% Note it is essential that disk_positions and disk_indices line up with
% each other row-wise. Have to ensure this in caller.

% First, need to cut away ROIs around the given disks of interest. Build a
% function for this.
% The reason that the cell is needed rather than a traditional array is
% that sometimes we loose a disk in the registration process; thus this
% might be four rather than five, for instance. 
%
% We will make the roicuts a little bit thicker after the fact. This is so
% that alignment will be easier once we go to actually do the alignment
% later in the function. 


num_DPs = size(DPs,3);
% fliplr here for the sake of the testing driver; note that in test cut
% roi this was also done originally.
[ roicut1, Irangecut1, Jrangecut1 ] = cutROI( DPs(:,:,1), fliplr(disk_positions_cell{1}), disk_indices_cell{1}, ...
    index_choice, window_width );
DPROIs = zeros(size(roicut1,1)*3,size(roicut1,2)*3,num_DPs);
pad = zeros(size(roicut1));
for i = 1:num_DPs
    [ roicut, Irangecut, Jrangecut ] = cutROI( DPs(:,:,i), fliplr(disk_positions_cell{i}), disk_indices_cell{i}, ...
        index_choice, window_width );
    DPROIs(:,:,i) = [pad,pad,pad;
                     pad,roicut,pad;
                     pad,pad,pad];
end

[ average_probe_convolution ] = getAveragedDiskShape_innards( DPROIs );


end


