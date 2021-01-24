function [false_subpixel_mask,conv_map,intensities_subpixel,coords_subpixel,intensities_pixel,...
    coords_pixel] = braggDiskRegistrationRecipe_innards(DP,probe_kernel,sigma,subpixel)
% Mostly refactored code from braggDiskRegistrationRecipe_1.m
%
% This code has the goal of performing the disk registration.
%
% Nathanael Kazmierczak, 08/16/2019
% UC Berkeley, Bediako Lab

% In the event that no sigma is passed, assume that no gaussian
% interpolation on the convolution map is desired.
%
% Really false_subpixel_mask is all that we will be using in the event that
% we are convolving something other than Bragg disks, like interference
% rods.
if nargin < 3
    sigma = 0;
end

TEST_FIGURES = 0;
DP_dim = size(DP);

% This line was really rather dumb.
% corr = xcorr2(double(uint8(DP)),double(probe_kernel));
conv_map = conv2(double(DP),double(probe_kernel),'same');  % restr
if TEST_FIGURES
    figure; image(conv_map,'CDataMapping','scaled'); shading flat; title 'Cross correlation matrix';
end
% Need to map back onto the original image. This line deserves further

if TEST_FIGURES
    figure; image(uint8(DP),'CDataMapping','scaled'); title 'Overlay';
    hold all; contour(conv_map/100);
    colormap jet;
end
%ax = gca;


%% Find maxima
if sigma ~= 0
        smootheddata = imgaussfilt(conv_map,sigma);
    else
        smootheddata = conv_map;
    end
if strcmp(subpixel,'False')
    false_subpixel_mask = imregionalmax(smootheddata);
    intensities_subpixel = [];    % Hacky: output arguments that are only used if we are using subpixel.
    coords_subpixel = [];
    intensities_pixel = [];
    coords_pixel = [];
else
    % For some reason, it seems like the import will only succeed if we are
    % in the correct folder, even if we have added it to the matlab path.
    current_dir = pwd;
    % directory change is contingent upon which system we are working on.
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/ClusteringStudies');
    [ platform_id ] = getPlatform();
    if platform_id == 1
        cd '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/py4DSTEM/py4DSTEM/process/utils'
    else
        cd 'C:\Users\Bediako Lab\Desktop\PythonWorkspace\py4DSTEM\py4DSTEM\process\utils'
    end
    py.importlib.import_module('utils');
    cd(current_dir);
    minans_subpixel = py.utils.get_maxima_2D(smootheddata,pyargs('subpixel',logical(1)));
    minans_pixel = py.utils.get_maxima_2D(smootheddata,pyargs('subpixel',logical(0)));
    % basically we will just use the non-subpixel algorithm with the mask
    % to determine if the peak is in the mask, and if so, we will use it.
    % Need to perform the following sequence of type conversions:
    % Python tuple > numpy n-D > python list > Matlab cell > Matlab double
    % array.
    y_coords_subpixel = cell2mat(cell(minans_subpixel{1}.tolist()));
	x_coords_subpixel = cell2mat(cell(minans_subpixel{2}.tolist()));
    intensities_subpixel = cell2mat(cell(minans_subpixel{3}.tolist()))';
    coords_subpixel = [y_coords_subpixel'+1, x_coords_subpixel'+1]; % keep these in indexing format, in which we will need them
    y_coords_pixel = cell2mat(cell(minans_pixel{1}.tolist()));
	x_coords_pixel = cell2mat(cell(minans_pixel{2}.tolist()));
    intensities_pixel = cell2mat(cell(minans_pixel{3}.tolist()))';
    coords_pixel = [y_coords_pixel'+1, x_coords_pixel'+1]; % keep these in indexing format, in which we will need them
end


end

