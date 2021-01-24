% beamstop_zoomed_in_center.m
%
% Script for ascertaining how well and consistently we can power law fit
% the center of a diffraction pattern where the camera was zoomed in on the
% side.
%
% Nathanael Kazmierczak, 08/12/2019
% Bediako Lab, UC Berkeley

%% Load dataset 7 for these tests.
addpath Maddiedata07292019
filename = '7_5x5_ss10nm_8s_spot4_alpha=0p14_bin=2_CL=380_80kV_C2=10um.h5';
info = h5info(filename);
h5disp(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data');

% Read in one of these larger, 1024x1024 diffraction patterns. Start here
% to avoid the weird triangle artifact that I recall having seen on some of
% the initial DPs.
start = [1,1,2,2];
count = [1024,1024,1,1];
DP = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);

    imsize = size(DP);

%% Visualize
figure;
image(uint8(DP),'CDataMapping','scaled');
colormap(pink);
title('Uint8 visualization')

figure;
image(DP,'CDataMapping','scaled');
colormap(pink);
title('Uint16 visualization')

figure;
image(uint16(double(DP).^(1/5)),'CDataMapping','scaled');
colormap(pink);
title('Uint16 scaled visualization')

%% Build out the kernel from an ROI
ROIkernel = 0;
if ROIkernel
    disp('Please draw a polygon surrounding your kernel.');
    mypoly = drawpolygon;
    disp('Please select the center of the kernel.');
    [center_coords] = ginput(1);
    % Think about this further if it seems useful -- there could be ways of
    % defining the kernel center that don't require the pixel rounding.
    global_kernel_center_coords = round(center_coords);
    % Need to convert the specified kernel into an image of the same original
    % size (an apodization if nothing else).
    polymask = createMask(mypoly);
    masked_kernel = DP;
    masked_kernel(~polymask) = 0;
    [Icoords,Jcoords] = ind2sub(size(polymask),find(polymask));
    Irng = [min(Icoords),max(Icoords)];
    Jrng = [min(Jcoords),max(Jcoords)];
    trunc_polymasked = masked_kernel(Irng(1):Irng(2),Jrng(1):Jrng(2));
    % now we need to add zeros in each dimension to fill up the full measure of
    % the image size.
    % I_left_needed = imsize(1)/2 - 1; % because the "one" is occupied by the kernel center.
    % I_right_needed = imsize(1)/2;
    % J_top_needed = imsize(2)/2;
    % J_bottom_needed = imsize(2)/2 - 1;
    % I_left_current = global_kernel_center_coords(1) - Irng(1);
    % I_right_current = Irng(2) - global_kernel_center_coords(1);
    % J_top_current = global_kernel_center_coords(2) - Jrng(1);
    % I_bottom_current = Jrng(2) - global_kernel_center_coords(2);
    I_offset = global_kernel_center_coords(2) - Irng(1);  % because ginput comes out as xy, not ij
    J_offset = global_kernel_center_coords(1) - Jrng(1);
    I_length = Irng(2) - Irng(1);
    J_length = Jrng(2) - Jrng(1);
    
    kernel = zeros(imsize(1),imsize(2));
    new_center_position = [imsize(1)/2,imsize(2)/2];
    new_Ispan = (new_center_position(1)-I_offset):((new_center_position(1)-I_offset)+I_length);
    new_Jspan = (new_center_position(2)-J_offset):((new_center_position(2)-J_offset)+J_length);
    
    kernel(new_Ispan,new_Jspan) = trunc_polymasked;
    figure;
    image(kernel,'CDataMapping','scaled');
    % This seems like it is built correctly.
    % Though come to think of it, it isn't actually needed for fitting the
    % center, but we will want it when we try to get a strain map.
end

%% Build masks on the center.
use_teeny_mask = 0;
teeny_mask_radius_factor = 0;
draw_hard_delete_region = 1;
% First need to apodize so that the circle can be drawn, it seems.
ypadding = 1024;
padzeros = zeros(ypadding,imsize(2));
paddedDP = vertcat(DP,padzeros);

[beamstop_masks,Irng_beamstopmask,Jrng_beamstopmask,this_maskeddata] = makeBeamstopCenterMasks(paddedDP,...
    use_teeny_mask,teeny_mask_radius_factor,draw_hard_delete_region);


