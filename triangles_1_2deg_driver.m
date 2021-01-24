% triangles_1_2deg_driver.m
%
% Script which attempts an unconventional analysis of the zoomed in data
% where we only have two Bragg peaks to work with per lattice.
% 
% Contains code moved from beamstop_zoomedin_center.m, because that one was
% not named quite so well.
%
% Nathanael Kazmierczak, 08/16/2019
% Bediako Lab, UC Berkeley

addpath('C:\Users\Bediako Lab\Desktop\4DSTEM data');
addpath('C:\Users\Bediako Lab\Desktop\4DSTEM data\Maddiedata07292019');
addpath('/Volumes/Lexar/BediakoLabLexar');

%% Load
start = [1,1,1,1];
% start = [1,1,30,50];  % Trying to test one of the ones that was failing.
count = [1024,1024,5,5];  % the full dimensions of this cropped piece of data.
% data = h5read('6__25x25_ss4nm_2s_spot4_alpha=0p14_bin=4_CL=380_80kV_C2=10um.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
% data = ReadDMFile('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um.dm4');
% filename = '7_5x5_ss10nm_8s_spot4_alpha=0p14_bin=2_CL=380_80kV_C2=10um.h5';  % While only 5x5, this has the binning of two rather than the binning by four.
filename = '7_5x5_ss10nm_8s_spot4_alpha=0p14_bin=2_CL=380_80kV_C2=10um.h5';  % While only 5x5, this has the binning of two rather than the binning by four.
data = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
firstDP = data(:,:,1,1); %data(:,:,1);

%% Upsample the image, because I think we are going to need more resolution.
upsample_factor = 4;
firstDP = imresize(firstDP,upsample_factor);
imsize = size(firstDP);
%% Visualize
figure;
image(uint8(firstDP),'CDataMapping','scaled');
colormap(pink);
title('Uint8 visualization')

figure;
image(firstDP,'CDataMapping','scaled');
colormap(pink);
title('Uint16 visualization')

figure;
power = 1;
image(uint16(double(firstDP).^power),'CDataMapping','scaled');
colormap(pink);
title('Uint16 scaled visualization')

%% Build out the kernel from an ROI
ROIkernel = 1;
flat_kernel = 1;  % Rather than taking the actual pixels from the ROI for the kernel, use the mask.
boundary_kernel = 1;
boundary_thickness = 1; % number of boundarymask calls which are employed
if ROIkernel
    %%%%%% NOTE 09032019 this code was refactored into getROIKernel.m
    %%%%%% function, so that it could be employed in the attempts to
    %%%%%% quantitatively register the interference patterns.
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
    if flat_kernel
        masked_kernel = polymask;
    else
        masked_kernel = firstDP;
        masked_kernel(~polymask) = 0;
    end
    if boundary_kernel
        bmasks = zeros(size(firstDP,1),size(firstDP,2),boundary_thickness);
        for i = 1:boundary_thickness
            if i == 1
                bmasks(:,:,1) = boundarymask(masked_kernel);
            else
                bmasks(:,:,i) = boundarymask(bmasks(:,:,i-1));
            end
        end
        % When done, add them all up to make a new mask.
        logicalmask = sum(bmasks,3);
        logicalmask(logicalmask > 1) = 1;
        masked_kernel(~logicalmask) = 0;  % So under this implementation, the kernel will never grow; just cutting a hole out of the middle.
    end
    % This seemed too confusing.
% %     elseif ~flatkernel
% %         % If we are doing boundary masks, we want to have consistent
% %         % behavior between flat and non-flat. In both cases you are
% %         % basically drawing the middle of the mask.
% %         
% %     end
        
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

%% Build masks for triangle registration.
numberofsites = 2;
kernel_radius = kernel;  % the hacky way of getting the ROI kernel into the function
sigma = 0;
min_disk_distance = 20;
num_to_find_per_window = 2;
make_figures = 1;
min_intens_rel_thres = 0.01;
max_intens_rel_thres = 100;
subpixel = 'False';
[hexagonal_window_masks] = getFreehandMasks(firstDP,numberofsites);
[disk_locations, disk_intensities, com_coordinates] = braggDiskRegistrationRecipie_1(firstDP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,min_intens_rel_thres,max_intens_rel_thres,subpixel);


%% Build masks on the center.
use_teeny_mask = 0;
teeny_mask_radius_factor = 0;
draw_hard_delete_region = 1;
% First need to apodize so that the circle can be drawn, it seems.
ypadding = imsize(1)*2;
padzeros = zeros(ypadding,imsize(2));
paddedDP = vertcat(firstDP,padzeros);

THRESHOLD = 150;
[beamstop_masks,Irng_beamstopmask,Jrng_beamstopmask,this_maskeddata] = makeBeamstopCenterMasks(paddedDP,...
    use_teeny_mask,teeny_mask_radius_factor,draw_hard_delete_region,THRESHOLD);



