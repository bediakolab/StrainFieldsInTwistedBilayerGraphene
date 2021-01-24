function [ kernel ] = getROIKernel( DP, flat_kernel, boundary_kernel )
% Function developed 09032019, refactored out of triangles_1_2deg_driver.m
%
% Nathanael Kazmierczak, Bediako Lab, UC Berkeley

figure; image(uint8(DP),'CDataMapping','Scaled'); colormap(hsv); colorbar; pbaspect([1 1 1]);
disp('Please draw a polygon surrounding your kernel.');
polymask = roipoly;
[ center_coords ] = getCOMofMASK( polymask );
% % linlogical = polymask(:);
% % lininds = (1:numel(linlogical))';
% % [Iinds,Jinds] = ind2sub(size(polymask),lininds(linlogical))
% % center_coords = mean([Iinds,Jinds],1);  % center of mass.

% Think about this further if it seems useful -- there could be ways of
% defining the kernel center that don't require the pixel rounding.
global_kernel_center_coords = round(center_coords);
% Need to convert the specified kernel into an image of the same original
% size (an apodization if nothing else).
if flat_kernel
    masked_kernel = polymask;
else
    masked_kernel = DP;
    masked_kernel(~polymask) = 0;
end
if boundary_kernel
    bmasks = zeros(size(DP,1),size(DP,2),boundary_thickness);
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

imsize = size(DP);
kernel = zeros(imsize(1),imsize(2));
if ~mod(imsize(1),2) && ~mod(imsize(2),2)
    new_center_position = [imsize(1)/2,imsize(2)/2];
elseif mod(imsize(1),2) && ~mod(imsize(2),2)
    new_center_position = [(imsize(1)+1)/2,imsize(2)/2];
elseif ~mod(imsize(1),2) && mod(imsize(2),2)
    new_center_position = [imsize(1)/2,(imsize(2)+1)/2];
else
    new_center_position = [(imsize(1)+1)/2,(imsize(2)+1)/2];
end
new_Ispan = (new_center_position(1)-I_offset):((new_center_position(1)-I_offset)+I_length);
new_Jspan = (new_center_position(2)-J_offset):((new_center_position(2)-J_offset)+J_length);

kernel(new_Ispan,new_Jspan) = trunc_polymasked;
figure;
image(kernel,'CDataMapping','scaled'); colormap(hsv); colorbar; pbaspect([1 1 1]);

end

