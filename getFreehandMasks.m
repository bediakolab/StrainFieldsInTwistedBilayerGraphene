function [freehand_window_masks] = getFreehandMasks(plotdata,number_of_sites)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Do this as six for now for consistency, but may want only five later on.
freehand_window_masks = zeros(size(plotdata,1),size(plotdata,2),number_of_sites);
figure;
image(uint8(plotdata),'CDataMapping','scaled');
colormap pink

disp('This function is for developing windows around clusters of Bragg disks for detection.');
fprintf('You will be prompted to draw %d polygons for the masks.\n',number_of_sites);

for i = 1:number_of_sites
    fprintf('Please draw polyhon %d\n',i);
    mypoly = drawpolygon;
    freehand_window_masks(:,:,i) = createMask(mypoly);
end


end

