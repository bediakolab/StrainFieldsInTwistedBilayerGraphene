function [ roi_data ] = roiBoilerplate(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

overlay = sum(obj.disk_averages,3);
figure;
imagesc(overlay);
title('Overlay of all disks');
colormap(gray);
BW = false(obj.datacube_size(1),obj.datacube_size(2));
while true
    disp('Please select a region of known AB stacking.');
    thisBW = roipoly();
    BW = BW & thisBW;
    yn = input('Do you wish to select another AB region? 1/0');
    if ~yn
        break
    end
end
roi_data = zeros(nnz(BW),12,'uint16');
for i = 1:12
    disk_data = obj.disk_averages(:,:,i);
    roi_data(:,i) = disk_data(BW);
end

end

