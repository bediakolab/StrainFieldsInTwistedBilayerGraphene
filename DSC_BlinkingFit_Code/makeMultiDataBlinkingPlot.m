function figh = makeMultiDataBlinkingPlot( disk_averages,shadingtype,normalize,xaxis,yaxis )
%UNTITLED Summary of this function goes here
%   disk_averages should be a 3D array with size3 = 6, for the six unique
%   disks and the proper subplotting dimensions.
%
%   assume these are disks 1:3, 7:9
%
%   Still need to center these because of elastic scattering

figh = figure;
k = size(disk_averages,3);
assert(k == 6);
diskstor = [1:3,7:9];

for i = 1:6
    this_disk_average = disk_averages(:,:,i);
    % NPK took the following line out 02/02/2020 because, while it makes
    % for prettier plots, it was not how the fitting code worked.
    % Reinstated 02/09/2020
    if normalize
        this_disk_average = this_disk_average - min(min(this_disk_average));
        this_disk_average = this_disk_average./repmat(max(max(this_disk_average)),size(this_disk_average,1),size(this_disk_average,2));
    end
%     this_disk_average = this_disk_average
    subplot(2,3,i)
    disknum = diskstor(i);
%     disk_averages = disk_average./max(disk_average);
    
    if strcmp(shadingtype,'interp')
        contourf(xaxis,yaxis,this_disk_average,50,'LineStyle','None');
        shading(gca,shadingtype);
    elseif strcmp(shadingtype,'flat')
        pcolor(xaxis,yaxis,this_disk_average);
        axis tight
        shading flat;
    end
    colormap(gray);
    colorbar
    axis equal
    title(sprintf('Disk %d',disknum))
end


end

