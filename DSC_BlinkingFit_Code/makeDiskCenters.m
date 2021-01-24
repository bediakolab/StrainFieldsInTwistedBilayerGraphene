function [ disk_centers, radius, skipped_disks ] = makeDiskCenters( meanDP, exponent, hBN_mask )
% Nathanael Kazmierczak, Dec 2019
%
% Assume diffraction_patterns is a 4D array such as would come out of the
% original.
%
% The definition of the disks is such where the topmost disk in the inner
% ring should be done first and proceed clockwise, while the outer ring
% should go the same way except the first outer disk should be the one
% trailing the first inner one on the clockwise rotation.
if nargin < 3
    
    use_mask = false;
else 
    if isempty(hBN_mask)
        use_mask = false;
    else
        use_mask = true;
    end
end

DPtoplot = meanDP; %mean(mean(diffraction_patterns(:,:,:,:),3),4);
if use_mask
    DPtoplot(hBN_mask) = 0;
end
fighandle = plotDP(DPtoplot,exponent);
% xlim([2.263892128279883e+02 3.706107871720116e+02]);
% ylim([2.272580174927113e+02 3.680043731778425e+02]);
totalsize = size(meanDP);
disk_centers = zeros(6,2);
skipped_disks = [];

for i = 1:12
    fprintf('\n----- %dth disk -----.\n',i);
    input(sprintf('Zoom graph to prepare to click the %dth disks, and press any key when ready.',i));
    
    % NPK removed this section on 03/15/2020 because the averaging of the
    % diffraction patterns makes it obsolete.
%     while true
%         thisans = input('If the desired disk is not visible, enter 0; otherwise any entrance will continue.');
%         if thisans == 0
%             close(fighandle);
%             disp('Displaying a randomly different diffraction pattern.');
%             randind1 = randi(totalsize(3));
%             randind2 = randi(totalsize(4));
%             DPtoplot = meanDP; %diffraction_patterns(:,:,randind1,randind2);
%             if use_mask
%                 DPtoplot(hBN_mask) = 0;
%             end
%             fighandle = plotDP(DPtoplot,exponent);
%         else
%             break
%         end
%     end
    
    if i == 1
        disp('Click once for the first disk (top vertical position) and a second time to set the circular radius.');
        setcirc = ginput(2);
        radius = sum((setcirc(1,:) - setcirc(2,:)).^2).^0.5;
        disk_centers(1,:) = setcirc(1,:);
    else
        
        tf = input(sprintf('1/0: Is the %dth disk to be skipped? ',i));
        if tf
            disk_centers(i,:) = nan;
            skipped_disks = [skipped_disks,i];
        else
            fprintf('Click once to set the center of the %dth disk.\n',i);
            disk_centers(i,:) = ginput(1);
        end
        
    end
end
viscircles(disk_centers,radius*ones(size(disk_centers,1),1));



end

