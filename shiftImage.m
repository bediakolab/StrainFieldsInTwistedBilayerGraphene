function [ shiftedDPROI ] = shiftImage(this_DPROI,shiftx,shifty)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

shiftedDPROI = zeros(size(this_DPROI));
for i = 1:size(this_DPROI,1)
    for j = 1:size(this_DPROI,2)
        
        test1 = i+shifty > 0 && i+shifty <= size(shiftedDPROI,1);
        test2 = j+shiftx > 0 && j+shiftx <= size(shiftedDPROI,2);
        if test1 && test2
            shiftedDPROI(i+shifty,j+shiftx) = this_DPROI(i,j);
        end
        
% % %         try
% % %             shiftedDPROI(i+shifty,j+shiftx);  % Is there an element to be accessed
% % %             shiftedDPROI(i+shifty,j+shiftx) = this_DPROI(i,j);
% % %         catch  % in this case we have gone beyond the bounds of the shiftedDPROI
% % %             % Don't have to do anything! it was initialized to zero.
% % %         end
    end
end

% disp('Ending shifting function');

end

