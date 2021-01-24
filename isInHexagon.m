function [tf] = isInHexagon(x,y,x0,y0,d)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
latticepoints = getHexagon(x0,y0,d);
[xsort] = sortrows(latticepoints);
finalsort = [xsort(1:2,:); xsort(4,:); xsort(6,:); xsort(5,:); xsort(3,:); xsort(1,:)];

cond1 = (x < x0 & x > (x0 - d)) & ((y < 1/sqrt(3)*(x-xsort(2,1)) + xsort(2,2)) & (y > -1/sqrt(3)*(x-xsort(1,1)) + xsort(1,2)));
cond2 = (x >= x0 & x <= (x0 + d)) & ((y > 1/sqrt(3)*(x-xsort(3,1)) + xsort(3,2)) & (y < -1/sqrt(3)*(x-xsort(4,1)) + xsort(4,2)));
% 
% if x < x0 & x > (x0 - d)
%     tf = ;
% elseif x > x0 & x <= (x0 + d)
%     tf = ((y > 1/sqrt(3)*(x-xsort(3,1)) + xsort(3,2)) & (y < -1/sqrt(3)*(x-xsort(4,1)) + xsort(4,2)));
% else
%     tf = logical(0);
% end
tf = cond1 | cond2;

end

