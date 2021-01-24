function [mask] = generateHexagonalLatticeMask(xspace,yspace,x0s,y0s,d)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% mask = isInHexagon(space(:,1),space(:,2),0,0,1);
assert(numel(x0s) == numel(y0s));
masks = zeros(size(xspace,1),size(xspace,2),numel(x0s));
for i = 1:numel(x0s)
    masks(:,:,i) = isInHexagon(xspace,yspace,x0s(i),y0s(i),d);
end
mask = sum(masks,3);

end

