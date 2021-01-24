function [ thickened_mask ] = thickenMask( mask, thickening_iterations )
% thickenMask.m

% need to convert shape to true binary mask if it isn't already.
mask = logical(mask); % this seems to work quickly.
bmstor = zeros(size(mask,1),size(mask,2),thickening_iterations);
for i = 1:thickening_iterations
    lastmat = max(bmstor,[],3);  % A safeguard to make sure nothing falls out. zeros are fine
    bmstor(:,:,i) = boundarymask(lastmat | mask);
end
thickened_mask = max(bmstor,[],3) | mask;

end

