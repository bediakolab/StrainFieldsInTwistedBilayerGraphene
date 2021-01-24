function [ com_coords ] = getCOMofMASK( mask )
% Refactored out of getROIKernel.m

linlogical = mask(:);
lininds = (1:numel(linlogical))';
[Iinds,Jinds] = ind2sub(size(mask),lininds(linlogical));
com_coords = mean([Iinds,Jinds],1);  % center of mass.


end

