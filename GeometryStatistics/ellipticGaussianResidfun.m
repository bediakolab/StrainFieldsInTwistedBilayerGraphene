function [ residvals, J ] = ellipticGaussianResidfun( c, space, amplitudes, mask )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[ predvals,J ] = ellipticGaussianPredfun( space, c(1),c(2),c(3),c(4),c(5), c(6), c(7) );
residvals = mask(:).*(amplitudes(:) - predvals);
J = -J.*repmat(mask(:),1,7);
end

