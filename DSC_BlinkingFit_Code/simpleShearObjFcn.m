function [ normal_strain_residuals ] = simpleShearObjFcn( theta, tensormat )
% Nathanael Kazmierczak, 06/15/2020
% Prepares an optimization to rotate coordinate axes such that the amount
% of normal strain is minimized.

[ rotated_tensormat ] = simpleShearPredfun( theta, tensormat );
exx = rotated_tensormat(1,1);
eyy = rotated_tensormat(2,2);
% normal_mag = sqrt(exx^2 + eyy^2);
normal_strain_residuals = [exx,eyy];  % for lsqnonlin syntax

end

