function [ rotated_tensormat ] = simpleShearPredfun( theta, tensormat )
% Nathanael Kazmierczak, 06/15/2020

Q = [cos(theta),-sin(theta);sin(theta),cos(theta)];
rotated_tensormat = Q'*tensormat*Q;

end

