function [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors()
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

hexagon_lattice_constant = 2.461;
t = 60/180*pi;
rotmat = [cos(t) sin(t); -sin(t) cos(t)];
v1 = [0,hexagon_lattice_constant/sqrt(3)];
v2 = v1*rotmat';

end

