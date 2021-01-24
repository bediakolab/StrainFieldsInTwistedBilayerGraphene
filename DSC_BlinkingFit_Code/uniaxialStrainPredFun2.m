function [ moire_wavelengths ] = uniaxialStrainPredFun2( theta_moire, theta_strain, epsilon_strain_percent )
% Fitting function for computing uniaxial heterostrain of a twisted bilayer
% graphene sample.
%
% This function is identical to uniaxialStrainPredFun() except it doesn't
% have the sort before the end. So there is not a 60 degree angle
% redundancy, but rather a 180 degree angle redundancy.
%
% Nathanael Kazmierczak, 06/02/2020

theta_moire = deg2rad(theta_moire);
theta_strain = deg2rad(theta_strain);
epsilon_strain = epsilon_strain_percent/100;

[ ~, ~, a0 ] = getDSCBasisVectors();
k = 4*pi/(sqrt(3)*(a0/10));  % reciprocal lattice vector, in nm-1 because nm are the units of the Moire wavelength measurement.
del = 0.16;  % The Poisson ratio for graphene
R = @(t) [cos(t), -sin(t);sin(t) cos(t)];
S = @(theta_n,eps) R(-theta_n)*diag([1/(1+eps),1/(1-del*eps)])*R(theta_n);
k1_template = [k;0];
k2_template = R(pi/3)*k1_template;
k3_template = R(2*pi/3)*k1_template;
% Consider uniaxial heterostrain. k1s, k2s, k3s constitute a
% strained lattice with strain applied at an angle theta_s
% relative to the x-axis. k1p, k2p, k3p constitute an
% unstrained lattice rotated by the Moire angle theta_t. We are
% calculating the magnitude of the resulting Moire wavelengths,
% so the direction of the lattice doesn't matter here -- just
% be sure to sort the values before fitting.
k1s = S(theta_strain,epsilon_strain)*k1_template;
k2s = S(theta_strain,epsilon_strain)*k2_template;
k3s = S(theta_strain,epsilon_strain)*k3_template;
k1p = R(theta_moire)*k1_template;
k2p = R(theta_moire)*k2_template;
k3p = R(theta_moire)*k3_template;
K1s = k1p - k1s;
K2s = k2p - k2s;
K3s = k3p - k3s;
M1mag = 4*pi/(sqrt(3)*norm(K1s));
M2mag = 4*pi/(sqrt(3)*norm(K2s));
M3mag = 4*pi/(sqrt(3)*norm(K3s));
% moire_wavelengths = sort([M1mag,M2mag,M3mag]);
moire_wavelengths = [M1mag,M2mag,M3mag];

end

