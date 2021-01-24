function [strain_transform] = buildStrainTransform(e_xx,e_yy,e_xy,theta)
% Function for building the linear transformation matrix for performing
% linear transformation of the bragg vectors under infinitesimal strain
% theory.
%
% Theta should presumably be in radians. This is strain in reciprocal
% space.
%
% See Pekin et al., Optimizing Disk Registration Algorithms for Nanobeam
% Electron Diffraction Strain Mapping, Ultramicroscopy 2017.

strain_transform = zeros(2,2);
strain_transform(1,1) = e_xx;
strain_transform(1,2) = 0.5*(e_xy - theta);
strain_transform(2,1) = 0.5*(e_xy + theta);
strain_transform(2,2) = e_yy;
strain_transform = strain_transform + eye(2);

end

