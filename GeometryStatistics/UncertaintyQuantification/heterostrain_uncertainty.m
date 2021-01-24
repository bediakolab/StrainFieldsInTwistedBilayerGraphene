% heterostrain_uncertainty.m
%
% Probing the effect of AA domain registration uncertainty on the
% distribution of heterostrain % values calculated from the triangulation.
% Nathanael Kazmierczak, 10/23/2020

% thetam_deg = 0.63;
thetam_deg = 1.37;

thetam_rad = deg2rad(thetam_deg);
[ ~, ~, hexagon_lattice_constant ] = getDSCBasisVectors();
a = hexagon_lattice_constant/10;  % in nm
moire_wavelength = a/(2*sin(thetam_rad/2));

b1 = [moire_wavelength;0];
ra = pi/3;
rotmat = [cos(ra), -sin(ra); sin(ra), cos(ra)];
b2 = rotmat*b1;
% rectbasis_b1 = -50:100;
% rectbasis_b2 = -50:100;
rectbasis_b1 = -25:50;
rectbasis_b2 = -25:50;
[allpoints_b1,allpoints_b2] = meshgrid(rectbasis_b1,rectbasis_b2);
rectbasis = [allpoints_b1(:),allpoints_b2(:)];

hexlattice = ([b1,b2]*rectbasis')';
% crop
maxval = 500;
hexlattice(hexlattice(:,1) > maxval | hexlattice(:,1) < 0 | hexlattice(:,2) > maxval | hexlattice(:,2) < 0,:) = [];

figure
scatter(hexlattice(:,1),hexlattice(:,2));
axis equal;
title('Zero-noise AA lattice');

% Perform a random perturbation of the centers:
% sigma = 0.3; % nm
sigma = 0.1;
hexlattice_noisy = hexlattice + normrnd(0,sigma,size(hexlattice));

figure
scatter(hexlattice_noisy(:,1),hexlattice_noisy(:,2));
axis equal;
title(sprintf('Noisy AA lattice, sigma = %.2f',sigma));

%% Triangulate and get Moire wavelengths
% Strategy here is to load into instance variables.
m4 = FourDSTEM_Analysis_Engine;
m4.datacube_size = [100,100,100,100];
m4.scan_stepsize = 1;
DT = delaunayTriangulation(hexlattice_noisy);
m4.AAtriangulation = DT;
m4.setFacetTwistAngles();
m4.repairTriangulation();
% [average_angle,std_angle,max_angle,min_angle,n_triangles_used,angles_used] = m4.plotFacetedTwistAngles([0,inf],parula,false,false,false);
angle_guess = thetam_deg;
average_angle_range = [];
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);
make_plots = 1;
[av_moire_angle,std_moire_angle,min_moire_angle,max_moire_angle,...
                av_strain_percent,std_strain_percent,min_strain_percent,max_strain_percent,copypaste] = ...
                m4.getUniaxialStrainStatistics(make_plots)%,moire_binedges,strainangle_binedges,strainval_binedges);

removeLines = true;     
cmap = magma;
m4.plotFacetedUniaxialStrain(cmap,removeLines);

% average_angle
% std_angle
% max_angle
% min_angle
% n_triangles_used
