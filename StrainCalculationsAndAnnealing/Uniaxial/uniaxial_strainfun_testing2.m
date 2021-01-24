% uniaxial_strainfun_testing2.m
%
% Function to ascertain how the uniaxial function responds to different
% angles.
%
% Nathanael Kazmierczak, 06/02/2020

angles = 0:360;
wavelength_stor = zeros(numel(angles),3);
theta_moire = 1;
epsilon_strain_percent = 0.3;
for i = 1:numel(angles)
    this_theta_strain = angles(i); 
    [ moire_wavelengths ] = uniaxialStrainPredFun2( theta_moire, this_theta_strain, epsilon_strain_percent );
    wavelength_stor(i,:) = moire_wavelengths;
end

figure;
plot(angles',wavelength_stor,'-o');
legend('Side 1', 'Side 2','Side 3');
