% extract_Stephen_strain_values.m
%
% Getting heterostrain values for Stephen's simulations.
%
% Nathanael Kazmierczak, 07/28/2020

%% Stephen heterostrained #1
[v1,v2,a] = getDSCBasisVectors();
thetam_1 = 1.03;
sl = a/deg2rad(thetam_1)/10;

hetstrain1 = [17.8947,20.5263,13.6842];
m4 = FourDSTEM_Analysis_Engine();
m4.triangle_sidelengths = hetstrain1;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain([0,inf],1);
m4.uniaxial_moire_angles
% obj.uniaxial_strain_angles = uniaxial_strain_angle_storage;
m4.uniaxial_strain_percents
% obj.uniaxial_strain_residuals = bestresidstor;

%% My atomic heterstrain simulator (ThetaM = 2.0)
tm = 2;
slpred = a/deg2rad(tm)/10;
ticks = [24 18 16];
nms = ticks*15/39;

m4 = FourDSTEM_Analysis_Engine();
m4.triangle_sidelengths = nms;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain([0,inf],2);
m4.uniaxial_moire_angles
m4.uniaxial_strain_percents

%% Stephen #2 (Labeled as 0.63)
thetam_2 = 0.63;
slpred2 = a/deg2rad(thetam_2)/10;

ticks = [22.5 26 17];
nms = ticks*20/15;
m4 = FourDSTEM_Analysis_Engine();
m4.triangle_sidelengths = nms;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain([0,inf],1);
m4.uniaxial_moire_angles
% obj.uniaxial_strain_angles = uniaxial_strain_angle_storage;
m4.uniaxial_strain_percents


%% Stephen #3 (labeled as 0.26)
thetam_3 = 0.26;
slpred3 = a/deg2rad(thetam_3)/10;

ticks = [16 28 32];
nms3 = ticks*50/15;
m4 = FourDSTEM_Analysis_Engine();
m4.triangle_sidelengths = nms3;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain([0,inf],0.2);
m4.uniaxial_moire_angles
% obj.uniaxial_strain_angles = uniaxial_strain_angle_storage;
m4.uniaxial_strain_percents


%% Replicate Kerelsky paper
% Case 1
s1 = [13.72,12.7,10.18];
m4 = FourDSTEM_Analysis_Engine();
m4.triangle_sidelengths = s1;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain([0,inf],0.2);
m4.uniaxial_moire_angles
% obj.uniaxial_strain_angles = uniaxial_strain_angle_storage;
m4.uniaxial_strain_percents
slav = a/deg2rad(m4.uniaxial_moire_angles)/10
% This checks out


% Case 2
s1 = [14.5,13.2,10.84];
m4 = FourDSTEM_Analysis_Engine();
m4.triangle_sidelengths = s1;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain([0,inf],0.2);
m4.uniaxial_moire_angles
% obj.uniaxial_strain_angles = uniaxial_strain_angle_storage;
m4.uniaxial_strain_percents
% This also checks out
slav = a/deg2rad(m4.uniaxial_moire_angles)/10

