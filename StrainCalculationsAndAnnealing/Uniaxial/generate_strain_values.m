% generate_strain_values.m
%
% Script to obtain precise sidelengths for refined FEM simulations.
%
% Nathanael Kazmierczak, 07/28/2020.

ticks = [25 26 19];
nms = sort(ticks*20/32)

% tm = 0.95;
tm = 1;
ts = 17;
eta = 0.7;

[ moire_wavelengths ] = sort(uniaxialStrainPredFun( tm, ts, eta ))
% So under 0.95 deg the first set of wavelength conditions to use would be [11.8,15.8,
% 16.6] nm.

% Under 1 deg moire, the set would be 11.3, 15.0, 15.7. nm

%% Second set of conditions: 
ticks2 = [40, 41, 29];
nms2 = sort(ticks2*20/32)

tm2 = 0.63;
ts2 = 20;
eta2 = 0.45;
[ moire_wavelengths2 ] = sort(uniaxialStrainPredFun( tm2, ts2, eta2 ))

% Use 17.9 nm, 24.3 nm, 24.4 nm


%% Revised 1.1 degree, unstrained:
tm = 1.1;
ts = 17;
eta = 0.7;
[ moire_wavelengths_str ] = sort(uniaxialStrainPredFun( tm, ts, eta ))

eta = 0;
[ moire_wavelengths_unstr ] = sort(uniaxialStrainPredFun( tm, ts, eta ))

%% Extra stuff Bediako wants
tm = 1.5;
ts = 17;
eta = 0.7;
[ moire_wavelengths_str ] = sort(uniaxialStrainPredFun( tm, ts, eta ))

eta = 0;
[ moire_wavelengths_unstr ] = sort(uniaxialStrainPredFun( tm, ts, eta ))

tm = 2.0;
ts = 17;
eta = 0.7;
[ moire_wavelengths_str ] = sort(uniaxialStrainPredFun( tm, ts, eta ))

eta = 0;
[ moire_wavelengths_unstr ] = sort(uniaxialStrainPredFun( tm, ts, eta ))

