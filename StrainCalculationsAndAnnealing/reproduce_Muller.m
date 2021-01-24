% reproduce_Muller.m
%
% Script for reproducing the strain calculations from Muller's 2013 PNAS
% paper.
%
% Nathanael Kazmierczak

%% Determined from STEM experiments.
FWHM_tensile = 10.1;  % +/- 1.4nm
FWHM_shear = 6.2;  % +/- 0.6nm

%% Use calibration from multislice simulations to get soliton widths
A1 = 1.458;
A0 = 0.099;
soliton_width_rln = @(FWHM) A1*FWHM + A0;
soliton_shear_width = soliton_width_rln(FWHM_shear);
soliton_tensile_width = soliton_width_rln(FWHM_tensile);

%% Use sine-gordon equation to get displacement (u(x)):
syms x w
delu = 1 + 2/pi*atan(exp(pi*x/w));
strain = diff(delu,x);
% Plot strain for a given ws, say 2;
strain_ws2 = subs(strain,w,2);
strainfun = matlabFunction(strain_ws2);
xbase = -5:0.01:5;
strainvals = strainfun(xbase);
figure
plot(xbase,strainvals,'-o');
% Go for a symbolic solution to the maximum strain, assuming it's always at
% x = 0 with this equation
strain_center = subs(strain,x,0);
strain_w = matlabFunction(strain_center);

%% strain function
units_strain_percent = @(w_soliton) 0.07104/w_soliton * 100
tensile_max_strain = units_strain_percent(soliton_tensile_width)
shear_max_strain = units_strain_percent(soliton_shear_width)

