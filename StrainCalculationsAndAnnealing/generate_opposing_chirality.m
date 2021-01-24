% generate_opposing_chirality.m
%
% This script attempts to build a dataset that has saddle point chirality
% matching that of our experimental datasets. The goal is to shed light on
% why the varying sign conventions arise.
%
% Nathanael Kazmierczak, 05/29/2020, Bediako Lab, UC Berkeley.

%% Dataset construction
moire_angle = -1.1;
recon_angle = -1.0;
recon_distance = 45;
scan_dimensions = [200,200];
scan_step = 0.5;
scan_angle = 55;
add_noise_amount = 0;
plot_sample = true;
raster_start_pos = [5,7];


[ displacement_field ] = makeSimulatedDisplacementDataset(moire_angle, recon_angle, recon_distance, scan_dimensions,...
    scan_step, scan_angle, add_noise_amount, plot_sample, raster_start_pos);

