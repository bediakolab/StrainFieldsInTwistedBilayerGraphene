% filter_validation.m
%
% 05/20/2020 Nathanael Kazmierczak

% [ interferometry_dataset, G1_atomic_coordinates, G2_atomic_coordinates ] = ...
%     generateTrainingData01( moire_angle, recon_angle, recon_distance, scan_dimensions,...
%     scan_step, scan_angle, beam_FWHM, max_intensity, add_noise_flag, make_plots_flag );

moire_angle = 1.2;
recon_angle = 0;
recon_distance = 45;
scan_dimensions = [120,120];
scan_step = 0.5;
scan_angle = 55;
beam_FWHM = 1.3;
max_intensity = 5*30;
add_noise_flag = true;
make_plots_flag = true;
raster_start_pos = [3,7];

use_old_reconfun = -1;  % NPK note: uses gaussRotateLattice03 -- should be robust to larger recon angles.
tblg = TwistedBilayerGrapheneAugmented(moire_angle,recon_angle,recon_distance,use_old_reconfun);
alg_flag = 0; % NPK note: The old and slow way, but more robust.
tblg.computeDSCField(alg_flag); 
tblg.computeNoModDSCField();
% G1_atomic_coordinates = tblg.getLayer1();
% G2_atomic_coordinates = tblg.getLayer2();
scan_centers = tblg.getRasterScan(raster_start_pos,scan_step,scan_angle,scan_dimensions);
scan_centers_mod = zeros(size(scan_centers));

% Do a scattered interpolant on the nomod DSC lattice
si_xcomp = scatteredInterpolant(tblg.layer1_fornomod(:,1),tblg.layer1_fornomod(:,2),tblg.noModDSCLat(:,1),'nearest');
si_ycomp = scatteredInterpolant(tblg.layer1_fornomod(:,1),tblg.layer1_fornomod(:,2),tblg.noModDSCLat(:,2),'nearest');
scan_centers_mod(:,1) = mod(scan_centers(:,1),tblg.cellDimXY(1));
scan_centers_mod(:,2) = mod(scan_centers(:,2),tblg.cellDimXY(2));
xnomoddisps = si_xcomp(scan_centers_mod(:,1),scan_centers_mod(:,2));
ynomoddisps = si_ycomp(scan_centers_mod(:,1),scan_centers_mod(:,2));

% [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( [xnomoddisps(:),ynomoddisps(:)] );
% reduced_xdisp = reshape(reduced_zone_disps(:,1),scan_dimensions);
% reduced_ydisp = reshape(reduced_zone_disps(:,2),scan_dimensions);
% 
y = FourDSTEM_Analysis_Engine;
% y.makeOutsideCustomDisplacementColorPlot(reduced_xdisp,reduced_ydisp);


%% Perform filters to simulate their impact
% First, blank filter to visualize dataset in a manner comparable to the
% other experiments.
filterstruct = buildFilterStructPreset( 'blank filter' );
plotflag = 3;
xnomoddisps = reshape(xnomoddisps,scan_dimensions);
ynomoddisps = reshape(ynomoddisps,scan_dimensions);
xnomoddisps = xnomoddisps + rand(size(xnomoddisps))*0.05;
ynomoddisps = ynomoddisps + rand(size(xnomoddisps))*0.05;
[ xdisp_filt,ydisp_filt ] = filterDisplacement( xnomoddisps,ynomoddisps,filterstruct,plotflag,y );

