function [ displacement_field ] = makeSimulatedDisplacementDataset(moire_angle, recon_struct, heterostrain_struct, scan_dimensions,...
    scan_step, scan_angle, add_noise_amount, plot_sample, raster_start_pos, customCellDimXY)
% This function does a scattered interpolation of the displacement field to
% generate a reduced-zone displacement field. There is no beamwidth
% averaging here.
%
% Compare to generateTrainingData01.m in the machine learning folder, which
% is designed to generate interferometry data and atomic coordinates, and a
% forthcoming function makeBeamwidthAveragedDisplacementDataset.m, which
% will then fit that interferometry to the cosine^2 model to produce a
% displacement dataset, thereby completing the full simulation process.
%
% Here noise will be added to the displacement vectors themselves according
% to a normal distribution. This does not capture biasing effects in the
% same way that makeBeamwidthAveragedDisplacementDataset.m will.
%
% Nathanael Kazmierczak, 05/30/2020, Bediako Lab, UC Berkeley

% use_old_reconfun = -1;  % NPK note: uses gaussRotateLattice03 -- should be robust to larger recon angles.
if nargin < 10
    customCellDimXY = [];
end
tblg = TwistedBilayerGrapheneAugmented(moire_angle,recon_struct,heterostrain_struct,customCellDimXY);
alg_flag = 0; % NPK note: The old and slow way, but more robust.
tblg.computeDSCField(alg_flag); 
[~,scan_centers_mod] = tblg.getRasterScan(raster_start_pos,scan_step,scan_angle,scan_dimensions);

% Do a scattered interpolant on the reduced zone displacement field
% Note: This is different than in filter_validation.m. There, the
% interpolation was on the extended zone DSC lattice. A problem there,
% though, is that the extended zone field computation is prone to
% artifacts, whereas the reduced zone one is less so.
layer1coords = tblg.getLayer1();
layer2coords = tblg.getLayer2();
DSC_lat = tblg.getDSClat();
% NPK caught a minor bug here 07/13/2020. May have been responsible for
% some of the asymmetries observed in simulated displacement data.
% si_xcomp = scatteredInterpolant(layer1coords(:,1),layer2coords(:,2),DSC_lat(:,1),'nearest');
% si_ycomp = scatteredInterpolant(layer1coords(:,1),layer2coords(:,2),DSC_lat(:,2),'nearest');
si_xcomp = scatteredInterpolant(layer1coords(:,1),layer1coords(:,2),DSC_lat(:,1),'nearest');
si_ycomp = scatteredInterpolant(layer1coords(:,1),layer1coords(:,2),DSC_lat(:,2),'nearest');
xdisps = si_xcomp(scan_centers_mod(:,1),scan_centers_mod(:,2));
ydisps = si_ycomp(scan_centers_mod(:,1),scan_centers_mod(:,2));

if add_noise_amount > 0
    [r,c] = size(xdisps);
    xdisps = xdisps + normrnd(0,add_noise_amount,r,c);
    ydisps = ydisps + normrnd(0,add_noise_amount,r,c);
end
% Because some of the displacement vectors might now be outside of the
% reduced zone, bring them back. Furthermore, the class generates
% displacement vectors in the full hexagon, not half hexagon.
edisps = [xdisps(:),ydisps(:)];
[ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( edisps );
xdisps = reduced_zone_disps(:,1);
ydisps = reduced_zone_disps(:,2);
xdisps = reshape(xdisps,scan_dimensions);
ydisps = reshape(ydisps,scan_dimensions);
displacement_field = cat(3,xdisps,ydisps);

if plot_sample
    y = FourDSTEM_Analysis_Engine;
    y.makeOutsideCustomDisplacementColorPlot(displacement_field(:,:,1),displacement_field(:,:,2));
    title(sprintf('%.2f moire angle, [%.2f, %d] reconstruction, no beamwidth averaging',moire_angle,recon_struct.recon_angle,recon_struct.recon_distance));
end

end

