% reconstruction_function_testing_script.m
%
% Nathanael Kazmierczak, 06/08/2020

if false
    moire_angle_deg = 1.1;
    recon_angle = 0;
    recon_distance = 0;
    use_old_reconfun = -1;
    tblg = TwistedBilayerGrapheneAugmented(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun);
    
    xy_lat = tblg.getLayer1();
    cellDimXY = tblg.cellDimXY;
    
    AA_angle = 1.1;
    AB_angle = 1.5;
    AA_distance = 45;
    AB_buffer_distance = 10;
    AB_smooth_distance = 5;
    
    rotated_lattice = reconstructionFunctionNPK( AA_angle, AA_distance, AB_angle, AB_buffer_distance, AB_smooth_distance, xy_lat, cellDimXY );
    figure; scatter(rotated_lattice(:,1),rotated_lattice(:,2)); axis equal; set(gca,'yDir','normal');
end

if true
    moire_angle_deg = 1.1;
    recon_struct.type = 'NPK';
%     recon_struct.AA_angle = 0;
%     recon_struct.AA_distance = 45;
%     recon_struct.AB_angle = -0.6;
%     recon_struct.AB_buffer = 15;
%     recon_struct.AB_smooth = 10;

% Case 1
%     recon_struct.AA_angle = 0;
%     recon_struct.AA_distance = 45;
%     recon_struct.AB_angle = 0;
%     recon_struct.AB_buffer = 15;
%     recon_struct.AB_smooth = 10;

% Case 2
%     recon_struct.AA_angle = 1.0;
%     recon_struct.AA_distance = 45;
%     recon_struct.AB_angle = 0;
%     recon_struct.AB_buffer = 0;
%     recon_struct.AB_smooth = 5;

% recon_struct.recon_angle = 1.0;
%     recon_struct.recon_distance = 45;
%     recon_struct.type = 'Tadmor';

% Case 3
%     recon_struct.AA_angle = 0;
%     recon_struct.AA_distance = 45;
%     recon_struct.AB_angle = -1;
%     recon_struct.AB_buffer = 15;
%     recon_struct.AB_smooth = 10;

% Case 4
    recon_struct.AA_angle = 1;
    recon_struct.AA_distance = 45;
    recon_struct.AB_angle = -1;
    recon_struct.AB_buffer = 15;
    recon_struct.AB_smooth = 10;


    
    tblg = TwistedBilayerGrapheneAugmented(moire_angle_deg,recon_struct);
    tblg.plot();
    tblg.computeDSCField(0);
    tblg.plotDSCField();
    tblg.plotDSCLattice();
    tblg.plotDSCQuiver();
    tblg.computeNoModDSCField();
    remove_mean = 0;
    tblg.computeStrainField(remove_mean);
    type = 'all';
    tblg.plotStrainField([],type);
    figh = tblg.plotDSCField();
    [figh,perc_AA,perc_AB,perc_SP] = tblg.assignPsuedostacking(figh)
end
