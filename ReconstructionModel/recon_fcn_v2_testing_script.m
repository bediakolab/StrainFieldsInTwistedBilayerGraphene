% recon_fcn_v2_testing_script.m

% reconstruction_function_script0p5.m
%
% Nathanael Kazmierczak, 06/08/2020

moire_angle_deg = 0.5;
recon_struct.type = 'NPK_v2';
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
    recon_struct.AA_angle = 0.8;
%     recon_struct.AA_angle = 0;
    recon_struct.AA_distance = 40;
    recon_struct.AB_angle = -0.350;
    recon_struct.AB_buffer = 10;
    recon_struct.AB_smooth = 20;
    
    % Logic behind new reconstruction function: interpolate angle to what
    % would be parallel to the soliton wall. 
    recon_struct.boundary_rectification.type = 'interpolate';
    recon_struct.boundary_rectification.gamma = 10;
    recon_struct.boundary_rectification.corner_angle_deg = 20;
    recon_struct.plotrotfield = false;

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
% recon_struct.AA_angle = 1;
% recon_struct.AA_distance = 45;
% recon_struct.AB_angle = -1;
% recon_struct.AB_buffer = 15;
% recon_struct.AB_smooth = 10;



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


