% heterostrain_atomic_simulator.m
%
% Eventually load this into the 

heterostrain_struct.type = 'real space';
heterostrain_struct.amount = 2;  % percent
heterostrain_struct.amount = 0;  % percent
% heterostrain_struct.shear_amount %%%
heterostrain_struct.angle = 75;  % no rotation before applying the heterostrain 
heterostrain_struct.usePoisson = false;

recon_struct.type = 'NPK';
recon_struct.AA_angle = 0;
recon_struct.AA_distance = 45;
recon_struct.AB_angle = 0;
recon_struct.AB_buffer = 10;
recon_struct.AB_smooth = 5;
recon_struct = struct;
recon_struct.type = 'none';

moire_angle_deg = 2;
custom_cell_size = [300,150];
tblg = TwistedBilayerGrapheneAugmented(moire_angle_deg,recon_struct,heterostrain_struct,custom_cell_size);
tblg.plot(figure,true);

tblg.computeDSCField(0);
tblg.plotDSCField();
